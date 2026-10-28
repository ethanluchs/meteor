/* meteor_detect.c - forward-scatter meteor detector, 16384-point version
 *
 * Changes from the 4096 build:
 *   - FFT 16384: bins 15.26 Hz (was 61.0), noise bandwidth 6 dB narrower.
 *   - Parabolic peak interpolation -> sub-bin frequency, ~1-2 Hz resolution.
 *     This is what makes aircraft Doppler (tens of Hz) visible at all.
 *   - Pilot tracker: follows slow transmitter drift, frozen during events.
 *   - Full-rate raw logging (15.26 rows/s) so drift can be found post-hoc.
 *   - Non-zero exit on device loss so a wrapper can restart.
 *
 * Build: gcc -O2 -o meteor_detect meteor_detect.c -lrtlsdr -lfftw3f -lm
 * Run:   ./meteor_detect              normal
 *        ./meteor_detect --cal        calibrate and exit
 *        ./meteor_detect --bin 1280   force a carrier bin
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <signal.h>
#include <stdint.h>
#include <rtl-sdr.h>
#include <fftw3.h>

/* ---------------------------------------------------------------- radio -- */
#define CENTER_FREQ     54289909u        /* puts the pilot on a bin centre */
#define CARRIER_OFFSET  19531            /* 54309440.6 - 54289909 */
#define SAMPLE_RATE     250000u
#define FFT_SIZE        16384
#define GAIN_TENTH_DB   496              /* sweep showed max is correct */

/* ------------------------------------------------------------ detection -- */
#define SEARCH_HZ       6000
#define GUARD_BINS      6                /* ~92 Hz either side of the peak */
#define BIN_TOL         2                /* carrier +/- 2 bins = +/-30 Hz */
#define WIN_FRAMES      152              /* ~9.96 s occupancy window */
#define TRIG_Z_FLOOR    8.0
#define TRIG_MARGIN     4.0              /* must beat hottest decoy by this */
#define REL_Z           3.0
#define MIN_TRIG_FRAMES 15               /* ~1 s held above the bar */
#define P0_TAU_S        300.0
#define N_DECOY         16
#define CARRIER_GUARD   24               /* keep decoys clear of the carrier */

/* ------------------------------------------------------------- tracking -- */
#define TRACK_WINDOW    900              /* frames of history (~59 s) */
#define TRACK_SPAN      40               /* search +/- bins (+/-610 Hz) */
#define TRACK_PERIOD    150              /* re-evaluate every ~10 s */

#define CAL_SECONDS     120
#define RAW_DECIMATION  1                /* log every frame: 15.26 rows/s */
#define STATUS_SECONDS  300

#define FRAME_HZ        ((double)SAMPLE_RATE / (double)FFT_SIZE)
#define BIN_HZ          ((double)SAMPLE_RATE / (double)FFT_SIZE)

/* --------------------------------------------------------------- state --- */
static rtlsdr_dev_t  *dev     = NULL;
static fftwf_complex *fft_in  = NULL, *fft_out = NULL;
static fftwf_plan     plan;
static float         *hann    = NULL;
static float         *pwr     = NULL;
static float         *noisebuf= NULL;

static int       win_lo, win_hi, nominal_bin, carrier_bin = -1;
static int       cal_only = 0;
static long long frame_idx = 0;
static time_t    t_wall0;
static FILE     *fev = NULL, *fraw = NULL;
static volatile sig_atomic_t stop_flag = 0;
static int       dev_lost = 0;

static int      *cal_hist = NULL;
static long long cal_frames = 0, cal_need = 0;

typedef struct {
    int      bin;
    uint8_t  buf[WIN_FRAMES];
    int      idx, sum, filled;
    double   p0, last_z, max_z;
} occ_t;

static occ_t carrier_occ;
static occ_t decoy_occ[N_DECOY];
static int   n_decoy_used = 0;

static uint16_t trk_buf[TRACK_WINDOW];
static int      trk_idx = 0, trk_filled = 0;

static int       in_event = 0, trig_run = 0;
static long long ev_start_frame = 0;
static double    ev_peak_z = 0.0, ev_peak_snr = 0.0;
static double    ev_peak_frac = 0.0, ev_decoy_at_peak = 0.0;
static double    ev_bin_first = 0.0, ev_bin_last = 0.0;
static int      *ev_hist = NULL;

static void on_sig(int s) { (void)s; stop_flag = 1; }

static int cmpf(const void *a, const void *b)
{
    float x = *(const float *)a, y = *(const float *)b;
    return (x > y) - (x < y);
}

static double now_sec(void) { return (double)frame_idx / FRAME_HZ; }
static double bin_offset_hz(double b) { return (b - nominal_bin) * BIN_HZ; }

/* ------------------------------------------------------------ occupancy -- */
static void occ_init(occ_t *o, int bin, double p0)
{
    memset(o, 0, sizeof(*o));
    o->bin = bin;
    o->p0  = p0 > 1e-6 ? p0 : 1e-6;
}

static double occ_push(occ_t *o, int peak_bin)
{
    int hit = (abs(peak_bin - o->bin) <= BIN_TOL) ? 1 : 0;

    o->sum -= o->buf[o->idx];
    o->buf[o->idx] = (uint8_t)hit;
    o->sum += hit;
    o->idx = (o->idx + 1) % WIN_FRAMES;
    if (o->filled < WIN_FRAMES) o->filled++;

    if (o->last_z <= REL_Z) {                 /* freeze while elevated */
        double a = 1.0 / (P0_TAU_S * FRAME_HZ);
        o->p0 += a * ((double)hit - o->p0);
        if (o->p0 < 1e-6) o->p0 = 1e-6;
    }
    if (o->filled < WIN_FRAMES) { o->last_z = 0.0; return 0.0; }

    double p = o->p0;
    if (p > 0.999999) p = 0.999999;
    double mu = (double)WIN_FRAMES * p;
    double sd = sqrt((double)WIN_FRAMES * p * (1.0 - p));
    double z  = (sd > 0.0) ? ((double)o->sum - mu) / sd : 0.0;

    o->last_z = z;
    if (z > o->max_z) o->max_z = z;
    return z;
}

static double decoy_z_now(void)
{
    double m = -1e9;
    for (int i = 0; i < n_decoy_used; i++)
        if (decoy_occ[i].last_z > m) m = decoy_occ[i].last_z;
    return (n_decoy_used > 0) ? m : 0.0;
}

static double decoy_z_max(void)
{
    double m = 0.0;
    for (int i = 0; i < n_decoy_used; i++)
        if (decoy_occ[i].max_z > m) m = decoy_occ[i].max_z;
    return m;
}

static void place_decoys(int carrier, double p0_init)
{
    int lo = win_lo + BIN_TOL, hi = win_hi - BIN_TOL, span = hi - lo;
    n_decoy_used = 0;
    for (int i = 0; i < N_DECOY + 8 && n_decoy_used < N_DECOY; i++) {
        int b = lo + (int)((double)span * (i + 1) / (double)(N_DECOY + 9));
        if (abs(b - carrier) <= CARRIER_GUARD) continue;
        occ_init(&decoy_occ[n_decoy_used++], b, p0_init);
    }
}

/* -------------------------------------------------------- pilot tracker -- */
static void track_push(int peak)
{
    trk_buf[trk_idx] = (uint16_t)peak;
    trk_idx = (trk_idx + 1) % TRACK_WINDOW;
    if (trk_filled < TRACK_WINDOW) trk_filled++;
}

static void track_update(void)
{
    if (in_event || trig_run > 0) return;          /* never chase an echo */
    if (trk_filled < TRACK_WINDOW / 2) return;

    int lo = carrier_bin - TRACK_SPAN, hi = carrier_bin + TRACK_SPAN;
    if (lo < win_lo) lo = win_lo;
    if (hi > win_hi) hi = win_hi;
    int n = hi - lo + 1;
    static int h[2 * TRACK_SPAN + 1];
    memset(h, 0, sizeof(int) * (size_t)n);

    for (int i = 0; i < trk_filled; i++) {
        int b = trk_buf[i];
        if (b >= lo && b <= hi) h[b - lo]++;
    }
    int best = carrier_bin, bestc = -1;
    for (int k = 0; k < n; k++) {
        int c = 0;
        for (int d = -BIN_TOL; d <= BIN_TOL; d++)
            if (k + d >= 0 && k + d < n) c += h[k + d];
        if (c > bestc) { bestc = c; best = lo + k; }
    }
    if (best != carrier_bin) {
        fprintf(stderr, "[trk] %d -> %d (%+.1f Hz, %d/%d frames)\n",
                carrier_bin, best, bin_offset_hz(best), bestc, trk_filled);
        fflush(stderr);
        carrier_bin = best;
        carrier_occ.bin = best;
    }
}

/* ---------------------------------------------------------- calibration -- */
static void finish_calibration(void)
{
    int best = win_lo; long long bestc = -1;
    for (int k = win_lo; k <= win_hi; k++)
        if (cal_hist[k] > bestc) { bestc = cal_hist[k]; best = k; }

    int    nbins   = win_hi - win_lo + 1;
    double share   = (double)bestc / (double)cal_frames;
    double uniform = 1.0 / (double)nbins;

    long long tol = 0;
    for (int k = best - BIN_TOL; k <= best + BIN_TOL; k++)
        if (k >= win_lo && k <= win_hi) tol += cal_hist[k];
    double p0 = (double)tol / (double)cal_frames;

    carrier_bin = best;

    fprintf(stderr, "[cal] %lld frames over %.1f s\n",
            cal_frames, (double)cal_frames / FRAME_HZ);
    fprintf(stderr, "[cal] top 8 bins:\n");
    for (int n = 0; n < 8; n++) {
        int tb = win_lo; long long tc = -1;
        for (int k = win_lo; k <= win_hi; k++)
            if (cal_hist[k] > tc) { tc = cal_hist[k]; tb = k; }
        if (tc <= 0) break;
        fprintf(stderr, "[cal]   bin %5d  %+9.1f Hz  %.4f  (%.1fx uniform)\n",
                tb, bin_offset_hz(tb), (double)tc / (double)cal_frames,
                (double)tc / (double)cal_frames / uniform);
        cal_hist[tb] = -1;
    }
    cal_hist[best] = (int)bestc;

    fprintf(stderr, "[cal] carrier bin %d, %+.1f Hz from nominal %d\n",
            best, bin_offset_hz(best), nominal_bin);
    fprintf(stderr, "[cal] share %.4f vs uniform %.4f = %.1fx   p0(+/-%d) = %.4f\n",
            share, uniform, share / uniform, BIN_TOL, p0);
    if (share / uniform < 2.0)
        fprintf(stderr, "[cal] WARNING: not beating a random bin.\n");
    if (cal_only) { stop_flag = 1; return; }

    occ_init(&carrier_occ, carrier_bin, p0);
    place_decoys(carrier_bin, uniform * (2 * BIN_TOL + 1));
    fprintf(stderr, "[cal] %d decoys. armed: window %d frames (%.2f s), "
                    "bar z > max(%.1f, decoy+%.1f)\n",
            n_decoy_used, WIN_FRAMES, WIN_FRAMES / FRAME_HZ,
            TRIG_Z_FLOOR, TRIG_MARGIN);
    fflush(stderr);
}

/* --------------------------------------------------------------- frame --- */
static void process_frame(const uint8_t *iq)
{
    for (int i = 0; i < FFT_SIZE; i++) {
        float re = ((float)iq[2 * i]     - 127.4f) / 128.0f;
        float im = ((float)iq[2 * i + 1] - 127.4f) / 128.0f;
        fft_in[i][0] = re * hann[i];
        fft_in[i][1] = im * hann[i];
    }
    fftwf_execute(plan);

    int   peak = win_lo;
    float peak_p = -1.0f;
    for (int k = win_lo; k <= win_hi; k++) {
        float re = fft_out[k][0], im = fft_out[k][1];
        float p  = re * re + im * im;
        pwr[k] = p;
        if (p > peak_p) { peak_p = p; peak = k; }
    }

    /* parabolic interpolation in the log domain -> sub-bin frequency */
    double dsub = 0.0;
    if (peak > win_lo && peak < win_hi) {
        double a = 10.0 * log10((double)pwr[peak - 1] + 1e-30);
        double b = 10.0 * log10((double)pwr[peak]     + 1e-30);
        double c = 10.0 * log10((double)pwr[peak + 1] + 1e-30);
        double den = a - 2.0 * b + c;
        if (fabs(den) > 1e-12) dsub = 0.5 * (a - c) / den;
        if (dsub >  0.5) dsub =  0.5;
        if (dsub < -0.5) dsub = -0.5;
    }
    double fbin = (double)peak + dsub;

    int n = 0;
    for (int k = win_lo; k <= win_hi; k++)
        if (abs(k - peak) > GUARD_BINS) noisebuf[n++] = pwr[k];
    qsort(noisebuf, (size_t)n, sizeof(float), cmpf);
    float noise = noisebuf[n / 2];
    if (noise < 1e-20f) noise = 1e-20f;
    double snr_db = 10.0 * log10((double)peak_p / (double)noise);

    frame_idx++;

    if (carrier_bin < 0) {
        cal_hist[peak]++;
        if (++cal_frames >= cal_need) finish_calibration();
        return;
    }
    if (cal_only) return;

    track_push(peak);
    if ((frame_idx % TRACK_PERIOD) == 0) track_update();

    double z = occ_push(&carrier_occ, peak);
    for (int i = 0; i < n_decoy_used; i++) occ_push(&decoy_occ[i], peak);

    double dz   = decoy_z_now();
    double frac = (double)carrier_occ.sum / (double)WIN_FRAMES;
    double bar  = (TRIG_Z_FLOOR > dz + TRIG_MARGIN) ? TRIG_Z_FLOOR
                                                    : dz + TRIG_MARGIN;

    if (fraw && (frame_idx % RAW_DECIMATION) == 0) {
        fprintf(fraw, "%lld,%.4f,%ld,%.2f,%d,%.3f,%.2f,%d,%.4f,%.3f,%.5f,%.3f,%d\n",
                frame_idx, now_sec(), (long)(t_wall0 + (time_t)now_sec()),
                snr_db, peak, fbin, bin_offset_hz(fbin),
                carrier_occ.sum, frac, z, carrier_occ.p0, dz, carrier_bin);
    }

    if (!in_event) {
        if (z > bar) {
            if (trig_run == 0) {
                ev_start_frame = frame_idx;
                ev_peak_z = ev_peak_snr = ev_peak_frac = ev_decoy_at_peak = 0.0;
                ev_bin_first = fbin;
                memset(ev_hist, 0, sizeof(int) * FFT_SIZE);
            }
            trig_run++;
            if (trig_run == MIN_TRIG_FRAMES) {
                in_event = 1;
                fprintf(stderr, "[evt] open  t=%.2fs z=%.2f decoy=%.2f frac=%.3f\n",
                        now_sec(), z, dz, frac);
                fflush(stderr);
            }
        } else trig_run = 0;
    }

    if (trig_run > 0 || in_event) {
        if (z > ev_peak_z)      { ev_peak_z = z; ev_decoy_at_peak = dz; }
        if (snr_db > ev_peak_snr) ev_peak_snr  = snr_db;
        if (frac > ev_peak_frac)  ev_peak_frac = frac;
        ev_bin_last = fbin;
        ev_hist[peak]++;
    }

    if (in_event && z < REL_Z) {
        int dom = carrier_bin, domc = -1;
        for (int k = win_lo; k <= win_hi; k++)
            if (ev_hist[k] > domc) { domc = ev_hist[k]; dom = k; }

        double span = (double)(frame_idx - ev_start_frame) / FRAME_HZ;
        double est  = span - (double)WIN_FRAMES / FRAME_HZ;
        if (est < 0.0) est = 0.0;

        fprintf(fev, "%ld,%lld,%.3f,%.2f,%.2f,%.2f,%.2f,%.4f,%.4f,%.2f,%d,%.1f,%.1f,%.2f\n",
                (long)(t_wall0 + (time_t)((double)ev_start_frame / FRAME_HZ)),
                ev_start_frame, (double)ev_start_frame / FRAME_HZ,
                span, est, ev_peak_z, ev_decoy_at_peak, ev_peak_frac,
                carrier_occ.p0, ev_peak_snr, dom, bin_offset_hz(dom),
                (ev_bin_last - ev_bin_first) * BIN_HZ, decoy_z_max());
        fflush(fev);

        fprintf(stderr, "[evt] close span=%.2fs est=%.2fs z=%.2f (decoy %.2f) "
                        "frac=%.3f snr=%.2f bin=%d drift=%+.1f Hz\n",
                span, est, ev_peak_z, ev_decoy_at_peak, ev_peak_frac,
                ev_peak_snr, dom, (ev_bin_last - ev_bin_first) * BIN_HZ);
        fflush(stderr);
        in_event = 0; trig_run = 0;
    }

    if ((frame_idx % (long long)(STATUS_SECONDS * FRAME_HZ)) == 0) {
        fprintf(stderr, "[sts] t=%.0fs bin=%d p0=%.4f frac=%.3f z=%.2f | "
                        "carrier_max_z=%.2f decoy_max_z=%.2f\n",
                now_sec(), carrier_bin, carrier_occ.p0, frac, z,
                carrier_occ.max_z, decoy_z_max());
        fflush(stderr);
        if (fraw) fflush(fraw);
    }
}

/* ------------------------------------------------------------- capture --- */
static uint8_t carry[FFT_SIZE * 2];
static int     carry_len = 0;

static void rtl_cb(unsigned char *buf, uint32_t len, void *ctx)
{
    (void)ctx;
    const uint32_t FB = FFT_SIZE * 2;
    uint32_t off = 0;

    if (len == 0) { dev_lost = 1; stop_flag = 1; }

    if (carry_len > 0) {
        uint32_t need = FB - (uint32_t)carry_len;
        if (len < need) { memcpy(carry + carry_len, buf, len);
                          carry_len += (int)len; return; }
        memcpy(carry + carry_len, buf, need);
        process_frame(carry);
        carry_len = 0; off = need;
    }
    while (off + FB <= len) { process_frame(buf + off); off += FB; }
    carry_len = (int)(len - off);
    if (carry_len > 0) memcpy(carry, buf + off, (size_t)carry_len);

    if (stop_flag) rtlsdr_cancel_async(dev);
}

/* ---------------------------------------------------------------- main --- */
int main(int argc, char **argv)
{
    signal(SIGINT, on_sig);
    signal(SIGTERM, on_sig);

    int forced_bin = -1;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--cal")) cal_only = 1;
        else if (!strcmp(argv[i], "--bin") && i + 1 < argc)
            forced_bin = atoi(argv[++i]);
    }

    hann     = malloc(sizeof(float) * FFT_SIZE);
    pwr      = malloc(sizeof(float) * FFT_SIZE);
    noisebuf = malloc(sizeof(float) * FFT_SIZE);
    cal_hist = calloc(FFT_SIZE, sizeof(int));
    ev_hist  = calloc(FFT_SIZE, sizeof(int));
    if (!hann || !pwr || !noisebuf || !cal_hist || !ev_hist) return 1;

    nominal_bin = (int)lround((double)CARRIER_OFFSET / BIN_HZ);
    int half    = (int)((double)SEARCH_HZ / BIN_HZ);
    win_lo = nominal_bin - half;  if (win_lo < 1) win_lo = 1;
    win_hi = nominal_bin + half;  if (win_hi > FFT_SIZE - 2) win_hi = FFT_SIZE - 2;
    cal_need = (long long)(CAL_SECONDS * FRAME_HZ);

    if (rtlsdr_open(&dev, 0) < 0) { fprintf(stderr, "no rtlsdr device\n"); return 2; }
    rtlsdr_set_sample_rate(dev, SAMPLE_RATE);
    if (rtlsdr_set_center_freq(dev, CENTER_FREQ) < 0)
        fprintf(stderr, "[rf ] WARNING: set_center_freq reported an error\n");
    rtlsdr_set_agc_mode(dev, 0);
    rtlsdr_set_tuner_gain_mode(dev, 1);
    rtlsdr_set_tuner_gain(dev, GAIN_TENTH_DB);
    fprintf(stderr, "[rf ] %.4f MHz, %u sps, gain %.1f dB, FFT %d\n",
            CENTER_FREQ / 1e6, SAMPLE_RATE, GAIN_TENTH_DB / 10.0, FFT_SIZE);
    rtlsdr_reset_buffer(dev);

    fft_in  = fftwf_alloc_complex(FFT_SIZE);
    fft_out = fftwf_alloc_complex(FFT_SIZE);
    fftwf_import_wisdom_from_filename("fftw.wisdom");
    plan = fftwf_plan_dft_1d(FFT_SIZE, fft_in, fft_out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftwf_export_wisdom_to_filename("fftw.wisdom");
    for (int i = 0; i < FFT_SIZE; i++)
        hann[i] = 0.5f * (1.0f - cosf(2.0f * (float)M_PI * i / (FFT_SIZE - 1)));

    t_wall0 = time(NULL);

    if (forced_bin > 0) {
        carrier_bin = forced_bin;
        occ_init(&carrier_occ, carrier_bin, 0.035);
        place_decoys(carrier_bin, 0.015);
        fprintf(stderr, "[cal] forced to bin %d (%+.1f Hz)\n",
                carrier_bin, bin_offset_hz(carrier_bin));
    }

    if (!cal_only) {
        fev = fopen("events.csv", "a");
        if (fev && ftell(fev) == 0)
            fprintf(fev, "unix_time,start_frame,start_sec,span_s,est_true_dur_s,"
                         "peak_z,decoy_z_at_peak,peak_frac,base_frac,peak_snr_db,"
                         "dominant_bin,offset_hz,drift_hz,decoy_max_z\n");
        fraw = fopen("raw_snr_log.csv", "a");
        if (fraw && ftell(fraw) == 0)
            fprintf(fraw, "frame_idx,t_sec,unix_time,snr_db,peak_bin,frac_bin,"
                          "offset_hz,win_hits,win_frac,z,p0,decoy_z,carrier_bin\n");
    }

    fprintf(stderr, "[cfg] nominal bin %d, search %d..%d (%d bins), "
                    "%.4f frames/s, %.4f Hz/bin\n",
            nominal_bin, win_lo, win_hi, win_hi - win_lo + 1, FRAME_HZ, BIN_HZ);
    if (carrier_bin < 0)
        fprintf(stderr, "[cal] calibrating %d s%s...\n",
                CAL_SECONDS, cal_only ? " (cal-only)" : "");

    rtlsdr_reset_buffer(dev);
    rtlsdr_reset_buffer(dev);
    rtlsdr_read_async(dev, rtl_cb, NULL, 0, 131072);

    if (!cal_only)
        fprintf(stderr, "[end] %lld frames, %.3f h, carrier_max_z %.2f, "
                        "decoy_max_z %.2f, device_lost=%d\n",
                frame_idx, now_sec() / 3600.0,
                carrier_occ.max_z, decoy_z_max(), dev_lost);

    if (fev)  fclose(fev);
    if (fraw) fclose(fraw);
    fftwf_destroy_plan(plan);
    fftwf_free(fft_in);
    fftwf_free(fft_out);
    rtlsdr_close(dev);
    return (dev_lost && !stop_flag) ? 3 : 0;
}
