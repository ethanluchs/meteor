/* meteor_detect.c - forward-scatter meteor detector, occupancy-based
 *
 * Detection statistic is "what fraction of the last WIN_FRAMES frames had
 * their peak at the carrier bin", not per-frame SNR. Occupancy averages down
 * as 1/sqrt(N); per-frame snr_db does not, which is why sigma thresholds on
 * snr_db can never fire on a real event.
 *
 * The trigger is gated on N_DECOY control bins running the identical
 * statistic in the same data. Analytic binomial z is not trustworthy here:
 * the bin-winning process is clumpy, and decoys empirically reach z~12 where
 * independent-Bernoulli theory predicts ~4.3.
 *
 * Build: gcc -O2 -o meteor_detect meteor_detect.c -lrtlsdr -lfftw3f -lm
 * Run:   ./meteor_detect              normal, auto-locks the carrier bin
 *        ./meteor_detect --cal        calibrate and exit (antenna/retune tests)
 *        ./meteor_detect --bin 328    skip calibration, force a bin
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
#define CENTER_FREQ     90680000u
#define CARRIER_OFFSET  20000            /* target = center + offset = 54.310 MHz */
#define SAMPLE_RATE     250000u
#define FFT_SIZE        4096
#define GAIN_TENTH_DB   -1               /* -1 = highest gain the tuner offers */

/* ------------------------------------------------------------ detection -- */
#define SEARCH_HZ       6000             /* +/- around nominal carrier */
#define GUARD_BINS      3                /* excluded around peak for noise median */
#define BIN_TOL         1                /* carrier +/- this many bins is a "hit" */
#define WIN_FRAMES      600              /* ~9.83 s occupancy window */
#define TRIG_Z_FLOOR    8.0              /* absolute floor; raise after night 1 */
#define TRIG_MARGIN     4.0              /* must beat the hottest decoy by this */
#define REL_Z           3.0
#define MIN_TRIG_FRAMES 60               /* ~1 s held above trigger to open */
#define P0_TAU_S        300.0            /* baseline occupancy EMA time constant */
#define N_DECOY         16
#define CARRIER_GUARD   6                /* keep decoys this far from the carrier */

#define CAL_SECONDS     120
#define RAW_DECIMATION  10
#define STATUS_SECONDS  300

#define FRAME_HZ        ((double)SAMPLE_RATE / (double)FFT_SIZE)
#define BIN_HZ          ((double)SAMPLE_RATE / (double)FFT_SIZE)

/* --------------------------------------------------------------- state --- */
static rtlsdr_dev_t  *dev     = NULL;
static fftwf_complex *fft_in  = NULL, *fft_out = NULL;
static fftwf_plan     plan;
static float          hann[FFT_SIZE];
static float          pwr[FFT_SIZE];
static float          noisebuf[FFT_SIZE];

static int       win_lo, win_hi, nominal_bin, carrier_bin = -1;
static int       cal_only = 0;
static long long frame_idx = 0;
static time_t    t_wall0;
static FILE     *fev = NULL, *fraw = NULL;
static volatile sig_atomic_t stop_flag = 0;

static int       cal_hist[FFT_SIZE];
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

static int       in_event = 0, trig_run = 0;
static long long ev_start_frame = 0;
static double    ev_peak_z = 0.0, ev_peak_snr = 0.0;
static double    ev_peak_frac = 0.0, ev_decoy_at_peak = 0.0;
static int       ev_hist[FFT_SIZE];

static void on_sig(int s) { (void)s; stop_flag = 1; }

static int cmpf(const void *a, const void *b)
{
    float x = *(const float *)a, y = *(const float *)b;
    return (x > y) - (x < y);
}

static double now_sec(void) { return (double)frame_idx / FRAME_HZ; }

static double bin_offset_hz(int b) { return (double)(b - nominal_bin) * BIN_HZ; }

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

    /* freeze the baseline while this bin is elevated, or a real event drags
       its own baseline up and erases itself */
    if (o->last_z <= REL_Z) {
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
    int lo = win_lo + BIN_TOL, hi = win_hi - BIN_TOL;
    int span = hi - lo;
    n_decoy_used = 0;
    for (int i = 0; i < N_DECOY + 6 && n_decoy_used < N_DECOY; i++) {
        int b = lo + (int)((double)span * (i + 1) / (double)(N_DECOY + 7));
        if (abs(b - carrier) <= CARRIER_GUARD) continue;
        occ_init(&decoy_occ[n_decoy_used], b, p0_init);
        n_decoy_used++;
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
        int    tb = win_lo; long long tc = -1;
        for (int k = win_lo; k <= win_hi; k++) {
            int seen = 0;
            for (int m = 0; m < n; m++) (void)m;
            if (cal_hist[k] > tc) { tc = cal_hist[k]; tb = k; (void)seen; }
        }
        fprintf(stderr, "[cal]   bin %4d  %+8.1f Hz  %.4f  (%.1fx uniform)\n",
                tb, bin_offset_hz(tb),
                (double)tc / (double)cal_frames,
                (double)tc / (double)cal_frames / uniform);
        cal_hist[tb] = -1;                 /* consumed for this listing */
    }
    cal_hist[best] = (int)bestc;           /* restore the winner */

    fprintf(stderr, "[cal] carrier bin %d, offset %+.1f Hz from nominal bin %d\n",
            best, bin_offset_hz(best), nominal_bin);
    fprintf(stderr, "[cal] share %.4f vs uniform %.4f = %.1fx   "
                    "p0(+/-%d) = %.4f\n",
            share, uniform, share / uniform, BIN_TOL, p0);

    if (share / uniform < 2.0)
        fprintf(stderr, "[cal] WARNING: not beating a random bin. Check tuning "
                        "and antenna before trusting anything downstream.\n");
    if (fabs(bin_offset_hz(best)) > 3.0 * BIN_HZ)
        fprintf(stderr, "[cal] WARNING: winner is %.0f Hz off nominal. That is "
                        "%.1f ppm at the target frequency, which is far beyond a "
                        "TCXO's spec. Suspect a spur; rerun with the antenna "
                        "disconnected and compare.\n",
                bin_offset_hz(best),
                fabs(bin_offset_hz(best)) /
                ((CENTER_FREQ + CARRIER_OFFSET) / 1e6));

    if (cal_only) { stop_flag = 1; return; }

    occ_init(&carrier_occ, carrier_bin, p0);
    place_decoys(carrier_bin, uniform * (2 * BIN_TOL + 1));
    fprintf(stderr, "[cal] %d decoy bins placed. armed: window %d frames "
                    "(%.2f s), trigger z > max(%.1f, decoy_z + %.1f)\n",
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

    double z = occ_push(&carrier_occ, peak);
    for (int i = 0; i < n_decoy_used; i++) occ_push(&decoy_occ[i], peak);

    double dz   = decoy_z_now();
    double frac = (double)carrier_occ.sum / (double)WIN_FRAMES;
    double bar  = (TRIG_Z_FLOOR > dz + TRIG_MARGIN) ? TRIG_Z_FLOOR
                                                    : dz + TRIG_MARGIN;

    if (fraw && (frame_idx % RAW_DECIMATION) == 0) {
        fprintf(fraw, "%lld,%.4f,%ld,%.2f,%d,%.1f,%d,%.4f,%.3f,%.5f,%.3f\n",
                frame_idx, now_sec(), (long)(t_wall0 + (time_t)now_sec()),
                snr_db, peak, bin_offset_hz(peak),
                carrier_occ.sum, frac, z, carrier_occ.p0, dz);
    }

    if (!in_event) {
        if (z > bar) {
            if (trig_run == 0) {
                ev_start_frame = frame_idx;
                ev_peak_z = ev_peak_snr = ev_peak_frac = ev_decoy_at_peak = 0.0;
                memset(ev_hist, 0, sizeof(ev_hist));
            }
            trig_run++;
            if (trig_run == MIN_TRIG_FRAMES) {
                in_event = 1;
                fprintf(stderr, "[evt] open  t=%.2fs z=%.2f decoy=%.2f frac=%.3f\n",
                        now_sec(), z, dz, frac);
                fflush(stderr);
            }
        } else {
            trig_run = 0;
        }
    }

    if (trig_run > 0 || in_event) {
        if (z > ev_peak_z)        { ev_peak_z = z; ev_decoy_at_peak = dz; }
        if (snr_db > ev_peak_snr)   ev_peak_snr  = snr_db;
        if (frac > ev_peak_frac)    ev_peak_frac = frac;
        ev_hist[peak]++;
    }

    if (in_event && z < REL_Z) {
        int dom = carrier_bin, domc = -1;
        for (int k = win_lo; k <= win_hi; k++)
            if (ev_hist[k] > domc) { domc = ev_hist[k]; dom = k; }

        double span = (double)(frame_idx - ev_start_frame) / FRAME_HZ;
        double est  = span - (double)WIN_FRAMES / FRAME_HZ;
        if (est < 0.0) est = 0.0;

        fprintf(fev, "%ld,%lld,%.3f,%.2f,%.2f,%.2f,%.2f,%.4f,%.4f,%.2f,%d,%.1f,%.2f\n",
                (long)(t_wall0 + (time_t)((double)ev_start_frame / FRAME_HZ)),
                ev_start_frame, (double)ev_start_frame / FRAME_HZ,
                span, est, ev_peak_z, ev_decoy_at_peak, ev_peak_frac,
                carrier_occ.p0, ev_peak_snr, dom, bin_offset_hz(dom),
                decoy_z_max());
        fflush(fev);

        fprintf(stderr, "[evt] close span=%.2fs est=%.2fs z=%.2f (decoy %.2f) "
                        "frac=%.3f base=%.4f snr=%.2f bin=%d\n",
                span, est, ev_peak_z, ev_decoy_at_peak, ev_peak_frac,
                carrier_occ.p0, ev_peak_snr, dom);
        fflush(stderr);

        in_event = 0; trig_run = 0;
    }

    if ((frame_idx % (long long)(STATUS_SECONDS * FRAME_HZ)) == 0) {
        fprintf(stderr, "[sts] t=%.0fs p0=%.4f frac=%.3f z=%.2f | "
                        "carrier_max_z=%.2f decoy_max_z=%.2f\n",
                now_sec(), carrier_occ.p0, frac, z,
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

    nominal_bin = (int)lround((double)CARRIER_OFFSET / BIN_HZ);
    int half    = (int)((double)SEARCH_HZ / BIN_HZ);
    win_lo = nominal_bin - half;  if (win_lo < 1) win_lo = 1;
    win_hi = nominal_bin + half;  if (win_hi > FFT_SIZE - 2) win_hi = FFT_SIZE - 2;
    cal_need = (long long)(CAL_SECONDS * FRAME_HZ);

    if (rtlsdr_open(&dev, 0) < 0) { fprintf(stderr, "no rtlsdr device\n"); return 1; }
    rtlsdr_set_sample_rate(dev, SAMPLE_RATE);
    rtlsdr_set_center_freq(dev, CENTER_FREQ);
    rtlsdr_set_agc_mode(dev, 0);
    rtlsdr_set_tuner_gain_mode(dev, 1);
    {
        int gains[128];
        int ng = rtlsdr_get_tuner_gains(dev, gains);
        int g  = GAIN_TENTH_DB;
        if (g < 0 && ng > 0) g = gains[ng - 1];
        rtlsdr_set_tuner_gain(dev, g);
        fprintf(stderr, "[rf ] %.4f MHz center, %u sps, gain %.1f dB\n",
                CENTER_FREQ / 1e6, SAMPLE_RATE, g / 10.0);
    }
    rtlsdr_reset_buffer(dev);

    fft_in  = fftwf_alloc_complex(FFT_SIZE);
    fft_out = fftwf_alloc_complex(FFT_SIZE);
    plan    = fftwf_plan_dft_1d(FFT_SIZE, fft_in, fft_out, FFTW_FORWARD, FFTW_MEASURE);
    for (int i = 0; i < FFT_SIZE; i++)
        hann[i] = 0.5f * (1.0f - cosf(2.0f * (float)M_PI * i / (FFT_SIZE - 1)));

    t_wall0 = time(NULL);

    if (forced_bin > 0) {
        carrier_bin = forced_bin;
        occ_init(&carrier_occ, carrier_bin, 0.035);
        place_decoys(carrier_bin, 0.015);
        fprintf(stderr, "[cal] skipped, carrier forced to bin %d (%+.1f Hz)\n",
                carrier_bin, bin_offset_hz(carrier_bin));
    }

    if (!cal_only) {
        fev = fopen("events.csv", "a");
        if (fev && ftell(fev) == 0)
            fprintf(fev, "unix_time,start_frame,start_sec,span_s,est_true_dur_s,"
                         "peak_z,decoy_z_at_peak,peak_frac,base_frac,peak_snr_db,"
                         "dominant_bin,offset_hz,decoy_max_z\n");
        fraw = fopen("raw_snr_log.csv", "a");
        if (fraw && ftell(fraw) == 0)
            fprintf(fraw, "frame_idx,t_sec,unix_time,snr_db,peak_bin,"
                          "peak_offset_hz,win_hits,win_frac,z,p0,decoy_z\n");
    }

    fprintf(stderr, "[cfg] nominal bin %d, search %d..%d (%d bins), "
                    "%.3f frames/s, %.3f Hz/bin\n",
            nominal_bin, win_lo, win_hi, win_hi - win_lo + 1, FRAME_HZ, BIN_HZ);
    if (carrier_bin < 0)
        fprintf(stderr, "[cal] calibrating %d s%s...\n",
                CAL_SECONDS, cal_only ? " (cal-only mode)" : "");

    rtlsdr_read_async(dev, rtl_cb, NULL, 0, 131072);

    if (!cal_only)
        fprintf(stderr, "[end] %lld frames, %.2f h, carrier_max_z %.2f, "
                        "decoy_max_z %.2f\n",
                frame_idx, now_sec() / 3600.0,
                carrier_occ.max_z, decoy_z_max());

    if (fev)  fclose(fev);
    if (fraw) fclose(fraw);
    fftwf_destroy_plan(plan);
    fftwf_free(fft_in);
    fftwf_free(fft_out);
    rtlsdr_close(dev);
    return 0;
}
