#include <fftw3.h>
#include <math.h>
#include <rtl-sdr.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define SAMPLE_RATE 250000
#define CENTER_FREQ 54290000
#define CARRIER_OFFSET 20000
#define TUNER_BW 200000
#define TUNER_GAIN 200

#define FFT_N 4096
#define BUF_BYTES (FFT_N * 2)

// wider than the detector's window so we can see where the tone actually
// sits rather than assuming it lands where we expect
#define SEARCH_HZ 6000

// bins either side of the search window used to estimate the noise floor
#define NOISE_SPAN 900

// how many frames before the rolling baseline is trusted
#define WARMUP_FRAMES 2000

// how far above the rolling baseline (in std devs) counts as elevated
#define N_SIGMA 4.0f

// how slowly the rolling baseline adapts -- small = stable, resistant to
// being dragged around by a real event while it's happening
#define BASELINE_ALPHA 0.01f

// minimum sustained duration, in seconds, before an elevated stretch is
// logged as a real event rather than a single noisy frame
#define MIN_EVENT_SECONDS 1.0f

static const long RUN_SECONDS = 21600;

static int check(const char *what, int r) {
  if (r < 0)
    fprintf(stderr, "%s failed: %d\n", what, r);
  return r;
}

static int cmp_float(const void *a, const void *b) {
  float fa = *(const float *)a, fb = *(const float *)b;
  return (fa > fb) - (fa < fb);
}

int main(void) {
  setvbuf(stdout, NULL, _IOLBF, 0);

  rtlsdr_dev_t *dev = NULL;
  if (rtlsdr_open(&dev, 0) != 0) {
    fprintf(stderr, "Failed to open device\n");
    return 1;
  }

  check("set_sample_rate", rtlsdr_set_sample_rate(dev, SAMPLE_RATE));
  check("set_tuner_bandwidth", rtlsdr_set_tuner_bandwidth(dev, TUNER_BW));
  check("set_center_freq", rtlsdr_set_center_freq(dev, CENTER_FREQ));
  check("set_tuner_gain_mode", rtlsdr_set_tuner_gain_mode(dev, 1));
  check("set_tuner_gain", rtlsdr_set_tuner_gain(dev, TUNER_GAIN));
  check("set_agc_mode", rtlsdr_set_agc_mode(dev, 0));

  float bin_hz = (float)SAMPLE_RATE / FFT_N;
  int carrier_bin = (int)lrintf(CARRIER_OFFSET / bin_hz);
  int half_window = (int)lrintf(SEARCH_HZ / bin_hz);
  int lo_bin = carrier_bin - half_window;
  int hi_bin = carrier_bin + half_window;
  int n_window = hi_bin - lo_bin + 1;

  // frames-per-second, used to translate MIN_EVENT_SECONDS into a frame count
  float frames_per_sec = (float)SAMPLE_RATE / FFT_N;
  int min_event_frames = (int)lrintf(MIN_EVENT_SECONDS * frames_per_sec);
  if (min_event_frames < 1)
    min_event_frames = 1;

  fprintf(stderr, "rate=%u freq=%u gain=%.1f dB\n", rtlsdr_get_sample_rate(dev),
          rtlsdr_get_center_freq(dev), rtlsdr_get_tuner_gain(dev) / 10.0);
  fprintf(stderr, "bin=%.2f Hz  window bins %d..%d (%d bins, +/-%d Hz)\n",
          bin_hz, lo_bin, hi_bin, n_window, SEARCH_HZ);
  fprintf(stderr, "min_event_frames=%d (~%.2f sec at %.1f frames/sec)\n",
          min_event_frames, MIN_EVENT_SECONDS, frames_per_sec);

  fftwf_complex *in = fftwf_alloc_complex(FFT_N);
  fftwf_complex *out = fftwf_alloc_complex(FFT_N);
  fftwf_plan plan =
      fftwf_plan_dft_1d(FFT_N, in, out, FFTW_FORWARD, FFTW_MEASURE);

  float *window = malloc(FFT_N * sizeof(float));
  for (int i = 0; i < FFT_N; i++)
    window[i] = 0.5f * (1.0f - cosf(2.0f * (float)M_PI * i / (FFT_N - 1)));

  float *power = malloc(FFT_N * sizeof(float));
  float *scratch = malloc(FFT_N * sizeof(float));

  // histogram of which bin won, accumulated only while an event is open
  long *event_bin_hits = malloc((size_t)n_window * sizeof(long));
  memset(event_bin_hits, 0, (size_t)n_window * sizeof(long));

  check("reset_buffer", rtlsdr_reset_buffer(dev));

  uint8_t buf[BUF_BYTES];
  int n_read = 0;
  time_t start = time(NULL);
  long frame_count = 0;

  // rolling baseline of snr_db
  int baseline_initialized = 0;
  float baseline_mean = 0.0f;
  float baseline_var = 0.0f;

  // active-event tracking
  int in_event = 0;
  long event_frames = 0;
  float event_peak_snr = 0.0f;
  int event_peak_bin = 0;
  time_t event_start_time = 0;

  printf("start_iso,end_iso,duration_s,frames,peak_snr,peak_offset_hz,"
         "dominant_bin,dominant_offset_hz,dominant_fraction\n");

  while (time(NULL) - start < RUN_SECONDS) {
    if (rtlsdr_read_sync(dev, buf, sizeof(buf), &n_read) != 0) {
      fprintf(stderr, "read_sync error, aborting\n");
      break;
    }
    if (n_read < BUF_BYTES)
      continue;

    for (int i = 0; i < FFT_N; i++) {
      float bi = (float)buf[2 * i] - 127.5f;
      float bq = (float)buf[2 * i + 1] - 127.5f;
      in[i][0] = bi * window[i];
      in[i][1] = bq * window[i];
    }

    fftwf_execute(plan);

    for (int k = 0; k < FFT_N; k++)
      power[k] = out[k][0] * out[k][0] + out[k][1] * out[k][1];

    // noise floor from bins flanking the window, so the IF filter rolloff
    // at the band edges cannot drag the estimate down
    int n_noise = 0;
    for (int k = carrier_bin - NOISE_SPAN; k <= carrier_bin + NOISE_SPAN; k++) {
      if (k < 0 || k >= FFT_N)
        continue;
      if (k >= lo_bin && k <= hi_bin)
        continue;
      scratch[n_noise++] = power[k];
    }
    qsort(scratch, (size_t)n_noise, sizeof(float), cmp_float);
    float noise = scratch[n_noise / 2];

    float peak = 0.0f;
    int peak_bin = lo_bin;
    for (int k = lo_bin; k <= hi_bin; k++) {
      // >= rather than > so a flat/tied noise floor doesn't systematically
      // favor lo_bin every time -- ties now go to the last bin scanned
      if (power[k] >= peak) {
        peak = power[k];
        peak_bin = k;
      }
    }

    float snr_db = 10.0f * log10f((peak + 1e-12f) / (noise + 1e-12f));
    frame_count++;

    // update the rolling baseline unconditionally, every frame, regardless
    // of whether we're currently inside a candidate event -- alpha is small
    // enough that a real event barely moves it, but it still lets the
    // baseline track genuine long-term drift in the noise floor
    if (!baseline_initialized) {
      baseline_mean = snr_db;
      baseline_var = 0.0f;
      baseline_initialized = 1;
    } else {
      float diff = snr_db - baseline_mean;
      baseline_mean += BASELINE_ALPHA * diff;
      baseline_var = (1.0f - BASELINE_ALPHA) *
                     (baseline_var + BASELINE_ALPHA * diff * diff);
    }
    float baseline_std = sqrtf(baseline_var);
    float threshold = baseline_mean + N_SIGMA * baseline_std;

    int elevated = (frame_count > WARMUP_FRAMES) && (snr_db > threshold);

    if (elevated) {
      if (!in_event) {
        // starting a brand new candidate event
        in_event = 1;
        event_frames = 0;
        event_peak_snr = 0.0f;
        event_peak_bin = peak_bin;
        event_start_time = time(NULL);
        memset(event_bin_hits, 0, (size_t)n_window * sizeof(long));
      }
      event_frames++;
      event_bin_hits[peak_bin - lo_bin]++;
      if (snr_db > event_peak_snr) {
        event_peak_snr = snr_db;
        event_peak_bin = peak_bin;
      }
    } else if (in_event) {
      // the elevated stretch just ended -- decide whether it was long
      // enough to count as a real event, then log and reset either way
      if (event_frames >= min_event_frames) {
        time_t event_end_time = time(NULL);
        long best_hits = 0;
        int best_idx = 0;
        for (int k = 0; k < n_window; k++) {
          if (event_bin_hits[k] > best_hits) {
            best_hits = event_bin_hits[k];
            best_idx = k;
          }
        }
        int dominant_bin = lo_bin + best_idx;
        float dominant_offset = (dominant_bin - carrier_bin) * bin_hz;
        float duration_s = (float)event_frames / frames_per_sec;
        float dominant_fraction = (float)best_hits / (float)event_frames;
        float peak_offset = (event_peak_bin - carrier_bin) * bin_hz;

        struct tm tm_start, tm_end;
        char iso_start[32], iso_end[32];
        localtime_r(&event_start_time, &tm_start);
        localtime_r(&event_end_time, &tm_end);
        strftime(iso_start, sizeof(iso_start), "%Y-%m-%dT%H:%M:%S", &tm_start);
        strftime(iso_end, sizeof(iso_end), "%Y-%m-%dT%H:%M:%S", &tm_end);

        printf("%s,%s,%.2f,%ld,%.2f,%+.0f,%d,%+.0f,%.4f\n", iso_start, iso_end,
               duration_s, event_frames, event_peak_snr, peak_offset,
               dominant_bin, dominant_offset, dominant_fraction);
      }
      in_event = 0;
      event_frames = 0;
    }
  }

  fprintf(stderr, "done, %ld frames\n", frame_count);

  fftwf_destroy_plan(plan);
  fftwf_free(in);
  fftwf_free(out);
  free(window);
  free(power);
  free(scratch);
  free(event_bin_hits);
  rtlsdr_close(dev);
  return 0;
}
