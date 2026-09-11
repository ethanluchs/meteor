#include <fftw3.h>
#include <math.h>
#include <rtl-sdr.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define SAMPLE_RATE 250000
#define CENTER_FREQ 88280000
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

// roughly ten seconds of frames per summary row
#define FRAMES_PER_BUCKET 600

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

  fprintf(stderr, "rate=%u freq=%u gain=%.1f dB\n", rtlsdr_get_sample_rate(dev),
          rtlsdr_get_center_freq(dev), rtlsdr_get_tuner_gain(dev) / 10.0);
  fprintf(stderr, "bin=%.2f Hz  window bins %d..%d (%d bins, +/-%d Hz)\n",
          bin_hz, lo_bin, hi_bin, n_window, SEARCH_HZ);

  fftwf_complex *in = fftwf_alloc_complex(FFT_N);
  fftwf_complex *out = fftwf_alloc_complex(FFT_N);
  fftwf_plan plan =
      fftwf_plan_dft_1d(FFT_N, in, out, FFTW_FORWARD, FFTW_MEASURE);

  float *window = malloc(FFT_N * sizeof(float));
  for (int i = 0; i < FFT_N; i++)
    window[i] = 0.5f * (1.0f - cosf(2.0f * (float)M_PI * i / (FFT_N - 1)));

  float *power = malloc(FFT_N * sizeof(float));
  float *scratch = malloc(FFT_N * sizeof(float));

  // per bucket accumulators
  float *snr_samples = malloc(FRAMES_PER_BUCKET * sizeof(float));
  float *noise_samples = malloc(FRAMES_PER_BUCKET * sizeof(float));
  long *bin_hits = malloc((size_t)n_window * sizeof(long));

  check("reset_buffer", rtlsdr_reset_buffer(dev));

  uint8_t buf[BUF_BYTES];
  int n_read = 0;
  time_t start = time(NULL);
  int in_bucket = 0;
  long total_frames = 0;

  memset(bin_hits, 0, (size_t)n_window * sizeof(long));

  printf("unix_time,iso_time,frames,noise_db,snr_median,snr_p95,snr_max,"
         "mode_bin,mode_offset_hz,mode_fraction,clip_pct\n");

  while (time(NULL) - start < RUN_SECONDS) {
    if (rtlsdr_read_sync(dev, buf, sizeof(buf), &n_read) != 0) {
      fprintf(stderr, "read_sync error, aborting\n");
      break;
    }
    if (n_read < BUF_BYTES)
      continue;

    long clipped = 0;
    for (int i = 0; i < FFT_N; i++) {
      uint8_t bi = buf[2 * i];
      uint8_t bq = buf[2 * i + 1];
      if (bi <= 1 || bi >= 254 || bq <= 1 || bq >= 254)
        clipped++;
      in[i][0] = ((float)bi - 127.5f) * window[i];
      in[i][1] = ((float)bq - 127.5f) * window[i];
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
      if (power[k] > peak) {
        peak = power[k];
        peak_bin = k;
      }
    }

    snr_samples[in_bucket] = 10.0f * log10f((peak + 1e-12f) / (noise + 1e-12f));
    noise_samples[in_bucket] = 10.0f * log10f(noise + 1e-12f);
    bin_hits[peak_bin - lo_bin]++;
    in_bucket++;
    total_frames++;

    if (in_bucket < FRAMES_PER_BUCKET)
      continue;

    // find the bin that won most often this bucket
    long best_hits = 0;
    int best_idx = 0;
    for (int k = 0; k < n_window; k++) {
      if (bin_hits[k] > best_hits) {
        best_hits = bin_hits[k];
        best_idx = k;
      }
    }

    float noise_med;
    memcpy(scratch, noise_samples, (size_t)in_bucket * sizeof(float));
    qsort(scratch, (size_t)in_bucket, sizeof(float), cmp_float);
    noise_med = scratch[in_bucket / 2];

    memcpy(scratch, snr_samples, (size_t)in_bucket * sizeof(float));
    qsort(scratch, (size_t)in_bucket, sizeof(float), cmp_float);
    float snr_med = scratch[in_bucket / 2];
    float snr_p95 = scratch[(int)(in_bucket * 0.95f)];
    float snr_max = scratch[in_bucket - 1];

    time_t now = time(NULL);
    struct tm tmv;
    localtime_r(&now, &tmv);
    char iso[32];
    strftime(iso, sizeof(iso), "%Y-%m-%dT%H:%M:%S", &tmv);

    printf("%ld,%s,%d,%.2f,%.2f,%.2f,%.2f,%d,%+.0f,%.4f,%.3f\n", (long)now, iso,
           in_bucket, noise_med, snr_med, snr_p95, snr_max, lo_bin + best_idx,
           (lo_bin + best_idx - carrier_bin) * bin_hz,
           (double)best_hits / in_bucket, 100.0 * (double)clipped / FFT_N);

    in_bucket = 0;
    memset(bin_hits, 0, (size_t)n_window * sizeof(long));
  }

  fprintf(stderr, "done, %ld frames\n", total_frames);

  fftwf_destroy_plan(plan);
  fftwf_free(in);
  fftwf_free(out);
  free(window);
  free(power);
  free(scratch);
  free(snr_samples);
  free(noise_samples);
  free(bin_hits);
  rtlsdr_close(dev);
  return 0;
}
