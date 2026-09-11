#include <fftw3.h>
#include <math.h>
#include <rtl-sdr.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define SAMPLE_RATE 250000

// tuned 20 kHz low so the carrier of interest lands clear of the DC spike
// that sits at the center frequency
#define CENTER_FREQ 88480000
#define CARRIER_OFFSET 20000

#define TUNER_BW 200000
#define TUNER_GAIN 200

// 4096 point FFT gives 61 Hz bins and a 16.4 ms frame at 250 ksps
#define FFT_N 4096
#define BUF_BYTES (FFT_N * 2)

// how far either side of the expected carrier to look, covering station
// frequency error and the Doppler slide a trail produces
#define SEARCH_HZ 3000

// how far above the noise floor the peak bin has to sit to count
static const float THRESH_DB = 10.0f;

// consecutive frames required before calling it a detection
static const int MIN_FRAMES = 2;

static const int WARMUP_FRAMES = 100;
static const long RUN_SECONDS = 7200;
static const long HEARTBEAT_FRAMES = 60;

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

  fprintf(stderr, "actual rate=%u Hz  freq=%u Hz  gain=%.1f dB\n",
          rtlsdr_get_sample_rate(dev), rtlsdr_get_center_freq(dev),
          rtlsdr_get_tuner_gain(dev) / 10.0);

  float bin_hz = (float)SAMPLE_RATE / FFT_N;
  int carrier_bin = (int)lrintf(CARRIER_OFFSET / bin_hz);
  int half_window = (int)lrintf(SEARCH_HZ / bin_hz);
  int lo_bin = carrier_bin - half_window;
  int hi_bin = carrier_bin + half_window;

  fprintf(stderr, "bin=%.1f Hz  carrier bin=%d  search bins %d..%d\n", bin_hz,
          carrier_bin, lo_bin, hi_bin);

  fftwf_complex *in = fftwf_alloc_complex(FFT_N);
  fftwf_complex *out = fftwf_alloc_complex(FFT_N);
  fftwf_plan plan =
      fftwf_plan_dft_1d(FFT_N, in, out, FFTW_FORWARD, FFTW_MEASURE);

  // Hann window keeps a strong carrier from smearing across neighbouring bins
  float *window = malloc(FFT_N * sizeof(float));
  for (int i = 0; i < FFT_N; i++)
    window[i] = 0.5f * (1.0f - cosf(2.0f * (float)M_PI * i / (FFT_N - 1)));

  float *power = malloc(FFT_N * sizeof(float));
  float *scratch = malloc(FFT_N * sizeof(float));

  check("reset_buffer", rtlsdr_reset_buffer(dev));

  uint8_t buf[BUF_BYTES];
  int n_read = 0;
  time_t start = time(NULL);
  int hot_frames = 0;
  long frame_count = 0;

  while (time(NULL) - start < RUN_SECONDS) {
    if (rtlsdr_read_sync(dev, buf, sizeof(buf), &n_read) != 0) {
      fprintf(stderr, "read_sync error, aborting\n");
      break;
    }
    if (n_read < BUF_BYTES)
      continue;

    for (int i = 0; i < FFT_N; i++) {
      in[i][0] = ((float)buf[2 * i] - 127.5f) * window[i];
      in[i][1] = ((float)buf[2 * i + 1] - 127.5f) * window[i];
    }

    fftwf_execute(plan);

    for (int k = 0; k < FFT_N; k++)
      power[k] = out[k][0] * out[k][0] + out[k][1] * out[k][1];

    // the median across the spectrum is a robust noise floor: a narrow
    // carrier occupies too few bins to shift it
    memcpy(scratch, power, FFT_N * sizeof(float));
    qsort(scratch, FFT_N, sizeof(float), cmp_float);
    float noise = scratch[FFT_N / 2];

    float peak = 0.0f;
    int peak_bin = lo_bin;
    for (int k = lo_bin; k <= hi_bin; k++) {
      if (power[k] > peak) {
        peak = power[k];
        peak_bin = k;
      }
    }

    float snr_db = 10.0f * log10f((peak + 1e-12f) / (noise + 1e-12f));
    float doppler_hz = (peak_bin - carrier_bin) * bin_hz;

    frame_count++;

    if (frame_count > WARMUP_FRAMES && snr_db > THRESH_DB) {
      hot_frames++;
      if (hot_frames == MIN_FRAMES) {
        printf("[%ld] Detection: SNR=%.1f dB  doppler=%+.0f Hz\n",
               (long)time(NULL), snr_db, doppler_hz);
        fflush(stdout);
      }
    } else {
      hot_frames = 0;
    }

    if (frame_count % HEARTBEAT_FRAMES == 0) {
      fprintf(stderr, "frame %ld  peak SNR=%.1f dB  doppler=%+.0f Hz\n",
              frame_count, snr_db, doppler_hz);
    }
  }

  fftwf_destroy_plan(plan);
  fftwf_free(in);
  fftwf_free(out);
  free(window);
  free(power);
  free(scratch);
  rtlsdr_close(dev);
  return 0;
}
