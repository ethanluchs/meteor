#include <math.h>
#include <rtl-sdr.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define SAMPLE_RATE 250000
#define CENTER_FREQ 88500000
#define BUF_BYTES 16384

/* detector tuning knobs */
static const float THRESH_DB = 6.0f;
static const float ALPHA = 0.001f;
static const int MIN_BLOCKS = 5;
static const int WARMUP_BLOCKS = 200;

static int check(const char *what, int r) {
  if (r < 0)
    fprintf(stderr, "%s failed: %d\n", what, r);
  return r;
}

int main(void) {
  rtlsdr_dev_t *dev = NULL;

  if (rtlsdr_open(&dev, 0) != 0) {
    fprintf(stderr, "Failed to open device\n");
    return 1;
  }

  check("set_sample_rate", rtlsdr_set_sample_rate(dev, SAMPLE_RATE));
  check("set_center_freq", rtlsdr_set_center_freq(dev, CENTER_FREQ));
  check("set_tuner_gain_mode", rtlsdr_set_tuner_gain_mode(dev, 1));
  check("set_tuner_bandwidth", rtlsdr_set_tuner_bandwidth(dev, 200000));

  rtlsdr_set_agc_mode(dev, 0);

  /* confirm the hardware took what you asked for */
  fprintf(stderr, "actual rate=%u Hz  freq=%u Hz\n",
          rtlsdr_get_sample_rate(dev), rtlsdr_get_center_freq(dev));

  check("reset_buffer", rtlsdr_reset_buffer(dev));

  uint8_t buf[BUF_BYTES];
  int n_read = 0;
  time_t start = time(NULL);

  float baseline_db = 0.0f;
  int baseline_init = 0;
  int hot_blocks = 0;
  long block_count = 0;

  while (time(NULL) - start < 7200) {
    if (rtlsdr_read_sync(dev, buf, sizeof(buf), &n_read) != 0) {
      fprintf(stderr, "read_sync error, aborting\n");
      break;
    }
    if (n_read <= 0)
      continue;

    /* mean power across the whole block, not per-sample magnitude */
    double sum_pow = 0.0;
    long pairs = 0;
    for (int i = 0; i + 1 < n_read; i += 2) {
      float I = (float)buf[i] - 127.5f;
      float Q = (float)buf[i + 1] - 127.5f;
      sum_pow += (double)(I * I + Q * Q);
      pairs++;
    }
    if (pairs == 0)
      continue;

    float mean_pow = (float)(sum_pow / (double)pairs);
    float power_db = 10.0f * log10f(mean_pow + 1e-9f);

    block_count++;
    if (!baseline_init) {
      baseline_db = power_db;
      baseline_init = 1;
      continue;
    }

    int hot =
        (block_count > WARMUP_BLOCKS) && (power_db > baseline_db + THRESH_DB);

    if (hot) {
      hot_blocks++;
      if (hot_blocks == MIN_BLOCKS) {
        printf("[%ld] Detection: %.2f dB (baseline %.2f dB, +%.2f dB)\n",
               (long)time(NULL), power_db, baseline_db, power_db - baseline_db);
        fflush(stdout);
      }
    } else {
      hot_blocks = 0;
      /* only let the baseline move while nothing is happening */
      baseline_db = ALPHA * power_db + (1.0f - ALPHA) * baseline_db;
    }

    /* heartbeat so you can see real numbers instead of guessing */
    if (block_count % 250 == 0) {
      fprintf(stderr, "block %ld  power=%.2f dB  baseline=%.2f dB\n",
              block_count, power_db, baseline_db);
    }
  }

  rtlsdr_close(dev);
  return 0;
}
