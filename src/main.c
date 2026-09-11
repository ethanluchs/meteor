#include <math.h>
#include <rtl-sdr.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define SAMPLE_RATE 250000
#define CENTER_FREQ 88500000
#define TUNER_BW 200000

// tenths of a dB, so this is 20.0 dB
#define TUNER_GAIN 200

// 2048 IQ pairs, about 8.2 ms per block at 250 ksps
#define BUF_BYTES 4096

// how far above the baseline a block has to sit to count as a hit
static const float THRESH_DB = 6.0f;

// baseline smoothing, applied once per block rather than per sample
static const float ALPHA = 0.001f;

// consecutive hot blocks required before reporting a detection
static const int MIN_BLOCKS = 5;

// let the baseline settle before arming the detector
static const int WARMUP_BLOCKS = 200;

static const long RUN_SECONDS = 7200;
static const long HEARTBEAT_BLOCKS = 30;

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

  // order matters here: rate and bandwidth before tuning, gain after gain mode
  check("set_sample_rate", rtlsdr_set_sample_rate(dev, SAMPLE_RATE));
  check("set_tuner_bandwidth", rtlsdr_set_tuner_bandwidth(dev, TUNER_BW));
  check("set_center_freq", rtlsdr_set_center_freq(dev, CENTER_FREQ));
  check("set_tuner_gain_mode", rtlsdr_set_tuner_gain_mode(dev, 1));
  check("set_tuner_gain", rtlsdr_set_tuner_gain(dev, TUNER_GAIN));
  check("set_agc_mode", rtlsdr_set_agc_mode(dev, 0));

  // read the settings back so we know the hardware actually took them
  fprintf(stderr, "actual rate=%u Hz  freq=%u Hz  gain=%.1f dB\n",
          rtlsdr_get_sample_rate(dev), rtlsdr_get_center_freq(dev),
          rtlsdr_get_tuner_gain(dev) / 10.0);

  check("reset_buffer", rtlsdr_reset_buffer(dev));

  uint8_t buf[BUF_BYTES];
  int n_read = 0;
  time_t start = time(NULL);

  float baseline_db = 0.0f;
  int baseline_init = 0;
  int hot_blocks = 0;
  long block_count = 0;

  while (time(NULL) - start < RUN_SECONDS) {
    if (rtlsdr_read_sync(dev, buf, sizeof(buf), &n_read) != 0) {
      fprintf(stderr, "read_sync error, aborting\n");
      break;
    }
    if (n_read <= 0)
      continue;

    // average power over the whole block, and count how many samples are
    // slammed against the rails so we can tell when the front end saturates
    double sum_pow = 0.0;
    long pairs = 0;
    long clipped = 0;
    for (int i = 0; i + 1 < n_read; i += 2) {
      uint8_t bi = buf[i];
      uint8_t bq = buf[i + 1];
      if (bi <= 1 || bi >= 254 || bq <= 1 || bq >= 254)
        clipped++;
      float I = (float)bi - 127.5f;
      float Q = (float)bq - 127.5f;
      sum_pow += (double)(I * I + Q * Q);
      pairs++;
    }
    if (pairs == 0)
      continue;

    float mean_pow = (float)(sum_pow / (double)pairs);
    float power_db = 10.0f * log10f(mean_pow + 1e-9f);
    float clip_pct = 100.0f * (float)clipped / (float)pairs;

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

      // only let the baseline drift while nothing is happening, otherwise it
      // chases the signal and the detection disappears underneath it
      baseline_db = ALPHA * power_db + (1.0f - ALPHA) * baseline_db;
    }

    // heartbeat, so there are real numbers to look at instead of guesswork
    if (block_count % HEARTBEAT_BLOCKS == 0) {
      fprintf(stderr,
              "block %ld  power=%.2f dB  baseline=%.2f dB  clip=%.2f%%\n",
              block_count, power_db, baseline_db, clip_pct);
    }
  }

  rtlsdr_close(dev);
  return 0;
}
