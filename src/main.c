#include <math.h>
#include <rtl-sdr.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main(void) {

  rtlsdr_dev_t *dev;
  int status = rtlsdr_open(&dev, 0);
  time_t start = time(NULL);
  float baseline = 0;
  int baseline_initialized = 0;
  float alpha = 0.01f;

  if (status != 0) {
    fprintf(stderr, "Failed to open device\n");
    return 1;
  }
  rtlsdr_set_sample_rate(dev, 2048000);
  rtlsdr_set_center_freq(dev, 88500000);
  rtlsdr_set_tuner_gain_mode(dev, 1);
  rtlsdr_reset_buffer(dev);
  uint8_t buf[16384];
  int n_read;

  while (time(NULL) - start < 900) {
    rtlsdr_read_sync(dev, buf, sizeof(buf), &n_read);

    for (int i = 0; i + 1 < n_read; i += 2) {
      float i_val = (float)buf[i] - 127.5f;
      float q_val = (float)buf[i + 1] - 127.5f;
      float magnitude = sqrtf(i_val * i_val + q_val * q_val);
      if (!baseline_initialized) {
        baseline = magnitude;
        baseline_initialized = 1;
      } else {
        baseline = alpha * magnitude + (1 - alpha) * baseline;
      }

      if (baseline_initialized && magnitude > baseline * 5.0f) {
        printf("Possible detection: magnitude=%.2f baseline=%.2f\n", magnitude,
               baseline);
      }
    }
  }
  rtlsdr_close(dev);
  return 0;
}
