#include <rtl-sdr.h>
#include <stdio.h>
#include <stdlib.h>

int main(void) {

  rtlsdr_dev_t *dev;
  int status = rtlsdr_open(&dev, 0);
  if (status != 0) {
    fprintf(stderr, "Failed to open device\n");
    return 1;
  }
  rtlsdr_set_sample_rate(dev, 2048000);
  rtlsdr_set_center_freq(dev, 98000000);
  rtlsdr_set_tuner_gain_mode(dev, 1);
  rtlsdr_reset_buffer(dev);
  uint8_t buf[16384];
  int n_read;
  rtlsdr_read_sync(dev, buf, sizeof(buf), &n_read);

  for (int i = 0; i < 16 && i < n_read; i++) {
    printf("%x", buf[i]);
  }
  printf("\n");
  rtlsdr_close(dev);
  return 0;
}
