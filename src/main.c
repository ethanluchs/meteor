#include <rtl-sdr.h>
#include <stdio.h>
#include <stdlib.h>

int main(void) {
  uint32_t count = rtlsdr_get_device_count();
  uint8_t buf[16384]; // buffer to read raw rtlsdr bytes into
  int n_read;
  rtlsdr_read_sync(dev, buf, sizeof(buf), &n_read);

  if (count == 0) {
    fprintf(stderr, "No devices found\n");
    return 1;
  }

  for (uint32_t i = 0; i < count; i++) {
    const char *name = rtlsdr_get_device_name(i);
    printf("Device %u is a %s\n", i, name);
  }

  return 0;
}

void print_bytes(void *data, int len) {
  //*data means that we are given the value immediately
  // so we can just reference *data and print from that

  for (int i = 0; i < len; i++) {

    printf("%x", *arr[i]);
  }
  return 0;
}
