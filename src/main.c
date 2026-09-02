#include <rtl-sdr.h>
#include <stdio.h>

int main(void) {
  uint32_t count = rtlsdr_get_device_count();

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
