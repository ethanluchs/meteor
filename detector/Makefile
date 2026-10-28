CC ?= gcc
CFLAGS = -O2 -Wall -Wextra
LDLIBS = -lrtlsdr -lfftw3f -lm

meteor_detect: meteor_detect.c
	$(CC) $(CFLAGS) -o $@ $< $(LDLIBS)

clean:
	rm -f meteor_detect

.PHONY: clean
