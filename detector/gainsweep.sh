#!/bin/bash
set -u
CAL=45                       # shorter than 120 for a sweep
sed -i "s/^#define CAL_SECONDS .*/#define CAL_SECONDS     $CAL/" meteor_detect.c

for g in 496 421 364 297 229 166 87; do
  if pgrep -f 'rtl_power|meteor_detect' >/dev/null; then
    echo "device busy, aborting"; exit 1
  fi
  sed -i "s/^#define GAIN_TENTH_DB .*/#define GAIN_TENTH_DB   $g/" meteor_detect.c
  gcc -O2 -o meteor_detect meteor_detect.c -lrtlsdr -lfftw3f -lm || exit 1
  printf '=== gain %s.%s dB ===\n' "${g%?}" "${g: -1}"
  timeout $((CAL+25)) ./meteor_detect --cal 2>&1 \
    | grep -E 'carrier bin|share|WARNING'
  sleep 3
done
