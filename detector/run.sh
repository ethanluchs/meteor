#!/bin/bash
mkdir -p logs
for i in $(seq 1 200); do
  echo "--- attempt $i $(date -Is) ---" | tee -a logs/run.log
  ./meteor_detect 2>&1 | tee -a logs/run.log
  rc=${PIPESTATUS[0]}
  echo "--- exit $rc $(date -Is) ---" | tee -a logs/run.log
  [ $rc -eq 0 ] && break
  sleep 30
done
