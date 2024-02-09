#!/bin/bash

echo "Running test of quaqc output."

../quaqc --output-ext .test.txt -T tss.bed -t target.bed reads.bam

diff <(tail -n116 reads.quaqc.txt) <(tail -n116 reads.test.txt) > diff.txt

if [ -s diff.txt ] ; then
  echo "Test failed, found the following diff:"
  cat diff.txt
  exit 1
else
  echo "Test succeeded, no changes in output."
  rm -f reads.test.txt
  rm -f diff.txt
  exit 0
fi

