#!/bin/bash

for j in {0..59}
do
  let startIndex=j*186
  let endIndex=j*186+186
  echo "python compare_speed.py 413 8 $startIndex $endIndex;" >> size8_job1_cmd_file
done
