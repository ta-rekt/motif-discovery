#!/bin/bash

for j in {0..21}
do
  let startIndex=j*5
  let endIndex=j*5+5
  echo "python compare_speed_fns.py 402 6 $startIndex $endIndex;" >> size6_yu_job1_cmd_file
done
