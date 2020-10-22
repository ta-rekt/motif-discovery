#!/bin/bash

for i in 6 33 105 256
do
  echo "python compare_speed.py $i 5 0 5;" >> test_job_size5_cmd_file
done
