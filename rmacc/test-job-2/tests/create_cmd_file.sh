#!/bin/bash

for i in 31 6
do
  for j in {11..20}
  do
    echo "python test-job-2.py $i $j;" >> cmd_file
  done
done
