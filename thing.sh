#!/bin/bash
for i in 12 13 14 15 16
do
      ./mfc.sh run /storage/scratch1/6/bwilfong3/bodyForces/wenoEpsSurvey/GPU01/$i/case.py -N 1 -n 1 -e batch --case-optimization --gpu --dry-run -# WEG01$i -j 8
done
