#!/bin/bash
for i in 12 13 14 15 16
do
      ./mfc.sh run /storage/scratch1/6/bwilfong3/bodyForces/wenoEpsSurvey/CPU01/$i/case.py -N 3 -n 24 --no-gpu --dry-run -# WEC01$i
done
