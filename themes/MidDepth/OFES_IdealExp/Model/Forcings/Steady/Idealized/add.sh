#!/bin/bash
#
#   Add 90-dy forcing to 180-day and generate 180day+90 day forcing

dir="./netcdfs/"

for var in taux tauy; do
  for year in `seq 2000 2029`; do

    echo "year = ${year}"
    infile1="${dir}${var}_180dy_${year}.nc"
    infile2="${dir}${var}_90dy_${year}.nc"
    outfile="${dir}${var}_180dy+90dy_${year}.nc"
    cdo add ${infile1} ${infile2} ${outfile}

  done
done

exit
