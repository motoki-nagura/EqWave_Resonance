#!/bin/bash

dir0="/S/data01/G4006/y0330/data/ofes_exp/results/box_bounded/180dy+90dy/"

export dir0

#max_jobs=$(($(nproc) / 2))
max_jobs="15"                 #  conservative choice

seq 2000 2029 | xargs -n 1 -P ${max_jobs} bash -c '
  year=$1
  infile_u="${dir0}/${year}/u.nc"
  infile_v="${dir0}/${year}/v.nc"

  cdo add \
      -fldmean -mul "${infile_u}" "${infile_u}" \
      -fldmean -mul "${infile_v}" "${infile_v}" \
      "sum_mean_${year}.nc"
' _

