#!/bin/bash

indir="/S/data01/G4006/y0330/data/ofes_exp/results/box_bounded/"
exp_name="180dy+90dy"
years="2029"
#years="2000 2029"
vars="temp salinity"

for var in ${vars}; do
  for year in ${years}; do
    cdo timmean \
        -mermean \
        -zonmean \
        -sellonlatbox,-5,80,-5,5 \
        -seltimestep,-360/-1 \
        -selvar,${var} \
        ${indir}/${exp_name}/${year}/${var}.nc \
        ${var}_mean_${year}_${exp_name}.nc
  done
done

exit
