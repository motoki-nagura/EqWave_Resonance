#!/bin/bash
#               Get the initial condition of rho

exp_name="180dy+90dy"

infile="/S/data01/G4006/y0330/data/ofes_exp/results/box_bounded/"${exp_name}"/2000/rho.nc"

cdo seltimestep,1 ${infile} rho_2000_l1.nc

