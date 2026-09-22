#!/bin/bash
set -euo pipefail

exp_name="180dy+90dy"
indir_rho_model="/S/data01/G4006/y0330/data/ofes_exp/results/box_bounded/${exp_name}"
infile_BckGrnd="rho_2000_l1.nc"

export indir_rho_model infile_BckGrnd

ncores=$(nproc)
#max_jobs=$(( ncores > 1 ? ncores / 2 : 1 ))
max_jobs="15"                 #  conservative choice

seq 2000 2029 | xargs -n 1 -P "${max_jobs}" bash -c '
  set -e

  year=$1
  rho_in="${indir_rho_model}/${year}/rho.nc"

  cdo fldmean -mul \
      -sub "${rho_in}" "${infile_BckGrnd}" \
      -sub "${rho_in}" "${infile_BckGrnd}" \
      "rho_diff_sq_hmean_${year}.nc"
' _

