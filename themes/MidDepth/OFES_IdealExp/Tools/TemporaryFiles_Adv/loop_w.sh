#!/bin/bash

grid_name="box_bounded"

#exp_name="180dy+90dy"
#exp_name="180dy"
exp_name="BottomDamp/180dy"

years=`seq 2018 2029`

dir0="/A/data10/nagura/work/data/ofes_exp/results"

lons=0    #   Longitudinal range (full)
lone=70

lats=-0.1 #   Latitudinal range (equator)
lat0=0.0
late=0.1

for year_name in $years; do

  exp_name_1="${exp_name//\//_}"       #  Replace "/" with "_"
  ext_out="_${grid_name}_${exp_name_1}_${year_name}.nc"
  
  dir1="${dir0}/${grid_name}/${exp_name}/${year_name}/"
  
  in_w=${dir1}w.nc
  
  cdo sellonlatbox,${lons},${lone},${lats},${late} $in_w ./w_0p1S-0p1N${ext_out}
  
  python calc_w_interp.py ${ext_out}
done

