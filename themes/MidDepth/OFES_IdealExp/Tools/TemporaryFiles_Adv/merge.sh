#!/bin/bash

#ext="_eq_box_bounded_180dy+90dy"
#ext="_eq_box_bounded_180dy"
ext="_eq_box_bounded_BottomDamp_180dy"

vars="u v w_onU dudx dudy dudz"

for var_name in ${vars}; do
  cdo mergetime ${var_name}${ext}"_????.nc" ${var_name}${ext}".nc"
done

exit

