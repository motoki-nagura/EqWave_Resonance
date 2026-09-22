"""
     Compute variance of filtered zonal velocity obtained from Exp140day and Exp220day
"""
from rich import traceback
traceback.install()
import os

from themes.MidDepth.OFES_IdealExp.Analysis.Filtering.filtering.Below500m.main_filtering import (
        InputParams, compute_and_write, )

#======================================================================
#            Entry Point
#======================================================================

if __name__ == "__main__":

    grid_name   = "box_bounded"
    period_name = "2000-2029"
    time_width  = 1024         #  # of samples in the time direction
    zstrt, zend = 500., 5000.  #  Depth range
    lstrt, lend = -4096, None  #  Last 4096 days
    N2          = 5e-6         #  Buoyancy frequency squared [s^-2]; for the dispersion relation

    current_dir = os.getcwd()                      #  Get current directory path
    OUT_DIR = current_dir.replace("code", "data")  #  Replace "code" with "data"

    for T in [160., 200.]:
    #for T in [140., 160., 200., 220.]:

        print("\n"+f"  ======== {int(T):d} days ===========")

        input_params = \
                InputParams( grid_name, f"{int(T):d}dy", period_name, time_width,
                             zstrt,zend, lstrt,lend, N2,
                             T/3.,            #  T_target (Target period of ω3 [day])
                             OUT_DIR,         #  Output directory
                            )

        compute_and_write(input_params)

        txt = f"{int(T):d}dy"
        os.rename(OUT_DIR+"/tmp_psd.nc",          OUT_DIR+"/tmp_psd_"+txt+".nc")
        os.rename(OUT_DIR+"/tmp_variance.nc",     OUT_DIR+"/tmp_variance_"+txt+".nc")
        os.rename(OUT_DIR+"/tmp_psd_filtered.nc", OUT_DIR+"/tmp_psd_filtered_"+txt+".nc")
