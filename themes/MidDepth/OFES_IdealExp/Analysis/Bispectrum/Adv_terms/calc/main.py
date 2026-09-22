"""
  Compute cross bispectrum between u dudx u^*

     * Time series is obtained from OFES output
"""
from rich import traceback
traceback.install()

import numpy as np
#import xarray as xr
import time
#from dataclasses import dataclass, field

import sys, pathlib, os

MY_WORK_ROOT = str( pathlib.Path(os.environ["MY_WORK_ROOT"]) )
PATH_IMPORT = MY_WORK_ROOT+"/code/themes/MidDepth/OFES_IdealExp/Analysis/"

from themes.MidDepth.OFES_IdealExp.Tools.Bispectrum_tools.CrossBispectrum import(
        block_averaged_bispectrum, save_to_netcdf, )
from themes.MidDepth.OFES_IdealExp.Tools.Bispectrum_tools.sub_read_in_u_dudx_u import(
        Input_Parameters, Read_In_U_DUDX_U, )


#==================================================================#
#                            Parameters
#==================================================================#

def set_params():
    # ------   Fixed parameters  ----- #
    π       = np.pi 
    Erad    = 6371e3             #  Earth's radius [m]
    deg2met = 2. *π *Erad /360.  #  1 deg in meters
    day2sec = 86400.             #  1 day in seconds

    # ------   Input file name  ----- #
    #EXP_NAME = "180dy+90dy"
    #EXP_NAME = "180dy"
    EXP_NAME = "BottomDamp_180dy"

    dir_in = MY_WORK_ROOT+"/data/themes/MidDepth/OFES_IdealExp/Tools/TemporaryFiles_Adv/"
    extension = "_eq_box_bounded_"+EXP_NAME+".nc"

    # ------  Spatial range of the analysis  ----- #
    z0, z1 = 500., 5000.        #   Depth [m]

    x0, x1 = None, None         #  Longitude [deg]

    # ------  Temporal range of the analysis, and block size  ----- #
    block_size = 1024
    l0,l1  = -4096, None        #  Read the last "l0" timesteps

    #block_size = 128            #  Test
    #l0, l1 = 0, 366

    #l0, l1 = 6000, 10096        #  Day 6000 ~ 10095

    # ------  ω, m, k ranges to compute  ----- #
    ω_range = ( -2.*π/(70.*day2sec), -2.*π/(250.*day2sec) ) #  180 day forcing [rad s^-1]
    #ω_range = ( -2.*π/( 65.*day2sec), -2.*π/( 75.*day2sec),  #  140 day forcing [rad s^-1]
    #            -2.*π/(130.*day2sec), -2.*π/(150.*day2sec))
    #ω_range = ( -2.*π/(105.*day2sec), -2.*π/(115.*day2sec),  #  220 day forcing [rad s^-1]
    #            -2.*π/(200.*day2sec), -2.*π/(240.*day2sec))
    #ω_range = ( 2.*π/(250.*day2sec), 2.*π/(70.*day2sec) )  #  test

    m_range = ( -4.0e-2, 4.0e-2)                          #  [rad m^-1]
    k_range = ( -2.3e-5, 2.0e-5)                          #  [rad m^-1]

    #  TEST
    #ω_range = (-5e-6, 5e-6)   #  [rad s^-1]
    #m_range = (-1e-2, 1e-2)   #  [rad m^-1]
    #k_range = (-1e-5, 1e-5)   #  [rad m^-1]

    params = Input_Parameters(
                π=π, deg2met=deg2met, day2sec=day2sec,
                EXP_NAME=EXP_NAME, dir_in=dir_in, extension=extension,
                z0=z0, z1=z1, x0=x0, x1=x1, block_size=block_size, l0=l0, l1=l1,
                ω_range=ω_range, m_range=m_range, k_range=k_range,
        )

    return params

#==================================================================#
#                            Main
#==================================================================#

def main():
    elapsed_time_start = time.perf_counter()

    # ------- Output file name --------#
    p = set_params()
    
    OUT_DIR = \
      "/A/data10/nagura/work/data/themes/MidDepth/OFES_IdealExp/Analysis/Bispectrum/"+\
      "Adv_terms/calc/"+p.EXP_NAME+"/"

    #OUTFILE_NAME = "out_UUxU_PositiveOmega.nc"
    #OUTFILE_NAME = "out_UUxU_test.nc"

    #------------------------------------------------------------------#
    #                           Compute
    #------------------------------------------------------------------#
    #LIST_COMPONENTS = ["Zonal_Advection","Meridional_Advection", "Vertical_Advection"]
    LIST_COMPONENTS = ["Zonal_Advection"]

    for COMPONENT in LIST_COMPONENTS:

        print("\n"+"#####"+10*" "+"COMPONENT = ",COMPONENT+10*" "+"#####\n")
        OUTFILE_NAME = None

        match COMPONENT:
            case "Zonal_Advection":
                p.VAR_NAMES = ["u","dudx", "u"]
                OUTFILE_NAME = OUT_DIR+"out_UUxU.nc"

            case "Meridional_Advection":
                p.VAR_NAMES = ["v","dudy", "u"]
                OUTFILE_NAME = OUT_DIR+"out_VUyU.nc"

            case "Vertical_Advection":
                p.VAR_NAMES = ["w","dudz", "u"]
                OUTFILE_NAME = OUT_DIR+"out_WUzU.nc"

        #  Read in
        var1, var2, var3, dt, dz, dx = Read_In_U_DUDX_U(p)

        #  Compute
        ω, m, k, Bispec, Bicoherence = block_averaged_bispectrum(
            var1, var2, var3, dt, dz, dx,
            p.block_size,
            p.ω_range, p.m_range, p.k_range,
            save_lookup_table=True,
            bicoherence_option="Hinich_and_Wolinsky_2005",
        )

        #  Write out
        print("OUT >> ",OUTFILE_NAME)
        save_to_netcdf(OUTFILE_NAME, ω, m, k, Bispec, Bicoherence)

        print("")

    #------------------------------------------------------------------#
    elapsed_time_end = time.perf_counter()
    elapsed_time = (elapsed_time_end - elapsed_time_start) /60.
    print(f"Elapsed time: {elapsed_time:.2f} minutes")

#==================================================================#
#                           Entry Point 
#==================================================================#
if __name__ == "__main__":
    sys.exit(main())
