"""
  Compute synthetic cross bispectrum, for which the phase is randomized
     for the first variable

     * Time series is obtained from OFES output
     * Time required:
         Read in   : About  1 min
         Test      : About  7 min for NUMBER_OF_SAMPLES=2
                           34 min for NUMBER_OF_SAMPLES=20
                          150 min for NUMBER_OF_SAMPLES=100
                          434 min for NUMBER_OF_SAMPLES=300
                          723 min for NUMBER_OF_SAMPLES=500
                     Time requireed for test = 4 min + 1.5 min x NUMBER_OF_SAMPLES
         Write out : About 10 min
"""
from rich import traceback
traceback.install()

import numpy as np
import time, datetime
import sys, pathlib, os
import subprocess

MY_WORK_ROOT = str( pathlib.Path(os.environ["MY_WORK_ROOT"]) )
PATH_IMPORT = MY_WORK_ROOT+"/code/themes/MidDepth/OFES_IdealExp/Analysis/"

from themes.MidDepth.OFES_IdealExp.Tools.Bispectrum_tools.CrossBispectrum import(
        block_averaged_bispectrum, )
from themes.MidDepth.OFES_IdealExp.Tools.Bispectrum_tools.sub_read_in_u_dudx_u import(
        Input_Parameters, Read_In_U_DUDX_U, )
from common.sub_NetCDF_IO import write_variable_to_netcdf as write_NC

#==================================================================#
#                            Parameters
#==================================================================#


NUMBER_OF_SAMPLES = 500    #   Number of synthetic estimates


def set_params():
    # ------   Fixed parameters  ----- #
    π       = np.pi 
    Erad    = 6371e3             #  Earth's radius [m]
    deg2met = 2. *π *Erad /360.  #  1 deg in meters
    day2sec = 86400.             #  1 day in seconds

    # ------   Input file name  ----- #
    #EXP_NAME = "180dy+90dy"
    EXP_NAME = "180dy"

    dir_in = MY_WORK_ROOT+"/data/themes/MidDepth/OFES_IdealExp/Tools/TemporaryFiles_Adv/"
    extension = "_eq_box_bounded_"+EXP_NAME+".nc"

    # ------  Spatial range of the analysis  ----- #
    z0, z1 = 500., 5000.        #   Depth [m]

    x0, x1 = None, None         #  Longitude [deg]

    # ------  Temporal range of the analysis, and block size  ----- #
    block_size = 1024
    l0,l1  = -4096, None        #  Read the last "l0" timesteps

    # ------  ω, m, k ranges to compute  ----- #
    ω_range = ( -2.*π/(70.*day2sec), -2.*π/(250.*day2sec) ) #  180 day forcing [rad s^-1]

    m_range = ( -4.0e-2, 4.0e-2) #  [rad m^-1]
    k_range = ( -2.3e-5, 2.0e-5) #  [rad m^-1]

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
def show_elapsed_time(process_name,elapsed_time_start):
    elapsed_time_end = time.perf_counter()
    elapsed_time = (elapsed_time_end - elapsed_time_start) /60.
    print(f"****       "+process_name+f"   Elapsed time: {elapsed_time:.2f} minutes"+
          " "*5 + "Date: " + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
          )


def send_message_to_smart_phone():
    subprocess.run(
        ["curl", "-d", "Bispectrum, test using synthetic data", "https://ntfy.sh/EA"],
        check=True
    )


def main():
    elapsed_time_start = time.perf_counter()

    #------------------------------------------------------------------#
    #                       Output file name
    #------------------------------------------------------------------#
    p = set_params()

    OUT_DIR = \
      "/A/data10/nagura/work/data/themes/MidDepth/OFES_IdealExp/Analysis/Bispectrum/Adv_terms/"+\
      "calc_Synth/"+p.EXP_NAME+"/"

    #------------------------------------------------------------------#
    #                           Compute
    #------------------------------------------------------------------#
#    for COMPONENT in ["Zonal_Advection","Meridional_Advection","Vertical_Advection"]:
    for COMPONENT in ["Zonal_Advection",]:

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

        show_elapsed_time("Read In   ",elapsed_time_start)

        #    Carry out Test

        print(f"        NUMBER_OF_SAMPLES = {NUMBER_OF_SAMPLES:d}")

        ω, m, k, Bispec, Bicohe, Bispec_Re_vari = \
           block_averaged_bispectrum(
               var1, var2, var3, dt, dz, dx,
               p.block_size,
               p.ω_range, p.m_range, p.k_range,
               save_lookup_table=False,
               bicoherence_option="Hinich_and_Wolinsky_2005",
               test_synthetic=NUMBER_OF_SAMPLES,
           )
        show_elapsed_time(f"Test done ",elapsed_time_start)

        #------------------------------------------------------------------#
        #                           Write Out
        #------------------------------------------------------------------#
        c = { 0: {"name": "ω1", "values": ω},
              1: {"name": "m1", "values": m},
              2: {"name": "k1", "values": k},
              3: {"name": "ω2", "values": ω},
              4: {"name": "m2", "values": m},
              5: {"name": "k2", "values": k},
            }

        write_NC( Bispec.real,    OUTFILE_NAME, "BSR_mean",    coord_info=c, mode="w", )
        write_NC( Bicohe,         OUTFILE_NAME, "bicoherence", coord_info=c, mode="a", )
        write_NC( Bispec_Re_vari, OUTFILE_NAME, "BSR_vari",    coord_info=c, mode="a", )
        show_elapsed_time("Write out",elapsed_time_start)

    send_message_to_smart_phone()

#==================================================================#
#                           Entry Point 
#==================================================================#
if __name__ == "__main__":
    sys.exit(main())

#  ---------->  No longer used  <------------ {{{
#        #    Loop (in each loop, the phase is Ahat is randomized)
#
#        BispecR_acc    = None
#        Bicoher_acc    = None
#        BispecR_sq_acc = None
#        Bicoher_sq_acc = None
#
#        for iloop in range(NUMBER_OF_ITERATIONS):
#
#            #  Compute
#            ω, m, k, Bispec, Bicohe = block_averaged_bispectrum(
#                var1, var2, var3, dt, dz, dx,
#                p.block_size,
#                p.ω_range, p.m_range, p.k_range,
#                save_lookup_table=False,
#                bicoherence_option="Hinich_and_Wolinsky_2005",
#                test_synthetic=True,
#            )
#
#            #  Accumulate results
#            if BispecR_acc is None:
#                BispecR_acc    = Bispec.real
#                Bicoher_acc    = Bicohe
#                BispecR_sq_acc = Bispec.real **2
#                Bicoher_sq_acc = Bicohe **2
#            else:
#                BispecR_acc    += Bispec.real
#                Bicoher_acc    += Bicohe
#                BispecR_sq_acc += Bispec.real **2
#                Bicoher_sq_acc += Bicohe **2
#
#            show_elapsed_time(f"Iteration #{iloop+1:d}",elapsed_time_start)
#
#        #    Compute average
#
#        BispecR_mean = BispecR_acc    / NUMBER_OF_ITERATIONS
#        Bicoher_mean = Bicoher_acc    / NUMBER_OF_ITERATIONS
#        BispecR_vari = BispecR_sq_acc / NUMBER_OF_ITERATIONS
#        Bicoher_vari = Bicoher_sq_acc / NUMBER_OF_ITERATIONS
#
#        #------------------------------------------------------------------#
#        #                           Write Out
#        #------------------------------------------------------------------#
#        print("OUT >> ",OUTFILE_NAME)
#
#        coord_info = { 0: {"name": "ω1", "values": ω},
#                       1: {"name": "m1", "values": m},
#                       2: {"name": "k1", "values": k},
#                       3: {"name": "ω2", "values": ω},
#                       4: {"name": "m2", "values": m},
#                       5: {"name": "k2", "values": k},
#                     }
#
#        write_NC( BispecR_mean, OUTFILE_NAME, "BSR_mean", coord_info=coord_info, mode="w", )
#        write_NC( Bicoher_mean, OUTFILE_NAME, "BiC_mean", coord_info=coord_info, mode="a", )
#        write_NC( BispecR_vari, OUTFILE_NAME, "BSR_vari", coord_info=coord_info, mode="a", )
#        write_NC( Bicoher_vari, OUTFILE_NAME, "BiC_vari", coord_info=coord_info, mode="a", )
#        write_NC( np.array([NUMBER_OF_ITERATIONS]),
#                    OUTFILE_NAME, "number_of_interations",
#                    coord_info={0: {"name": "p", "values": np.array([0.])},},
#                    mode="a",
#                 )
#        show_elapsed_time("Write out",elapsed_time_start)
#  ---------->  No longer used  <------------ }}}
