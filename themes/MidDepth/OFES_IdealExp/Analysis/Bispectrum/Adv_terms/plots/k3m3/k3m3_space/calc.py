"""
   Map the real part of bispectrum onto k3-m3 space
"""
from rich import traceback
traceback.install()

import numpy as np
import sys, os, pathlib

from common.sub_NetCDF_IO import (write_variable_to_netcdf,)

from themes.MidDepth.OFES_IdealExp.Tools.Bispectrum_tools.tools_k3m3_region import(
        Read_In_Bispectrum, Read_In_LookUp_Tables,
        )

MY_WORK_ROOT = str( pathlib.Path(os.environ["MY_WORK_ROOT"]) )

π = np.pi
day2sec = 86400.   #  1 day in seconds
ω_180dy, ω_90dy, ω_60dy = \
    2.*π/(180.*day2sec), 2.*π/(90.*day2sec), 2.*π/(60.*day2sec)   #  [rad second^-1]

#==========================================================================================
#     Input Parameters
#==========================================================================================

#  Model run

EXP_NAME = { "std" :  "180dy+90dy",
             "sens" : "180dy",     }["std"]
#             "sens" : "180dy",     }["sens"]

#  Frequency set

FREQ_NAME = {"ABC":  "180-90-60day",
             "AAB":  "180-180-90day"}["ABC"]
#             "AAB":  "180-180-90day"}["AAB"]

#  Filter

FILTER = "3ptBox"
#FILTER = "NoSmth"

#  Component

COMPONENT = {"xyz": "xyz",  #  u ∂u/dx + v ∂u/∂y + w ∂u/∂z
             "xz":  "xz",   #  u ∂u/dx + w ∂u/∂z
             "x":   "x",    #  u ∂u/dx
             "y":   "y",    #  v ∂u/∂y
             "z":   "z",}   #  w ∂u/∂z
COMPONENT = COMPONENT["x"]

#==========================================================================================
#     Derived Parameters
#==========================================================================================

#  Input

DIR_INPUT = MY_WORK_ROOT+"/data/themes/MidDepth/OFES_IdealExp/Analysis/Bispectrum/"+\
             "Adv_terms/calc/"+EXP_NAME+"/"

DIR_INPUT_sync = \
        MY_WORK_ROOT+"/data/themes/MidDepth/OFES_IdealExp/Analysis/Bispectrum/"+\
        "Adv_terms/calc_Synth/"+EXP_NAME+"/"

FILE_INPUT_TABLE = MY_WORK_ROOT+"/code/themes/MidDepth/OFES_IdealExp/Analysis/"+\
                "Bispectrum/Adv_terms/calc/lookup_tables.nc"

#  Output

OUTFILE_EXT = { "NoSmth" : "NoSmth",
                "3ptBox" : "3ptBox",
              }[FILTER]

OUTDIR_SUB = "/results/Exp_"+EXP_NAME+"/"+FREQ_NAME+"/"+COMPONENT+"/"

current_dir = os.getcwd()                         #  Get current directory path
OUTDIR_DATA = current_dir.replace("code", "data") #  Replace "code" with "data"

OUTFILE_NETCDF = OUTDIR_DATA+OUTDIR_SUB+"out_"+OUTFILE_EXT+".nc"

#  Set ω1, ω2, and k3-m3 region for Rossby

if   FREQ_NAME == "180-90-60day":
    ω1, ω2, = ω_180dy, ω_90dy    #  ω3 is 2π/(60 days)

elif FREQ_NAME == "180-180-90day":
    ω1, ω2, = ω_180dy, ω_180dy,  #  ω3 is 2π/(90 days)

#==========================================================================================
def main():
    #---------------------------
    #   Real Bispectrum
    #---------------------------
    BSR, k_BSR,m_BSR, _ = \
            Read_In_Bispectrum(DIR_INPUT, ω1,ω2, FILTER, COMPONENT=COMPONENT,)

    m3,k3 = Read_In_LookUp_Tables(FILE_INPUT_TABLE, m_BSR,k_BSR)

    BSR_k3m3 = Map_onto_k3m3_space(BSR,k_BSR,m_BSR, m3,k3)

    coord = { 0: {"name": "dummy_axis"},
              1: {"name": "m_BSR", "values": m_BSR},
              2: {"name": "k_BSR", "values": k_BSR},
        }

    write_variable_to_netcdf(
        BSR_k3m3, OUTFILE_NETCDF, "BSR_k3m3", coord_info=coord, mode="w",)

    print("")

    #---------------------------
    #   Other variables
    #---------------------------
    list_var_name_in = \
      ["bicoherence",          "BSR_mean",   "bicoherence",           "BSR_vari", ]
    list_PATH = \
      [    DIR_INPUT,      DIR_INPUT_sync,  DIR_INPUT_sync,       DIR_INPUT_sync, ]
    list_var_name_out = \
      [    "BC_k3m3", "BSR_mean_sync_k3m3", "BC_sync_k3m3", "BSR_vari_sync_k3m3", ]

    #list_var_name_in  = ["bicoherence", ]
    #list_PATH         = [    DIR_INPUT, ]
    #list_var_name_out = [    "BC_k3m3", ]


    for var_name_in, IN_PATH, var_name_out in zip(
            list_var_name_in, list_PATH, list_var_name_out
            ):

        var, k_BSR,m_BSR, _ = \
                Read_In_Bispectrum(IN_PATH, ω1,ω2, FILTER, \
                                   VARIABLE_NAME=var_name_in, COMPONENT=COMPONENT,)

        m3,k3 = Read_In_LookUp_Tables(FILE_INPUT_TABLE, m_BSR,k_BSR)

        var_k3m3  = Map_onto_k3m3_space(var,k_BSR,m_BSR, m3,k3)

        write_variable_to_netcdf(
             var_k3m3, OUTFILE_NETCDF, var_name_out, coord_info=coord, mode="a",)

        print("")


def Map_onto_k3m3_space(var,k_BSR,m_BSR, m3,k3):
    var_mapped = np.full((3000,len(m_BSR),len(k_BSR)), np.nan)

    for i,k3_tar in enumerate(k_BSR):       #   target k3
        for j,m3_tar in enumerate(m_BSR):   #   target m3
            ik1,ik2 = np.argwhere( k3 == k3_tar ).T
            im1,im2 = np.argwhere( m3 == m3_tar ).T

            icnt = 0
            for ii in range(len(ik1)):
                for jj in range(len(im1)):
                    var_mapped[icnt,j,i] = var[ im1[jj], ik1[ii], im2[jj], ik2[ii] ]
                    icnt += 1

    return var_mapped


#  Script entry point
if __name__ == "__main__":
    sys.exit(main())
