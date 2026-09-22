"""
   The purpose of this program is to plot real part of bispectrum for a specific range of m3 and k3

   - Read the lookup table for (k3,m3) that satisfies k1+k2=k3 and m1+m2=m3
   - Map the bispectral values onto (k3,m3), keep the maximum values in each grid,and discard the rest
   - Limit k3 and m3 to the prescribed ranges
   - Cluster the bispectrum in (k1,m1) and (k2,m2) space
   - Write out results
   - Draw figures

   ω1 = 2π/(180 days), ω2 = 2π/(90 days), ω3 = 2π/(60 days) are fixed.
"""
from rich import traceback
traceback.install()

import numpy as np
import sys, pathlib, os

from themes.MidDepth.OFES_IdealExp.Tools.Bispectrum_tools.tools_k3m3_region import (
    Read_In_Bispectrum, Read_In_LookUp_Tables, Eliminate_BSR_Range,
    Keep_Max_Discard_Rest, Find_Top_N_Values,
    Substitute_Listed_Values_in_km_space, Pattern_Clustering,
    Save_list_to_NetCDF, Save_BSR_BC_to_NetCDF,
    Write_Out_To_TextFile, Draw_Figure_to_PDF,
    )

MY_WORK_ROOT = str( pathlib.Path(os.environ["MY_WORK_ROOT"]) )

#==========================================================================================
#     Options 
#==========================================================================================

OPTIONS_TARGET_WAVE3 = {"R":     "Rossby",
                        "R m<0": "Rossby_m<0",}

OPTIONS_FILTER = {"none":  "NoSmth",
                  "3box":  "3ptBox",
                  "121":   "121Filt",
                  "12321": "12321Filt",}

OPTIONS_COMPONENT = {"xyz": "xyz", #  u ∂u/dx + v ∂u/∂y + w ∂u/∂z
                     "xz":  "xz",  #  u ∂u/dx + w ∂u/∂z
                     "x":   "x",}  #  u ∂u/dx

EXPERIMENT = "180dy+90dy"

#----------------------

TARGET_WAVE3 = OPTIONS_TARGET_WAVE3["R"] #  Target region in k3-m3 space where top values are selected
FILTER       = OPTIONS_FILTER["3box"]    #  Filtering in k1-m1 and k2-m2 space
COMPONENT    = OPTIONS_COMPONENT["x"]    #  Components of advection terms

#---------------------

k_range = [-2.0e-5, 2.0e-5]     #  Wave number range for figures
m_range = [-4.0e-2, 4.0e-2]     #  Clustering is also done in this range

numb_TopValues           = 75   #  Number of top values to be chosed in (k3,m3) space
threshold_find_TopValues = 0.2  #  BSR smaller than this value is ignored

                                #  Parameters for clustering
cluster_eps             = 1.5   #    Larger values group less similar shapes together
cluster_threshold_ratio = 0.0   #    Threshold for values to be considered
cluster_crop_and_pad    = False #    If true, patterns are cropped and padded

#==========================================================================================
#     Fixed Parameters

π = np.pi
day2sec = 86400.   #  1 day in seconds

ω_180dy, ω_90dy, ω_60dy = 2.*π/(180.*day2sec), 2.*π/(90.*day2sec), 2.*π/(60.*day2sec)
                    #  [rad second^-1]
ω1, ω2 = ω_180dy, ω_90dy
ω3 = ω1 + ω2        #  ω3 is 2π/(60 days)

#-----------------------------------------------
#   Range for k3 and m3

match TARGET_WAVE3:
    case "Rossby":
        k3_target, m3_target = [-1.8e-5, 0.0e-5], [-0.4e-2, 0.4e-2]
    case "Rossby_m<0":
        k3_target, m3_target = [-1.5e-5, -0.5e-5], [-0.4e-2, -0.1e-2]

    ##  (From Filtering/filtering/500-5000m/main.py)
    #case "Rossby":
    #    k3_target, m3_target = [-2.0e-5, -0.2e-5], [-0.6e-2, 0.6e-2]
    #case "Rossby_LowK":
    #    k3_target, m3_target = [-1.1e-5, -0.2e-5], [-0.6e-2, 0.6e-2]
    #case "Rossby_HighK":
    #    k3_target, m3_target = [-2.0e-5, -1.1e-5], [-0.6e-2, 0.6e-2]
    #case "OffResonant":
    #    k3_target, m3_target = [-1.5e-5, -0.3e-5], [ 0.9e-2, 1.9e-2]
    #    #    Off resonant variability for k3 < 0 and m3 > 0

k3_target, m3_target = np.array(k3_target), np.array(m3_target)

#------------------------------------------------------------------
#   file name

INPUT_NETCDF = \
        MY_WORK_ROOT+"/data/themes/MidDepth/OFES_IdealExp/Analysis/Bispectrum/Adv_terms/calc/"+\
        EXPERIMENT+"/"

INPUT_TABLE = MY_WORK_ROOT+"/code/themes/MidDepth/OFES_IdealExp/Analysis/"+\
                "Bispectrum/Adv_terms/calc/lookup_tables.nc"

OUTPUT_DIR = "./results/"+COMPONENT+"/"
OUTPUT_TXT = OUTPUT_DIR+"out_60dy_"+TARGET_WAVE3+"_"+FILTER+".txt"
OUTPUT_NC  = OUTPUT_DIR+"out_60dy_"+TARGET_WAVE3+"_"+FILTER+".nc"
OUTPUT_PDF = OUTPUT_DIR+"fig_60dy_"+TARGET_WAVE3+"_"+FILTER+".pdf"

if TARGET_WAVE3 == "Rossby" and FILTER == "3ptBox":
    OUTPUT_LIST_NC = OUTPUT_DIR+"out_list_60dy_"+TARGET_WAVE3+"_"+FILTER+".nc"

#==========================================================================================
#==========================================================================================
def main():
    #======================
    #     Computation
    #======================
    BSR_ω1ω2, k_BSR,m_BSR, scale_BSR_ω1ω2 = \
            Read_In_Bispectrum(INPUT_NETCDF, ω1,ω2, FILTER, COMPONENT=COMPONENT)

    m3_table,k3_table = Read_In_LookUp_Tables(INPUT_TABLE, m_BSR,k_BSR)

    BSR_ω1ω2_cond = BSR_ω1ω2.copy()
    BSR_ω1ω2_cond, im3_cond, ik3_cond = \
            Eliminate_BSR_Range("outside",BSR_ω1ω2_cond, m3_table,k3_table, m3_target, k3_target)

    BSR_ω1ω2_cond = Keep_Max_Discard_Rest(BSR_ω1ω2_cond,m_BSR,k_BSR, m3_table,k3_table)

    list_TopNValues = Find_Top_N_Values(BSR_ω1ω2_cond, numb_TopValues, threshold_find_TopValues)

    BSR_ω1ω2_m1k1, BSR_ω1ω2_m2k2, BSR_ω1ω2_m3k3 = \
        Substitute_Listed_Values_in_km_space(BSR_ω1ω2,m_BSR,k_BSR, list_TopNValues,
                                             m3_table,k3_table, m3_target,k3_target)

    groups_cluster, proc_pat_cluster = \
        Pattern_Clustering(BSR_ω1ω2_m1k1,BSR_ω1ω2_m2k2, m_BSR,k_BSR,m_range,k_range,
                           cluster_eps, cluster_threshold_ratio, cluster_crop_and_pad)

    #============================================
    #    Output
    #============================================

    #  Save results to a NetCDF file
    if TARGET_WAVE3 == "Rossby" and FILTER == "3ptBox":
        Save_list_to_NetCDF(OUTPUT_LIST_NC, list_TopNValues,) #  Save list_TopNValues

        Save_BSR_BC_to_NetCDF(
                OUTPUT_NC, BSR_ω1ω2_m1k1,BSR_ω1ω2_m2k2,BSR_ω1ω2_m3k3, m_BSR,k_BSR,
                groups_cluster,scale_BSR_ω1ω2,m3_target,k3_target)
                                                           #  Save BSR and groups_cluster

    #   PDF output of processed patters for clustering
    #save_clusters_to_pdf(patterns=BSR_ω1ω2_m2k2_sub,
    #                     groups=groups_cluster, processed_patterns=proc_pat_cluster,
    #                     output_pdf="clusters.pdf",)

    Write_Out_To_TextFile(OUTPUT_TXT,
                          m3_target, k3_target,
                          m_BSR,k_BSR, m3_table,k3_table, im3_cond,ik3_cond, BSR_ω1ω2_cond,
                          scale_BSR_ω1ω2, groups_cluster, list_TopNValues, BSR_ω1ω2, 
                          numb_TopValues,
                          )

    Draw_Figure_to_PDF(OUTPUT_PDF,
                       m_BSR,k_BSR,groups_cluster,
                       BSR_ω1ω2_m1k1,BSR_ω1ω2_m2k2,BSR_ω1ω2_m3k3,scale_BSR_ω1ω2,
                       m_range,k_range,ω1,ω2,ω3,m3_target,k3_target,
                       ttl_periods=("180","90","60"),
                       )

#==========================================================================
#  Script entry point
#==========================================================================
if __name__ == "__main__":
    sys.exit(main())
