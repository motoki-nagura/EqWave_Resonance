import numpy as np
from collections import defaultdict
import sys, pathlib, os

from common.sub_NetCDF_IO import read_variable_from_netcdf
from themes.MidDepth.OFES_IdealExp.Tools.ResonantTrio_tools.sub_Read_ResonantTrios import (
        Read_ResonantTrios,)



def Read_Processed_Bispectrum_NaNAvg_Sum(INPUT_BS_FILE, cluster_considered=None):
    #-----------------
    #     Read in
    #-----------------
    BSR_m1k1, c = read_variable_from_netcdf(INPUT_BS_FILE, "BSR_ω1ω2_m1k1", return_coords=True)
    BSR_m2k2    = read_variable_from_netcdf(INPUT_BS_FILE, "BSR_ω1ω2_m2k2")
    BSR_m3k3    = read_variable_from_netcdf(INPUT_BS_FILE, "BSR_ω1ω2_m3k3")
    groups      = read_variable_from_netcdf(INPUT_BS_FILE, "groups_cluster")
    params      = read_variable_from_netcdf(INPUT_BS_FILE, "parameters")

    #print("BSR_m1k1.shape = ", BSR_m1k1.shape)
    # BSR_m1k1.shape =  (25, 57, 47)

    kwvn, mwvn = c["k"], c["m"]
    scale, k3_min, k3_max, m3_min, m3_max = params

    #   convert "groups" to a dictionary
    groups_cluster = defaultdict(list)
    cluster_ids = np.arange(groups.shape[0])
    for cluster, members in zip(cluster_ids, groups):
        groups_cluster[int(cluster)] = members[members != -999].tolist()
    print("groups_cluster = ", groups_cluster)

    #-------------------------------
    #        Average and Sum
    #-------------------------------
    size = (cluster_considered, BSR_m1k1.shape[1], BSR_m1k1.shape[2],)
    BSR_m1k1_avg = np.zeros(size)
    BSR_m2k2_avg = np.zeros(size)
    BSR_m3k3_sum = np.zeros(size)

    BSR_m1k1_NaN = np.where( BSR_m1k1 > 0., BSR_m1k1, np.nan )  #  Replace 0 with NaN
    BSR_m2k2_NaN = np.where( BSR_m2k2 > 0., BSR_m2k2, np.nan )

    for i_cluster, indices in groups_cluster.items():
        if i_cluster < cluster_considered:
            BSR_m1k1_avg[i_cluster,:,:] = np.nanmean(BSR_m1k1_NaN[indices,:,:], axis=0)
            BSR_m2k2_avg[i_cluster,:,:] = np.nanmean(BSR_m2k2_NaN[indices,:,:], axis=0)
            BSR_m3k3_sum[i_cluster,:,:] = np.sum(    BSR_m3k3[    indices,:,:], axis=0)

    return BSR_m1k1_avg, BSR_m2k2_avg, BSR_m3k3_sum, \
           kwvn, mwvn, scale, k3_min, k3_max, m3_min, m3_max



def Process_ResonantTrios_For_Cluster(
        set_N2, set_freq, sets_n1, sets_add, sets_cluster, km12_range_clus, PATH_OUTPUT):
    """
         Output

         ω123[3]          : ω1, ω2, ω3

         m_filtered[3,:]  : Lists of m and k for which (k1,m1) is in (k1_range,m1_range) and
         k_filtered[3,:]  :                            (k2,m2) is in (k2_range,m2_range)

         m_remaining[3,:] : Lists of m and k for which (k1,m1) is NOT in (k1_range,m1_range) or
         k_remaining[3,:] :                            (k2,m2) is NOT in (k2_range,m2_range)
    """

    for iloop in range(len(sets_cluster)):
        #    Region for large positive bispectrum (obtained from Analysis/CrossBispectrum)
        i_cluster = None
        match sets_cluster[iloop]:
            case "C1": i_cluster = 0
            case "C2": i_cluster = 1
            case "C3": i_cluster = 2

        k1_range = km12_range_clus["k1"][:,i_cluster]
        m1_range = km12_range_clus["m1"][:,i_cluster]
        k2_range = km12_range_clus["k2"][:,i_cluster]
        m2_range = km12_range_clus["m2"][:,i_cluster]

        #--------------------------------------------
        #      Read resonant k1,k2,k3, m1,m2,m3
        #--------------------------------------------

        INFILE_RESONANCE = \
            PATH_OUTPUT+"compute/results/"+set_N2+"/"+set_freq+"/"+sets_n1[iloop]+\
            "/out"+sets_add[iloop]+".nc"

        ω123, mwvn,kwvn, nmer = Read_ResonantTrios(INFILE_RESONANCE)

        #------------------------------------------------------------------------------------#
        #                               Mask results
        #------------------------------------------------------------------------------------#

        #  Eliminate triads with n2 = 2
        mask = (nmer[1,:] != 2)

        mwvn = mwvn[:,mask]
        kwvn = kwvn[:,mask]

        # Mask for a certain range of k and m
        mask = \
        (
          (kwvn[0,:] > np.min(k1_range)) & (kwvn[0,:] < np.max(k1_range)) & # (k1,m1) is in the prescribed range
          (mwvn[0,:] > np.min(m1_range)) & (mwvn[0,:] < np.max(m1_range)) &
          (kwvn[1,:] > np.min(k2_range)) & (kwvn[1,:] < np.max(k2_range)) & # (k2,m2) is in the prescribed range
          (mwvn[1,:] > np.min(m2_range)) & (mwvn[1,:] < np.max(m2_range))
        )

        m_filtered = mwvn[:,mask]       #   in the range
        k_filtered = kwvn[:,mask]

        m_remaining = mwvn[:,~mask]     #   outside the range
        k_remaining = kwvn[:,~mask]

        #--------------------------------------------
        #      Save
        #--------------------------------------------
        match sets_cluster[iloop]:
            case "C1":
                m_filtered_C1, m_remaining_C1 = m_filtered, m_remaining
                k_filtered_C1, k_remaining_C1 = k_filtered, k_remaining
            case "C2":
                m_filtered_C2, m_remaining_C2 = m_filtered, m_remaining
                k_filtered_C2, k_remaining_C2 = k_filtered, k_remaining
            case "C3":
                m_filtered_C3, m_remaining_C3 = m_filtered, m_remaining
                k_filtered_C3, k_remaining_C3 = k_filtered, k_remaining

    if   len(sets_cluster) == 3:
        return ω123, \
                m_filtered_C1, m_remaining_C1, k_filtered_C1, k_remaining_C1, \
                m_filtered_C2, m_remaining_C2, k_filtered_C2, k_remaining_C2, \
                m_filtered_C3, m_remaining_C3, k_filtered_C3, k_remaining_C3
    elif len(sets_cluster) == 2:
        return ω123, \
                m_filtered_C1, m_remaining_C1, k_filtered_C1, k_remaining_C1, \
                m_filtered_C2, m_remaining_C2, k_filtered_C2, k_remaining_C2
    else:
        raise ValueError("Number of cluster must be 2 or 3")


#   ------->>>  NO LONGER USED
def Read_Processed_Bispectrum_avg_sum(INPUT_BS_FILE, cluster_considered=None):
    """
    Superseded by "Read_Processed_Bispectrum_NaNAvg_Sum"
    """
    #-----------------
    #     Read in
    #-----------------
    BSR_m1k1, c = read_variable_from_netcdf(INPUT_BS_FILE, "BSR_ω1ω2_m1k1", return_coords=True)
    BSR_m2k2    = read_variable_from_netcdf(INPUT_BS_FILE, "BSR_ω1ω2_m2k2")
    BSR_m3k3    = read_variable_from_netcdf(INPUT_BS_FILE, "BSR_ω1ω2_m3k3")
    groups      = read_variable_from_netcdf(INPUT_BS_FILE, "groups_cluster")
    params      = read_variable_from_netcdf(INPUT_BS_FILE, "parameters")

    #print("BSR_m1k1.shape = ", BSR_m1k1.shape)
    # BSR_m1k1.shape =  (25, 57, 47)

    kwvn, mwvn = c["k"], c["m"]
    scale, k3_min, k3_max, m3_min, m3_max = params

    #   convert "groups" to a dictionary
    groups_cluster = defaultdict(list)
    cluster_ids = np.arange(groups.shape[0])
    for cluster, members in zip(cluster_ids, groups):
        groups_cluster[int(cluster)] = members[members != -999].tolist()
    print("groups_cluster = ", groups_cluster)

    #-------------------------------
    #        Average and Sum
    #-------------------------------
    size = (cluster_considered, BSR_m1k1.shape[1], BSR_m1k1.shape[2],)
    BSR_m1k1_avg = np.zeros(size)
    BSR_m2k2_avg = np.zeros(size)
    BSR_m3k3_sum = np.zeros(size)

    for i_cluster, indices in groups_cluster.items():
        if i_cluster < cluster_considered:
            for index in indices:
                BSR_m1k1_avg[i_cluster,:,:] += BSR_m1k1[index,:,:]
                BSR_m2k2_avg[i_cluster,:,:] += BSR_m2k2[index,:,:]
                BSR_m3k3_sum[i_cluster,:,:] += BSR_m3k3[index,:,:]

            BSR_m1k1_avg[i_cluster,:,:] /= len(indices)
            BSR_m2k2_avg[i_cluster,:,:] /= len(indices)

    return BSR_m1k1_avg, BSR_m2k2_avg, BSR_m3k3_sum, \
           kwvn, mwvn, scale, k3_min, k3_max, m3_min, m3_max
