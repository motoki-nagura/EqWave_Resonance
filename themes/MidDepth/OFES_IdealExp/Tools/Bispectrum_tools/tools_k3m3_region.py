from rich import traceback
traceback.install()

import numpy as np
import netCDF4
from itertools import product

from .sub_flip_omega_smth import sign_flip__omega_extract__smoothing
from .subs_lookup_table   import reorder, generate_k3, check_ω1ω2ω3
from .sub_find_index      import find_index_largest_values_NdimArray
from .sub_clustering.clustering_2dArray  import cluster_patterns
from common.sub_NetCDF_IO import write_variable_to_netcdf

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.patches as patches
import matplotlib.colors as colors

from common.new_page import new_page
from themes.MidDepth.OFES_IdealExp.Tools.Bispectrum_tools.subs_k3m3_region_PlotTools import (
    manage_page, set_ax, )
from themes.MidDepth.OFES_IdealExp.Tools.common.subs_ofes_idealexp import (
    overplot_DispersionRelation, )

#==========================================================================

def Read_In_Bispectrum(PATH_INPUT_NETCDF, ω1,ω2, FILTER, #{{{
                       VARIABLE_NAME='bispectrum_real', COMPONENT='xyz'):
    """ Read in real bispectrum and bicoherence
           COMPONENT : "xyz" = u ∂u/dx + v ∂u/∂y + w ∂u/∂z
                       "xz"  = u ∂u/dx + w ∂u/∂z
                       "x"   = u ∂u/dx
                       "y"   = v ∂u/dy
                       "z"   = w ∂u/dz

        Bispectrum components are summed.
        Bicoherence components are averaged.
    """

    if VARIABLE_NAME not in {"bispectrum_real", "bicoherence", "BSR_mean", "BSR_vari",}:
        raise ValueError(f"VARIABLE_NAME: {VARIABLE_NAME!r}. "
                         "Choose 'bispectrum_real', 'bicoherence', 'BSR_mean', "
                         "or 'BSR_vari'")

    if COMPONENT not in {"xyz", "xz", "x", "y", "z"}:
        raise ValueError(f"COMPONENT: {COMPONENT!r}. Choose 'xyz', 'xz' 'x', 'y', or 'z'")

    inf_x = PATH_INPUT_NETCDF+'out_UUxU.nc'
    inf_y = PATH_INPUT_NETCDF+'out_VUyU.nc'
    inf_z = PATH_INPUT_NETCDF+'out_WUzU.nc'

    f_x = netCDF4.Dataset(inf_x, 'r')

    if   "y" in COMPONENT:
        f_y = netCDF4.Dataset(inf_y, 'r')
    elif "z" in COMPONENT:
        f_z = netCDF4.Dataset(inf_z, 'r')

    with netCDF4.Dataset(inf_x, "r") as nc:
        if   ("ω" in nc.variables):
            ω = np.copy(f_x.variables['ω'])
            k = np.copy(f_x.variables['k'])
            m = np.copy(f_x.variables['m'])
        elif ("ω1" in nc.variables):
            ω = np.copy(f_x.variables['ω1'])
            k = np.copy(f_x.variables['k1'])
            m = np.copy(f_x.variables['m1'])

    if   "x" in COMPONENT:
        print("IN <<< "+inf_x)
        var_x = np.copy(f_x.variables[VARIABLE_NAME])
        f_x.close()
    elif "y" in COMPONENT:
        print("IN <<< "+inf_y)
        var_y = np.copy(f_y.variables[VARIABLE_NAME])
        f_y.close()
    elif "z" in COMPONENT:
        print("IN <<< "+inf_z)
        var_z = np.copy(f_z.variables[VARIABLE_NAME])
        f_z.close()

    if  VARIABLE_NAME in {"bispectrum_real", "BSR_mean", "BSR_vari", }:
       if   COMPONENT == "xyz":
           var = var_x + var_y + var_z
       elif COMPONENT == "xz":
           var = var_x + var_z
       elif COMPONENT == "x":
           var = var_x
       elif COMPONENT == "y":
           var = var_y
       elif COMPONENT == "z":
           var = var_z

    elif VARIABLE_NAME in {"bicoherence", }:
        if   COMPONENT == "xyz":
            var = (var_x + var_y + var_z) /3.
        elif COMPONENT == "xz":
            var = (var_x + var_z) /2.
        elif COMPONENT == "x":
            var = var_x
        elif COMPONENT == "y":
            var = var_y
        elif COMPONENT == "z":
            var = var_z

    smoothing = { "NoSmth":    None,
                  "3ptBox":    (3,3,3,3),
                  "121Filt":   "121",
                  "12321Filt": "12321" }[FILTER]

    if  VARIABLE_NAME in {"bispectrum_real", "BSR_mean", }:
        sign_flip = True
    elif VARIABLE_NAME in {"BSR_vari", "bicoherence", }:
        sign_flip = False

    var_ω1ω2 = sign_flip__omega_extract__smoothing(
            var, ω, ω1,ω2, sign_flip=sign_flip,  smoothing=smoothing)

    scale_var_ω1ω2 = np.max(var_ω1ω2)
    #print("scale = ", scale_var_ω1ω2)

    return var_ω1ω2, k,m, scale_var_ω1ω2
#}}}1

def Read_In_LookUp_Tables(INPUT_TABLE, m_BSR,k_BSR): #{{{1
    """ Read in Lookup Tables and Generate m3 and k3 Tables """

    print("IN <<< "+INPUT_TABLE)

    f1 = netCDF4.Dataset(INPUT_TABLE, 'r')
    k          = np.copy(f1.variables['k'])
    m          = np.copy(f1.variables['m'])
    k_r        = np.copy(f1.variables['k_r'])
    m_r        = np.copy(f1.variables['m_r'])
    ik3_lookup = np.copy(f1.variables['ik3_lookup'])
    im3_lookup = np.copy(f1.variables['im3_lookup'])
    f1.close()

    #  Reordering
    m_r, im3_lookup = reorder(m_r, im3_lookup)
    k_r, ik3_lookup = reorder(k_r, ik3_lookup)

    #  Check grid consistency
    check_consistency(m_r, m_BSR, 'm_r  ','m_BSR')
    check_consistency(k_r, k_BSR, 'k_r  ','k_BSR')

    #  Generate m3 and k3
    m3 = generate_k3(im3_lookup,m)   #  m3[im1,im2]
    k3 = generate_k3(ik3_lookup,k)   #  k3[ik1,ik2]

    #  Test if ω1[i1] + ω2[i2] = ω3[i1,i2] is satisfied
    check_ω1ω2ω3("k", ω1s=k_r,ω2s=k_r,ω3s=k3)
    check_ω1ω2ω3("m", ω1s=m_r,ω2s=m_r,ω3s=m3)

    return m3,k3
#}}}1

def Eliminate_BSR_Range(tag, BSR, m3,k3, m3_target,k3_target, fill_value=0.): #{{{
    """ Limit BSR to the range
          m3_target_min < m3 < m3_target_max,  k3_target_min < k3 < k3_target_max """

    #    Generate masks

    k3_target, m3_target = np.array(k3_target), np.array(m3_target)

    m_mask = (m3 >= np.min(m3_target)) & (m3 <= np.max(m3_target))
    k_mask = (k3 >= np.min(k3_target)) & (k3 <= np.max(k3_target))

    im = np.where(m_mask)
    ik = np.where(k_mask)

    #    4-dimensional mask
    mask = ( m_mask[:, None, :, None] &
             k_mask[None, :, None, :] )

    if   tag == "outside":
        #  Replace values outside the range
        BSR = np.where(mask, BSR, fill_value)

    elif tag == "inside":
        #  Replace values inside the range
        BSR[mask] = fill_value

    else:
        raise ValueError("Eliminate_BSR_Range: 'tag' must be 'outside' or 'inside'.")

    return BSR, im,ik
#}}}1

def Keep_Max_Discard_Rest(BSR_ω1ω2_cond,m_BSR,k_BSR, m3_table,k3_table): #{{{
    """ Keep the maximum and discard the rest for each (k3,m3) grid """

    arbitrary_negative_number = -1e10

    for k3 in k_BSR:
        for m3 in m_BSR:

            # Find indices satisfying the target conditions
            k_indices = list(zip(*np.where(k3_table == k3)))
            m_indices = list(zip(*np.where(m3_table == m3)))

            # Generate all valid index combinations
            index_list = [
                (im1, ik1, im2, ik2)
                for (ik1, ik2), (im1, im2) in product(k_indices, m_indices)
            ]

            # Convert the index list into NumPy arrays
            im1_idx = np.array([idx[0] for idx in index_list])
            ik1_idx = np.array([idx[1] for idx in index_list])
            im2_idx = np.array([idx[2] for idx in index_list])
            ik2_idx = np.array([idx[3] for idx in index_list])

            # Extract the corresponding values
            values = BSR_ω1ω2_cond[im1_idx, ik1_idx, im2_idx, ik2_idx]

            # Keep only the maximum values
            max_value = values.max()
            mask = values < max_value

            BSR_ω1ω2_cond[
                im1_idx[mask],
                ik1_idx[mask],
                im2_idx[mask],
                ik2_idx[mask]
            ] = arbitrary_negative_number

    return BSR_ω1ω2_cond
#}}}1

def Find_Top_N_Values(BSR_ω1ω2_cond, numb_TopValues, threshold_find_TopValues): #{{{1
    """ Find the top N values in the prescribed (k3,m3) range """

    #   Make a list for the index

    arr   = BSR_ω1ω2_cond              #  BSR limited to the (k3,m3) range
    ntop  = numb_TopValues
    thres = threshold_find_TopValues

    list_TopNValues = find_index_largest_values_NdimArray(arr, ntop, threshold=thres)
                                              #   list_TopNValues[{ntop},{m1,k1,m2,k2}]

    return list_TopNValues
#}}}1

def Substitute_Listed_Values_in_km_space( #{{{1
        arr,m_BSR,k_BSR, list_TopNValues,
        m3_table,k3_table, m3_target,k3_target):
    """ Substitute the top values in the list to (k1,m1), (k2,m2), and (k3,m3) grids
    """

    #   Substitute BSR in the list

    ntop = len(list_TopNValues)

    arr_ω1ω2_m1k1, arr_ω1ω2_m2k2, arr_ω1ω2_m3k3 = np.zeros((3,ntop,len(m_BSR),len(k_BSR)))

    for i_rank, idx in enumerate(list_TopNValues):
        im1,ik1,im2,ik2 = idx

        #  Error check

        #     Are m3 and k3 with in m3_target and k3_target?
        if np.logical_or( np.logical_or( m3_table[im1,im2] < np.min(m3_target),
                                         m3_table[im1,im2] > np.max(m3_target)),
                          np.logical_or( k3_table[ik1,ik2] < np.min(k3_target),
                                         k3_table[ik1,ik2] > np.max(k3_target)) ):
            raise ValueError(
                    "m3_table[im1,im2] or k3_table[ik1,ik2] is out of range!\n"
                    f"   (im1,im2)=({im1:d},{im2:d}),  m3_table[im1,im2]={m3_table[im1,im2]:.2e}\n"
                    f"   (ik1,ik2)=({ik1:d},{ik2:d}),  k3_table[ik1,ik2]={k3_table[ik1,ik2]:.2e}\n"
                    f"   {np.min(m3_target):.2e} <= m3_target <= {np.max(m3_target):.2e}\n"
                    f"   {np.min(k3_target):.2e} <= k3_target <= {np.max(k3_target):.2e}\n"
                    f"   BSR there = {arr[im1,ik1,im2,ik2]:.2e}"
                             )

        #     Are m3 and k3 on the grids of m_BSR and k_BSR?
        im3 = np.argmin( np.abs(m_BSR - m3_table[im1,im2]) )
        ik3 = np.argmin( np.abs(k_BSR - k3_table[ik1,ik2]) )

        if ( ( not np.any(np.isclose(m_BSR, m3_table[im1,im2])) ) or
             ( not np.any(np.isclose(k_BSR, k3_table[ik1,ik2])) )
           ):
            txt_m_BSR = ", ".join(f"{x: .2e}" for x in m_BSR[im3-1:im3+2])
            txt_k_BSR = ", ".join(f"{x: .2e}" for x in k_BSR[ik3-1:ik3+2])
            raise ValueError(
               f"m3_table[im1,im2]={m3_table[im1,im2]: .2e} or "
               f"k3_table[ik1,ik2]={k3_table[im1,im2]: .2e} does not coincide m or k grid\n"
                "m_BSR[-1:1]="+txt_m_BSR+"\n"
                "k_BSR[-1:1]="+txt_k_BSR+"\n"
                )

        #  Substitute

        val = arr[im1,ik1,im2,ik2]
        arr_ω1ω2_m1k1[i_rank,im1,ik1] = val
        arr_ω1ω2_m2k2[i_rank,im2,ik2] = val
        arr_ω1ω2_m3k3[i_rank,im3,ik3] = val

    return arr_ω1ω2_m1k1, arr_ω1ω2_m2k2, arr_ω1ω2_m3k3
#}}}1

def Pattern_Clustering(BSR_ω1ω2_m1k1,BSR_ω1ω2_m2k2, m_BSR,k_BSR, #{{{1
                       m_range,k_range,
                       cluster_eps, cluster_threshold_ratio, cluster_crop_and_pad):
    """ Pattern clustering in (k1,m1) and (k2,m2) space """

    #    Consider BSR only in "k_range" and "m_range"
    ymask = (m_range[0] <= m_BSR) & (m_BSR <= m_range[1])
    xmask = (k_range[0] <= k_BSR) & (k_BSR <= k_range[1])
    BSR_ω1ω2_m1k1_sub = BSR_ω1ω2_m1k1[:, ymask][:, :, xmask]
    BSR_ω1ω2_m2k2_sub = BSR_ω1ω2_m2k2[:, ymask][:, :, xmask]

    #    Combine BSR in (k1,m1) space and (k2,m2) space into one array
    BSR_ω1ω2_sub_concatenated = np.concatenate( (BSR_ω1ω2_m1k1_sub,BSR_ω1ω2_m2k2_sub), axis=1)
    #                                          The results seems to be the same if axis=2

    #   Clustering
    groups_cluster, proc_pat_cluster, _ = \
      cluster_patterns(BSR_ω1ω2_sub_concatenated,
        eps=cluster_eps, threshold_ratio=cluster_threshold_ratio, crop_and_pad=cluster_crop_and_pad,)

    return groups_cluster, proc_pat_cluster
#}}}1

def Save_list_to_NetCDF(OUTPUT_FILE,list,): #{{{1
    print("OUT >>> ",OUTPUT_FILE)
    data = np.array(list, dtype=np.int32)  # Convert to a NumPy array
    with netCDF4.Dataset(OUTPUT_FILE, mode="w", format="NETCDF4") as nc:
        nc.createDimension("record", len(data))
        i_m1 = nc.createVariable("index_m1", "i4", ("record",))
        i_k1 = nc.createVariable("index_k1", "i4", ("record",))
        i_m2 = nc.createVariable("index_m2", "i4", ("record",))
        i_k2 = nc.createVariable("index_k2", "i4", ("record",))
        i_m1[:] = data[:, 0]
        i_k1[:] = data[:, 1]
        i_m2[:] = data[:, 2]
        i_k2[:] = data[:, 3]
        nc.description = "List of m1, k1, m2, and k2 values"
#}}}1

def Save_BSR_BC_to_NetCDF(OUTPUT_NC,  #{{{1
                       BSR_ω1ω2_m1k1,BSR_ω1ω2_m2k2,BSR_ω1ω2_m3k3, m_BSR,k_BSR,
                       groups_cluster,scale_BSR_ω1ω2,m3_target,k3_target):
    """ Write out BSR to a NetCDF file """
    #  Save BSR
    write_variable_to_netcdf(
        BSR_ω1ω2_m1k1, OUTPUT_NC, "BSR_ω1ω2_m1k1",
        coord_info={ 0: {"name": "rank"},
                     1: {"name": "m", "values": m_BSR},
                     2: {"name": "k", "values": k_BSR}, },  mode="w" )
    write_variable_to_netcdf(
        BSR_ω1ω2_m2k2, OUTPUT_NC, "BSR_ω1ω2_m2k2",
        coord_info={ 0: {"name": "rank"}, 1: {"name": "m"}, 2: {"name": "k"}, },  mode="a" )
    write_variable_to_netcdf(
        BSR_ω1ω2_m3k3, OUTPUT_NC, "BSR_ω1ω2_m3k3",
        coord_info={ 0: {"name": "rank"}, 1: {"name": "m"}, 2: {"name": "k"}, },  mode="a" )

    #  Save groups_cluster
    groups = dict(groups_cluster)  # Convert defaultdict to dict

    cluster_ids = np.array(list(groups.keys()), dtype=np.int64)
    max_len = max(len(v) for v in groups.values())

    fill_value = -999
    arr = np.full((len(groups), max_len), fill_value, dtype=np.int64)

    for i, key in enumerate(cluster_ids):   #  Save groups_cluster in an array
        vals = groups[key]                  #  "fill_value" is the missing value
        arr[i, :len(vals)] = vals

    write_variable_to_netcdf(
        arr, OUTPUT_NC, "groups_cluster",
        coord_info={ 0: {"name": "cluster"}, 1: {"name": "member"}, },  mode="a" )

    #   Save parameters
    params = np.array([scale_BSR_ω1ω2, \
                       k3_target[0], k3_target[1], m3_target[0], m3_target[1]])
    write_variable_to_netcdf(
        params, OUTPUT_NC, "parameters", coord_info={ 0: {"name": "n"}, },  mode="a" )
#}}}1

#==========================================================================
#       Remaining functions

def check_consistency(a, b, txt_a, txt_b): #{{{1
    """ Check if ``a'' is almost identical to ``b'' """
    if not np.allclose(a, b):
        print(txt_a+' = ',a)
        print(txt_b+' = ',b)
        raise ValueError("np.allclose("+txt_a+", "+txt_b+") = False")
    return
#}}}1

def Map_Onto_k3m3(k12,m12,k3,m3,idx,arr): #{{{1
    """ Map arr in (k3,m3) space """
    im1,ik1,im2,ik2 = idx
    #
    ik3 = np.argmin( np.abs(k12 - k3) )
    im3 = np.argmin( np.abs(m12 - m3) )
    #
    arr_k3m3 = np.full( (len(m12),len(k12)), 0.)
    arr_k3m3[im3,ik3] = arr[im1,ik1,im2,ik2]
    #
    return arr_k3m3
#}}}1

def Write_Out_To_TextFile( #{{{1
        OUTPUT_TXT,
        m3_target, k3_target,
        m_BSR,k_BSR, m3_table,k3_table, im3_cond,ik3_cond, BSR_ω1ω2_cond,
        scale_BSR_ω1ω2, groups_cluster, list_m1k1m2k2, BSR_ω1ω2,
        numb_TopValues,
        ):
    text_out = []

    #  Limit_BSR_To_Range
    text_out.append(
          f"m3_target = ({np.min(m3_target): .2e}, {np.max(m3_target): .2e}); "+
          f" {np.min(m3_table[im3_cond]): .2e} <= m3[im3_cond] <= {np.max(m3_table[im3_cond]): .2e}")
    text_out.append(
          f"k3_target = ({np.min(k3_target): .2e}, {np.max(k3_target): .2e}); "+
          f" {np.min(k3_table[ik3_cond]): .2e} <= k3[ik3_cond] <= {np.max(k3_table[ik3_cond]): .2e}")

    numb_target_grids = \
            np.count_nonzero((k_BSR >= np.min(k3_target)) & (k_BSR <= np.max(k3_target))) \
          * np.count_nonzero((m_BSR >= np.min(m3_target)) & (m_BSR <= np.max(m3_target)))

    text_out.append(
          f"# of grids for (k3_target[0] <= k3 <= k3_target[1] &&"+
                          f"m3_target[0] <= m3 <= m3_target[1]) = {numb_target_grids:d}")
    text_out.append(f"numb_TopValues = {numb_TopValues:d}")
    text_out.append("")

    text_out.append("min/max of BSR_ω1ω2 in this (k3,m3) range = "+\
                f"{np.min(BSR_ω1ω2_cond):6.2e} / {np.max(BSR_ω1ω2_cond):6.2e} (not scaled)")
    text_out.append("                                          = "+\
                f"{np.min(BSR_ω1ω2_cond)/scale_BSR_ω1ω2:6.2e} / "+\
                f"{np.max(BSR_ω1ω2_cond)/scale_BSR_ω1ω2:6.2e} (scaled)")

    #  Top N values
    numb_tot = int(np.sum(BSR_ω1ω2_cond > 0))
    text_out.append(
        f"# of BSR values in the region = {numb_tot:d};  "
        f"# of chosen values = {numb_TopValues:d};  "
        f"Ratio = {numb_TopValues / numb_tot *100.: .1f} %")

    #  Clustering
    text_out.append("\n* Clustering")
    for label, indices in groups_cluster.items():
        text_out.append(f"  Cluster {label}: {indices} (Total {len(indices)} member(s))")

    text_out.append("\n* Each member (Index indicates rank)")
    text_out.append(f"    Δk={k_BSR[1]-k_BSR[0]: .2e}; Δm={m_BSR[1]-m_BSR[0]: .2e} m^-1")

    for label, indices in groups_cluster.items():
        for index in indices:
            im1,ik1,im2,ik2 = list_m1k1m2k2[index]

            text_out.append(
                    f"  Cluster {label}; Index {index:2}: "+
                    f"k1={k_BSR[ik1]: .2e}; "+f"m1={m_BSR[im1]: .2e}; "+
                    f"k2={k_BSR[ik2]: .2e}; "+f"m2={m_BSR[im2]: .2e}; "+
                    f"k3={k3_table[ik1,ik2]: .2e}; "+f"m3={m3_table[im1,im2]: .2e}; "+
                    f"BSR/scale = {BSR_ω1ω2[im1,ik1,im2,ik2]/scale_BSR_ω1ω2:.2f}")
        text_out.append("")

    #  Write out
    print('OUT >>> '+OUTPUT_TXT)
    with open(OUTPUT_TXT,'w') as fa:
        for text in text_out:
            fa.write(text)
            fa.write('\n')
            print(text)
#}}}1

def Draw_Figure_to_PDF(OUTPUT_PDF, #{{{1
                       m_BSR,k_BSR,groups_cluster,
                       BSR_ω1ω2_m1k1,BSR_ω1ω2_m2k2,BSR_ω1ω2_m3k3,scale_BSR_ω1ω2,
                       m_range,k_range,ω1,ω2,ω3,m3_target,k3_target,
                       ttl_periods=None,
                       ):

    print('OUT >>> '+OUTPUT_PDF)
    with PdfPages(OUTPUT_PDF) as pdf:

        nrows, ncols = 3, 3
        xcd,ycd, xax_range,yax_range = k_BSR,m_BSR, k_range,m_range
        max_rank = 20       #   Draw figures only when index < max_rank

        #-------------------------------------------------------------------------------------
        #   Real bispectrum in (k1,m1), (k2,m2), and (k3,m2) axes
        #-------------------------------------------------------------------------------------

        clab = 'Real Bispectrum'

        fig, axes, iplot = new_page(nrows, ncols)

        v0, cmap = 0.8, 'RdBu_r'
        norm = colors.SymLogNorm(linthresh=1e0, linscale=1e0, vmin=-v0, vmax=v0)

        for label, indices in groups_cluster.items():
            for index in indices:
                if index < max_rank:
                    for ttl0, ω, xlab,ylab, arr in zip(
                            ttl_periods,
                            (ω1,ω2,ω3,),
                            ("$k_1$","$k_2$","$k_3$"),
                            ("$m_1$","$m_2$","$m_3$"),
                            (BSR_ω1ω2_m1k1[index,:,:]/scale_BSR_ω1ω2,
                             BSR_ω1ω2_m2k2[index,:,:]/scale_BSR_ω1ω2,
                             BSR_ω1ω2_m3k3[index,:,:]/scale_BSR_ω1ω2),
                            ):

                        ttl = f"Cluster {label}, Index {index}, "+ttl0+" days"

                        ax = axes[iplot]
                        pcm = ax.pcolormesh(xcd,ycd,arr, cmap=cmap, shading='auto', norm=norm)
                        fig.colorbar(pcm, ax=ax, label=clab, location='bottom')
                        ax = set_ax(ax, ttl, xax_range, yax_range, xlab,ylab)
                        overplot_DispersionRelation(np.abs(ω), xcd, ax, line_width=0.5)

                        #   Add a rectangle to show k3-m3 range
                        if ttl0 == "60":
                            rect = patches.Rectangle(
                                   (np.min(k3_target), np.min(m3_target)),
                                   np.max(k3_target)-np.min(k3_target),
                                   np.max(m3_target)-np.min(m3_target),
                                   linestyle=':',linewidth=1,
                                   edgecolor='black', facecolor='none' )
                            axes[iplot].add_patch(rect)

                        fig,pdf,axes, iplot = manage_page(iplot,nrows,ncols,fig,pdf,axes)

        plt.tight_layout()
        pdf.savefig()
        plt.close(fig)
#}}}1

#==========================================================================
#           No Longer Used

def Map_arr_in_km_space_using_list(arr,m_BSR,k_BSR, list_TopNValues, #{{{1
                                   m3_table,k3_table, m3_target,k3_target):
    """
      Superseded by "Substitute_Listed_Values_in_km_space"

      Map given variable to (k1,m1), (k2,m2), and (k3,m3) space using list_TopNValues

     Worked parameters:
      cluster_eps             = 2.5   #    Larger values group less similar shapes together
      cluster_threshold_ratio = 0.2   #    Threshold for values to be considered
      cluster_crop_and_pad    = False #    If true, patterns are cropped and padded
    """

    #   Save BSR in the list

    ntop = len(list_TopNValues)

    arr_ω1ω2_m1k1, arr_ω1ω2_m2k2, arr_ω1ω2_m3k3 = np.zeros((3,ntop,len(m_BSR),len(k_BSR)))

    for i_rank, idx in enumerate(list_TopNValues):
        im1,ik1,im2,ik2 = idx

        #  Error check
        if np.logical_or( np.logical_or( m3_table[im1,im2] < np.min(m3_target),
                                         m3_table[im1,im2] > np.max(m3_target)),
                          np.logical_or( k3_table[ik1,ik2] < np.min(k3_target),
                                         k3_table[ik1,ik2] > np.max(k3_target)) ):
            raise ValueError(
                    "m3_table[im1,im2] or k3_table[ik1,ik2] is out of range!\n"
                    f"   (im1,im2)=({im1:d},{im2:d}),  m3_table[im1,im2]={m3_table[im1,im2]:.2e}\n"
                    f"   (ik1,ik2)=({ik1:d},{ik2:d}),  k3_table[ik1,ik2]={k3_table[ik1,ik2]:.2e}\n"
                    f"   {np.min(m3_target):.2e} <= m3_target <= {np.max(m3_target):.2e}\n"
                    f"   {np.min(k3_target):.2e} <= k3_target <= {np.max(k3_target):.2e}\n"
                    f"   BSR there = {arr[im1,ik1,im2,ik2]:.2e}"
                             )

        #  Save
        arr_ω1ω2_m1k1[i_rank,:,:] = arr[:,:,im2,ik2]
        arr_ω1ω2_m2k2[i_rank,:,:] = arr[im1,ik1,:,:]
        arr_ω1ω2_m3k3[i_rank,:,:] = Map_Onto_k3m3(k_BSR,m_BSR,
                                                  k3_table[ik1,ik2],m3_table[im1,im2],
                                                  idx,arr)

    return arr_ω1ω2_m1k1, arr_ω1ω2_m2k2, arr_ω1ω2_m3k3
#}}}1
