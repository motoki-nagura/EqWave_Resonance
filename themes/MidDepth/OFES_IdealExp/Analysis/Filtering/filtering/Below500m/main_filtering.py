"""
    Apply filtering in ω-k-m space to zonal velocity time series below 500 m and
      plot the variance of filtered time series in x-z space
"""

#============================================================================================#
#           Imports
#============================================================================================#
import numpy as np
import xarray as xr
from scipy import signal
from scipy.interpolate import interp1d

import string
import time as time_measure
from rich import traceback
import sys, pathlib, os

import matplotlib.pyplot as plt
from matplotlib import ticker
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import AutoMinorLocator
import matplotlib.patches as patches

from common.new_page import (new_page, )

from themes.MidDepth.OFES_IdealExp.Tools.common.Fixed_Parameters import (
        equatorial_beta, degrees_to_meters, )
from themes.MidDepth.OFES_IdealExp.Tools.common.Overplot_DispersionRelation import (
        Overplot_DispersionRelation,)
from common.fft_spectrum_3d import fft_spectrum_3d
from common.sub_NetCDF_IO import (write_variable_to_netcdf, read_variable_from_netcdf, )
from themes.MidDepth.OFES_IdealExp.Tools.common.Bandpass_Frequencies import (
        Bandpass_Frequencies, )
from themes.MidDepth.OFES_IdealExp.Analysis.Filtering.filtering.sub_FiltMultiDim import (
        MakeFilter, ApplyFilter, )

#============================================================================================#
#           Parameters    {{{1
#============================================================================================#
π = np.pi

β = equatorial_beta()         #  Meridional gradient of the Coriolis coefficient [m^-1 s^-1]
deg2met = degrees_to_meters() #  1 deg in meters

day_in_sec = 86400.           #  1 day in seconds
#}}}1

#============================================================================================#
#           Class   {{{1
#============================================================================================#
class InputParams:
    def __init__(self, grid_name, exp_name, period_name, time_width,
                       zstrt, zend, lstrt, lend, N2, T_target, InOutDirectory
                 ):
        self.grid_name     = grid_name
        self.exp_name      = exp_name
        self.period_name   = period_name
        self.time_width    = time_width
        self.zstrt         = zstrt
        self.zend          = zend
        self.lstrt         = lstrt
        self.lend          = lend
        self.N2            = N2
        self.T_target      = T_target
        self.DirectoryName = InOutDirectory
#}}}1

#============================================================================================#
#           Functions   {{{1
#============================================================================================#
#
#   main ─┬─── compute_and_write
#         │        ├─── get_filter_omega
#         │        │       └─── Bandpass_Filtering
#         │        ├─── get_filter_km_60days (or get_filter_km_90days)
#         │        ├─── MakeFilter
#         │        ├─── ApplyFilter
#         │        ├─── compute_psd --- 
#         │        │       └─── fft_spectrum_3d
#         │        └─── write_variable_to_netcdf
#         │
#         ├─── read_variable_from_netcdf
#         │
#         └─── make_figures

def get_filter_omega(T_target): #{{{2
    cpd_lowcut, cpd_highcut = Bandpass_Frequencies(T_target)
    ω_lowcut  = cpd_lowcut  *2.*π /86400.        #   cpd --> rad s^-1
    ω_highcut = cpd_highcut *2.*π /86400.
    ω_tap     =  0.5 *(ω_highcut -  ω_lowcut)
    print('ω_lowcut, ω_highcut = ',ω_lowcut, ω_highcut)
    print('2π/ω_lowcut/86400, 2π/ω_highcut/86400 = ',2*π/ω_lowcut/86400, 2*π/ω_highcut/86400)
    return ω_lowcut,ω_highcut,ω_tap
#}}}2

def get_filter_km(T_target): #{{{2
    #--------------------
    #  k and m [rad m^-1]

    if   T_target == 60:                #  60 days
        k_lowcut, k_highcut, k_tap, m_lowcut, m_highcut, m_tap, l_cut = get_filter_km_60days()

    elif T_target == 90:                #  90 days
        k_lowcut, k_highcut, k_tap, m_lowcut, m_highcut, m_tap, l_cut = get_filter_km_90days()

    elif T_target in [140./3., 160./3., 200./3., 220./3.]: #  Sensitivity experiments
        k_lowcut, k_highcut, k_tap, m_lowcut, m_highcut, m_tap, l_cut = get_filter_km_60days()

    else:
        raise ValueError(f"Invalid input of T_target (now ={T_target:.2f}) !")

    return k_lowcut, k_highcut, k_tap, m_lowcut, m_highcut, m_tap, l_cut
#}}}2

def get_filter_km_60days(): #{{{2
    k_lowcut, k_highcut, m_lowcut, m_highcut = np.zeros((4,4))
    k_tap, m_tap = np.zeros((2,4))

    k1l, k1r, m1l, m1u = -2.0e-5, -0.2e-5, -0.6e-2, 0.6e-2   #  Rossby
    k1_tap, m1_tap     =  4.5e-6,  0.4e-2

    k0l, k0r, m0l, m0u =  0.3e-5,  1.5e-5,  0.9e-2, 1.9e-2   #  Off resonance / Kelvin
    k0_tap, m0_tap     =  3.0e-6,  0.5e-3

    k_lowcut[0], k_highcut[0], m_lowcut[0], m_highcut[0] =  k1l,  k1r,  m1l,  m1u   #  Rossby
    k_lowcut[1], k_highcut[1], m_lowcut[1], m_highcut[1] = -k0r, -k0l,  m0l,  m0u   #  k < 0,  m > 0
    k_lowcut[2], k_highcut[2], m_lowcut[2], m_highcut[2] =  k0l,  k0r,  m0l,  m0u   #  k > 0,  m > 0
    k_lowcut[3], k_highcut[3], m_lowcut[3], m_highcut[3] =  k0l,  k0r, -m0u, -m0l   #  k > 0,  m < 0

    k_tap[0], m_tap[0] = k1_tap, m1_tap  #  Rossby
    k_tap[1], m_tap[1] = k0_tap, m0_tap  #  Off resonance / Kelvin
    k_tap[2], m_tap[2] = k0_tap, m0_tap
    k_tap[3], m_tap[3] = k0_tap, m0_tap

    l_cut = 112    #   Edge effect is non negligible from 0:lcut and -lcut:
                   #    See Analysis/Filtering/filtering/check/out_Rossby-waves_60dy.txt

    return k_lowcut, k_highcut, k_tap, m_lowcut, m_highcut, m_tap, l_cut
#}}}2

def get_filter_km_90days(): #{{{2
    k_lowcut, k_highcut, m_lowcut, m_highcut = np.zeros((4,4))
    k_tap, m_tap = np.zeros((2,4))

    k1l, k1r, m1l, m1u = -2.0e-5, -0.2e-5, -0.6e-2, 0.6e-2   #  Rossby
    k1_tap, m1_tap     =  4.5e-6,  0.4e-2

    k0l, k0r, m0l, m0u =  0.3e-5,  1.5e-5,  0.8e-2, 2.5e-2   #  Kelvin
    k0_tap, m0_tap     =  3.0e-6,  0.5e-3

    koffl, koffr, moffl, moffu = -1.5e-5, -0.3e-5,  0.8e-2, 2.0e-2  #  Off resonance

    k_lowcut[0], k_highcut[0], m_lowcut[0], m_highcut[0] = k1l,   k1r,   m1l,   m1u   #  Rossby
    k_lowcut[1], k_highcut[1], m_lowcut[1], m_highcut[1] = koffl, koffr, moffl, moffu #  Off resonant
    k_lowcut[2], k_highcut[2], m_lowcut[2], m_highcut[2] = k0l,   k0r,   m0l,   m0u   #  k > 0,  m > 0
    k_lowcut[3], k_highcut[3], m_lowcut[3], m_highcut[3] = k0l,   k0r,  -m0u,  -m0l   #  k > 0,  m < 0

    k_tap[0], m_tap[0] = k1_tap, m1_tap  #  Rossby
    k_tap[1], m_tap[1] = k0_tap, m0_tap  #  Off resonance / Kelvin
    k_tap[2], m_tap[2] = k0_tap, m0_tap
    k_tap[3], m_tap[3] = k0_tap, m0_tap

    l_cut = 168    #   Edge effect is non negligible from 0:lcut and -lcut:
                   #    See Analysis/Filtering/filtering/check/out_Rossby-waves_90dy.txt

    return k_lowcut, k_highcut, k_tap, m_lowcut, m_highcut, m_tap, l_cut
#}}}2

def compute_psd(var,x_in,z_in,t_in, time_width): #{{{2
    nt, nz, nx = len(t_in), len(z_in), len(x_in)

    hanning_1d = np.hanning(time_width)
    hanning_window = np.tile(hanning_1d[:, None, None], (1, nz,nx))

    #-------
    l0, l1, icnt = 0, time_width, 0.   #  initial l0, l1, and a counter

    psd = np.zeros((time_width,nz,nx))      #  [{ft},{fz},{fx}]

    while l1 < nt:

        var_seg =  var[l0:l1,:,:] * hanning_window
        t_seg   = t_in[l0:l1]

        tmp_psd, fx,fz,ft = fft_spectrum_3d(var_seg, x_in,z_in,t_seg)
        tmp_psd = tmp_psd *(8./3.)          #  Correction for the Hanning window in time

        psd = psd + tmp_psd

        icnt = icnt + 1.

        l0 += int(time_width /2)   #      advance by time_width/2
        l1 += int(time_width /2)

    psd = psd /icnt
    edof = icnt *2.

    return psd,fx,fz,ft
#}}}2

def f2T(x): #{{{2
    return 1./x
def T2f(x):
    return 1./x
#}}}2

def compute_and_write(input_params,): #{{{2

    time_elapsed_start = time_measure.perf_counter()

    #============================================================================================#
    #           Unpack input_params
    #============================================================================================#
    grid_name   = input_params.grid_name
    exp_name    = input_params.exp_name
    period_name = input_params.period_name
    time_width  = input_params.time_width
    zstrt       = input_params.zstrt
    zend        = input_params.zend
    lstrt       = input_params.lstrt
    lend        = input_params.lend
    #N2          = input_params.N2
    T_target    = input_params.T_target
    OUT_DIR     = input_params.DirectoryName

    #============================================================================================#
    #           Set ω, k, m ranges
    #============================================================================================#

    # ---------------  Set ω, k, m ranges  ------------------

    ω_lowcut, ω_highcut, ω_tap                                    = get_filter_omega(T_target)
    k_lowcut, k_highcut, k_tap, m_lowcut, m_highcut, m_tap, l_cut = get_filter_km(   T_target)

    #================================================================================================
    #                  Read In Time Series
    #================================================================================================
    print('*  Reading input ...')

    MY_WORK_ROOT = str( pathlib.Path(os.environ["MY_WORK_ROOT"]) )

    infile = MY_WORK_ROOT+"/data/ofes_exp/results/"+grid_name+"/compiled/"+\
             "u_eq_"+grid_name+"_"+exp_name+"_"+period_name+".nc"

    print('   IN <<< '+infile)

    var  = xr.open_dataset(infile)['u'].isel(lat=0,time=slice(lstrt,lend)).sel(lev=slice(zstrt,zend))

    lon  = np.array(var.lon)
    lev  = np.array(var.lev)
    var  = np.squeeze( np.array(var) )

    #
    time = np.arange(0., len(var[:,0,0]))   #  time in days

    #  Remove NaN
    i = np.argwhere(np.logical_not(np.isnan(var[0,0,:])))[:,0]
    var, lon = var[:,:,i], lon[i]

    #
    lon_met = lon *deg2met                 #  X coordinate, in meters

    lev_flipped = np.flip(lev *(-1.))      #  Z coordinate is positive upward
                                           #  Flipped in the vertical direction such that
                                           #     array is from -5000 m to 0 m

    var_flipped = np.flip(var, axis=1)

    time_seconds = time *day_in_sec        #  T coordinate, in seconds

    #
    time_elapsed_end = time_measure.perf_counter()
    print(f"      Elapsed time: {time_elapsed_end - time_elapsed_start:.6f} seconds")

    #================================================================================================
    #                  Bandpass filtering
    #================================================================================================
    print('*  Bandpass filtering ...')

    x,z,t = lon_met, lev_flipped, time_seconds

    Weight = np.zeros((len(k_lowcut), len(t),len(z),len(x)))

    for i in range(len(k_lowcut)):

        ω_cent  = 0.5 *(ω_lowcut    + ω_highcut   )
        ks_cent = 0.5 *(k_lowcut[i] + k_highcut[i])
        ms_cent = 0.5 *(m_lowcut[i] + m_highcut[i])

        ω_hwid  = 0.5 *(ω_highcut    - ω_lowcut   )
        ks_hwid = 0.5 *(k_highcut[i] - k_lowcut[i])
        ms_hwid = 0.5 *(m_highcut[i] - m_lowcut[i])

        ω_tap =  ω_hwid *0.5
        k_tap = ks_hwid *0.5
        m_tap = ms_hwid *0.5

        is_negative = (ω_hwid < 0) or np.any(ks_hwid < 0) or np.any(ms_hwid < 0)
        if is_negative:
            raise ValueError("### Invalid input of k_lowcut, k_highcut, m_lowcut or m_highcut. Stopped!")

        #

        ω0,ω1,ω2 =  ω_cent, ω_hwid, ω_tap
        k0,k1,k2 = ks_cent,ks_hwid, k_tap
        m0,m1,m2 = ms_cent,ms_hwid, m_tap

        tmp_weight_pos, ft_w, fz_w, fx_w = MakeFilter(x,z,t,  ω0,ω1,ω2,  k0,k1,k2,  m0,m1,m2)  # (ω,k,m)
        tmp_weight_neg, ft_w, fz_w, fx_w = MakeFilter(x,z,t, -ω0,ω1,ω2, -k0,k1,k2, -m0,m1,m2)  # (-ω,-k,-m)

        Weight[i,:,:,:] += tmp_weight_pos
        Weight[i,:,:,:] += tmp_weight_neg

    Weight[ Weight >= 1. ] = 1.   #   to avoid duplication

    #
    #idx_t = np.argsort(ft_w)
    #idx_z = np.argsort(fz_w)
    #idx_x = np.argsort(fx_w)

    #Weight_sort = np.copy(Weight)
    #Weight_sort = Weight_sort[:,idx_t,:,:]
    #Weight_sort = Weight_sort[:,:,idx_z,:]
    #Weight_sort = Weight_sort[:,:,:,idx_x]

    #ft_w_sort = ft_w[idx_t]
    #fz_w_sort = fz_w[idx_z]
    #fx_w_sort = fx_w[idx_x]

    #
    var_flipped_filtered = np.zeros((len(k_lowcut), len(t),len(z),len(x)))

    for i in range(len(k_lowcut)):
        var_flipped_filtered[i,:,:,:] = ApplyFilter(var_flipped, Weight[i,:,:,:])

    #
    time_elapsed_end = time_measure.perf_counter()
    print(f"      Elapsed time: {time_elapsed_end - time_elapsed_start:.6f} seconds")

    #========================================================================================
    #          Variance of filtered time series
    #========================================================================================
    print('*  Variance of filtered time series ...')

    var_filtered = np.flip(var_flipped_filtered, axis=2) #  flip in z

    print('var_filtered.shape = ',var_filtered.shape)
    var_filtered = var_filtered[:, l_cut:-l_cut, :,:]    #  Remove the edge effect in time
    print('var_filtered.shape = ',var_filtered.shape)

    var_filtered_variance = np.var(var_filtered, axis=1) #  about time

    #========================================================================================
    #          Power Spectrum Density in ω-k-m space
    #========================================================================================
    print('*  Computing spectrum ...')

    #nt,nz,nx = np.shape(var_flipped)     #   # of time steps in the data

    #-------
    var_flipped          = signal.detrend(var_flipped,          type='constant', axis=0)
    var_flipped_filtered = signal.detrend(var_flipped_filtered, type='constant', axis=1)

    #-------
    psd, fx_p,fz_p,ft_p = compute_psd(var_flipped, lon_met, lev_flipped, time, time_width)

    n1,n2,n3 = np.shape(psd)
    psd_filtered = np.zeros((len(k_lowcut), n1,n2,n3))

    for i in range(len(k_lowcut)):
        psd_filtered[i,:,:,:], fx_p,fz_p,ft_p = \
            compute_psd(var_flipped_filtered[i,:,:,:], lon_met, lev_flipped, time, time_width)

    #-------
    ft_p = 2.*π *ft_p /day_in_sec  #  Frequency [cycle per day] to [rad second^-1]
    fx_p = 2.*π *fx_p              #  Zonal wavenumber [cycle per meter] to [rad meter^-1]
    fz_p = 2.*π *fz_p              #  Vertical wavenumber [cycle per meter] to [rad meter^-1]

    #
    time_elapsed_end = time_measure.perf_counter()
    print(f"      Elapsed time: {time_elapsed_end - time_elapsed_start:.6f} seconds")

    #========================================================================================
    #          Write out to NetCDF files
    #========================================================================================
    print('*  Writing ...')

    write_variable_to_netcdf(
        var_filtered_variance,
        OUT_DIR+"/tmp_variance.nc",
        "var_filtered_variance",
        coord_info={
            0: {"name": "filter_type"},
            1: {"name": "lev", "values": lev},
            2: {"name": "lon", "values": lon},
        },
        mode="w"
    )

    write_variable_to_netcdf(
        psd,
        OUT_DIR+"/tmp_psd.nc",
        "psd",
        coord_info={
            0: {"name": "ft", "values": ft_p},
            1: {"name": "fz", "values": fz_p},
            2: {"name": "fx", "values": fx_p},
        },
        mode="w"
    )

    write_variable_to_netcdf(
        psd_filtered,
        OUT_DIR+"/tmp_psd_filtered.nc",
        "psd_filtered",
        coord_info={
            0: {"name": "filter_type"},
            1: {"name": "ft", "values": ft_p},
            2: {"name": "fz", "values": fz_p},
            3: {"name": "fx", "values": fx_p},
        },
        mode="w"
    )

    time_elapsed_end = time_measure.perf_counter()
    print(f"      Elapsed time: {time_elapsed_end - time_elapsed_start:.6f} seconds")
#}}}2

def make_figures( #{{{2
                 var_filtered_variance,lev,lon,
                 psd,ft_p,fz_p,fx_p, psd_filtered,
                 input_params):

    #========================================================================================

    T_target, exp_name, N2 = \
         input_params.T_target, input_params.exp_name, input_params.N2

    print('N2 = ',N2)

    outfile_fig = 'fig_Exp'+exp_name+'_T'+f'{int(T_target):d}'+'dy.pdf'

    #
    ω_lowcut, ω_highcut, _ = get_filter_omega(T_target)
    ω_cent = 0.5 *(ω_lowcut + ω_highcut)

    k_lowcut, k_highcut, _, m_lowcut, m_highcut, _, _ = get_filter_km(T_target)

    #========================================================================================
    #              Figure
    #========================================================================================
    print('*  Making figures...')

    #--------------------------------------------------------------------
    krange, mrange = [-2e-5,2e-5], [-4e-2,4e-2]  #  range in k and m to plot

    ft,fz,fx = ft_p,fz_p,fx_p
    ω_target = ω_cent
    print("2*π/ω_target/86400. = ", 2*π/ω_target/86400.)

    lon_Ticks = np.arange(0,61,20)
    lon_TickNames = [f"{int(x)}\N{degree sign}" for x in lon_Ticks]

    txt_labels = [f"({c})" for c in string.ascii_lowercase]

    txt_psd_unit = r'[cm$^2$ s$^{-2}$ (day$^{-1}$)$^{-1}$ (m$^{-1}$)$^{-2}$]'

    with PdfPages(outfile_fig) as pdf:

        #---------------------------------------------------------------------
        #          Variance of filtered time series
        #---------------------------------------------------------------------

        fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(8,5.5))
        axes = axes.flatten()     #   row major

        lvls1, lvls2, lvls3, lvls4, ctick1, ctick2, ctick3, ctick4 = \
            make_lvls_vars(T_target, var_filtered_variance)

        for iplot, (lvls, cbar_ticks) in enumerate(
            zip((lvls1, lvls2, lvls3, lvls4),
                (ctick1, ctick2, ctick3, ctick4))
             ):

            ax = axes[iplot]

            ax.set_title(txt_labels[iplot]+f"  Passband {iplot+1:d}")
            cs=ax.contourf(lon,lev,var_filtered_variance[iplot,:,:],lvls,
                           cmap='jet', extend='both')
            cbar = plt.colorbar(cs,ax=ax)
            cbar.set_label(r'Variance (cm$^2$ s$^{-2}$)')
            cbar.set_ticks(cbar_ticks)
            ax.invert_yaxis()

        for ax in axes:
            ax.set(xlabel='Longitude', ylabel='Depth (m)')
            ax.set_xticks(lon_Ticks)
            ax.set_xticklabels(lon_TickNames)
            ax.xaxis.set_minor_locator(AutoMinorLocator(4))
            ax.yaxis.set_ticks(np.arange(0., 5001., 1000.))
            ax.yaxis.set_minor_locator(AutoMinorLocator(5))
            ax.tick_params(right=True,top=True,which='both')

        plt.tight_layout()
        pdf.savefig()
        plt.close(fig)

        #---------------------------------------------------------------------
        #            Spectrum on k-m space
        #---------------------------------------------------------------------

        fig, axes, iplot = new_page(nrows=3, ncols=2, figsize=(8,10))

        for iloop in range(len(k_lowcut)+1):

            ax = axes[iplot]

            if iloop == len(k_lowcut):
                var = psd
            else:
                var = psd_filtered[iloop,:,:,:]

            #    Interpolate PSD onto the target frequency

            interp_fun = interp1d(ft,var, axis=0, kind='linear', fill_value="extrapolate")
            var_intp = interp_fun( ω_target )

            print("iloop, np.max(np.abs(var_intp)) = ", iloop, np.max(np.abs(var_intp)))

            #     title
            T = int( np.round( 2.*π/ω_target /86400. ) )

            title = txt_labels[iplot]+f'  {T:d} days'
            ax.text( 0.04, 0.96, title,
                        transform=ax.transAxes, ha="left", va="top",
                        bbox=dict(facecolor="white", edgecolor="black")
            )

            #     contour
            lvls = 10**np.arange(6.0, 10.6, 0.5)
            cs = ax.contourf(fx,fz,var_intp,lvls, cmap='jet',extend='both',
                             locator=ticker.LogLocator())

            #     color bar
            cbar = plt.colorbar(cs,ax=ax,location='right',pad=0.17)
            cbar.ax.tick_params(labelsize=9)
            cbar.set_label(txt_psd_unit,fontsize=9)

            #    axes
            ax = define_figure_axes(ax, krange, mrange)

            #    Dispersion relations
            Overplot_DispersionRelation(np.abs(ω_cent), fx, ax, N2=N2)

            #   Mark the target ω and m
            if iloop == len(k_lowcut):
               ax = Mark_PassBands(ax, k_lowcut,k_highcut, m_lowcut, m_highcut)
            else:
               i = iloop
               ax = Mark_PassBands(ax, k_lowcut[i],k_highcut[i], m_lowcut[i], m_highcut[i])

            iplot += 1

        plt.tight_layout()
        pdf.savefig()
        plt.close(fig)
        del(iplot)
#}}}2

def make_lvls_vars(T_target, var_filtered_variance):  #{{{2
    # Levels for Variance
    if  T_target == 60:       #  60 days
        lvls1, lvls2, lvls3, lvls4 = \
            np.arange(0., 0.51, 0.05), np.arange(0., 0.31, 0.05), \
            np.arange(0., 0.26, 0.05), np.arange(0., 0.61, 0.10)

        ctick1, ctick2, ctick3, ctick4 = \
            np.arange(0., 0.51, 0.1 ), np.arange(0., 0.31, 0.1 ), \
            np.arange(0., 0.24, 0.1 ), np.arange(0., 0.61, 0.1 )

    elif T_target == 90:      #  90 days
        lvls1, lvls2, lvls3, lvls4 = \
            np.arange(0., 3.01, 0.50), np.arange(0., 0.31, 0.05), \
            np.arange(0., 0.31, 0.05), np.arange(0., 0.51, 0.10)

        ctick1, ctick2, ctick3, ctick4 = \
            np.arange(0., 0.51, 0.1 ), np.arange(0., 0.31, 0.1 ), \
            np.arange(0., 0.31, 0.1 ), np.arange(0., 0.51, 0.1 )

    else:
        for iplot in range(4):
            max_val = np.max(var_filtered_variance[iplot,:,:])
            #
            if  max_val > 0.1:
                max_val = round(max_val, 1)
            elif max_val > 0.01:
                max_val = round(max_val, 2)
            #
            tmp_lvls  = np.linspace(0., max_val, 11)
            tmp_ctick = np.linspace(0., max_val, 4)
            #
            if   iplot == 0:
                lvls1,ctick1 = tmp_lvls, tmp_ctick
            elif iplot == 1:
                lvls2,ctick2 = tmp_lvls, tmp_ctick
            elif iplot == 2:
                lvls3,ctick3 = tmp_lvls, tmp_ctick
            elif iplot == 3:
                lvls4,ctick4 = tmp_lvls, tmp_ctick

    return  lvls1, lvls2, lvls3, lvls4, ctick1, ctick2, ctick3, ctick4
#}}}2

def define_figure_axes(ax, krange, mrange): #{{{2

    #     X and Y axes
    ax.axhline(y=0., c='k',ls='-',lw=0.5)
    ax.axvline(x=0., c='k',ls='-',lw=0.5)

    ax.set_xlim(krange)
    ax.set_ylim(mrange)

    ax.set_xlabel(r'$k$ (m$^{-1}$)',fontsize=10)
    ax.set_ylabel(r'$m$ (m$^{-1}$)',fontsize=10)

    ax.tick_params(axis='both', which='major', labelsize=9)
    ax.yaxis.set_ticks(np.arange(-0.04, 0.041, 0.02))
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.yaxis.set_minor_locator(AutoMinorLocator(5))

    ##     Secondary axes
    #xax2 = ax.secondary_xaxis('top', functions=(f2T,T2f))    #  Secondary axis on X
    #xax2.set_xlabel('Zonal Wavelength (km)',fontsize=10)
    #xax2.tick_params(axis='both', which='major', labelsize=9)

    #yax2 = ax.secondary_yaxis('right', functions=(f2T,T2f))  #  Secondary axis on Y
    #yax2.set_ylabel('Vertical Wavelength (km)',fontsize=10)
    #yax2.tick_params(axis='both', which='major', labelsize=9)

    #tmp = [-0.4, -0.75, 0.75, 0.4]
    #tick_labels = ['-400', '-750', '750', '400']
    #tmp1 = np.array(tmp)*1e6 /(2.*π)
    #xax2.xaxis.set_ticks(tmp1)
    #xax2.xaxis.set_ticklabels(tick_labels)

    #tmp = [-250, -1000, 1000, 250]
    #tick_labels = ['-0.25', '-1', '1', '0.25']
    #tmp1 = np.array(tmp) /(2.*π)
    #yax2.yaxis.set_ticks(tmp1)
    #yax2.yaxis.set_ticklabels(tick_labels)

    return ax
#}}}2

def Mark_PassBands(ax,k_lowcut,k_highcut, m_lowcut, m_highcut): #{{{2
    k_lowcut  = np.atleast_1d(k_lowcut)
    k_highcut = np.atleast_1d(k_highcut)
    m_lowcut  = np.atleast_1d(m_lowcut)
    m_highcut = np.atleast_1d(m_highcut)

    for k0,k1,m0,m1 in zip(k_lowcut, k_highcut, m_lowcut, m_highcut):
        rect = patches.Rectangle(
            (k0,m0),     # bottom-left corner
            (k1-k0),     # width
            (m1-m0),     # height
            linewidth=1,
            edgecolor='red',
            facecolor='none'
        )
        ax.add_patch(rect)
    return ax
#}}}2

#}}}1

#======================================================================
#             Main
#======================================================================
def main(input_params, Computation=True): #{{{1

    #   Compute and write out to temporary files

    if Computation:
        compute_and_write(input_params,)

    #   Read in from temporary output files

    var_filtered_variance, coords = read_variable_from_netcdf(
        input_params.DirectoryName+"/tmp_variance.nc",
        "var_filtered_variance",
        return_coords=True
    )
    lev, lon = coords["lev"], coords["lon"]
    del(coords)

    psd, coords = read_variable_from_netcdf(
        input_params.DirectoryName+"/tmp_psd.nc",
        "psd",
        return_coords=True
    )
    ft_p,fz_p,fx_p = coords["ft"], coords["fz"], coords["fx"]

    psd_filtered = read_variable_from_netcdf(
        input_params.DirectoryName+"/tmp_psd_filtered.nc",
        "psd_filtered",
    )

    #    Make figures

    make_figures(var_filtered_variance,lev,lon,
                 psd,ft_p,fz_p,fx_p,
                 psd_filtered,
                 input_params,
                 )
#}}}1

#======================================================================
#            Entry Point
#======================================================================

if __name__ == "__main__":
    traceback.install()

    current_dir = os.getcwd()                      #  Get current directory path
    OUT_DIR = current_dir.replace("code", "data")  #  Replace "code" with "data"

    #T_target = 60.  #  Target period [day]; 60 or 90
    T_target = 90.

    input_params = \
        InputParams(
           "box_bounded",  #  grid_name
           "180dy+90dy",   #  exp_name
           "2000-2029",    #  period_name
            1024,          #  time_width  (# of samples in the time direction)
            500.,          #  zstrt  (Depth range)
            5000.,         #  zend
            -4096,         #  lstrt  (Last 4096 days)
            None,          #  lend
            5e-6,          #  N2 (Buoyancy frequency squared [s^-2]; for the dispersion relation)
            T_target,      #  Target period
            OUT_DIR,       #  output directory
                  )

    #---
    sys.exit(main(input_params, Computation=True))
