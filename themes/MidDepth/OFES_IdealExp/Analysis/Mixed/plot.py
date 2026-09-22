"""
    Plot 1) spectrum in ω, 2) spectrum in k and m, and 3) variance of filtered time series
      for 140- and 220-day periods (or 160 and 200 days)
"""
from rich import traceback
traceback.install()

import numpy as np
import xarray as xr
import string, sys, pathlib, os
from scipy.interpolate import interp1d
from scipy import signal

import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from matplotlib import ticker
import matplotlib.patches as patches

from common.new_page import new_page
from common.sub_NetCDF_IO import read_variable_from_netcdf
from common.spec_err import spec_err
from themes.MidDepth.OFES_IdealExp.Tools.common.Overplot_DispersionRelation import (
        Overplot_DispersionRelation, )
from themes.MidDepth.OFES_IdealExp.Tools.common.EqWaves_Rays import (
        ray_K_LongRossby, draw_arrows_mid_head, )
from themes.MidDepth.OFES_IdealExp.Analysis.Filtering.filtering.Below500m.main_filtering import(
        InputParams, get_filter_omega, get_filter_km, )

#=================================================================
#=================================================================
def Averaged_Spectrum(infile,exp_name,MarkPeriods,alpha_spec, ax, ybot=1e-10, yerr=1.): #{{{
    #===================================================================================
    #                    Parameters
    #===================================================================================
    z0, z1 = 500., 5000.              #   Depth range of the analysis

    l0, l1 = -4096, None              #   Last 4096 days
    #l0, l1 = 6000, 10096               #   Day 6000 ~ 10095
    #l0, l1 = -3650, None              #   Read the last "l0" timesteps

    nperseg = 1024                    #  segment length

    #===================================================================================
    #                     Compute
    #===================================================================================
    Pxx,freq, edof, relerr_low, relerr_high = \
            Compute_Averaged_Spectrum(infile,l0,l1,z0,z1,nperseg,alpha_spec)

    #===================================================================================
    #                    Plots
    #===================================================================================
    #ax.set_title(panel_title)
    #ax.set_title(f'Average over longitudes and\n'+\
    #                     f'{int(z0):4d}-{int(z1):4d} m depths', fontsize=10)
    ax.plot(freq, Pxx, 'k-', label=exp_name)

    ax.set(xlabel='Frequency (cpd)', \
                   ylabel=r'PSD (cm$^2$ s$^{-2}$ cpd$^{-1}$)')

    #   Axis
    ax.tick_params(right=True,which='both')

    ax.set_xlim([1./500.,1./10.])
    ax.set_xscale('log')
    xax2 = ax.secondary_xaxis('top', functions=(f2T,T2f))
    xax2.set_xlabel('Period (day)')
    #
    xax2.xaxis.set_ticks(MarkPeriods)
    xax2.xaxis.set_ticklabels([f"{int(x):d}" for x in MarkPeriods])

    ax.set_yscale('log')
    ax.set_ylim(bottom=ybot)

    #  Mark periods
    for T in MarkPeriods:
        ax.axvline(x=1./T, c='k', ls=':', lw=0.5)

    #   Error bar
    x1, y1 = 3e-2, yerr
    err_low, err_hig = y1 *relerr_low, y1 *relerr_high
    ax.plot(x1,y1,'ko')
    ax.plot([x1,x1],[err_low,err_hig],c='k')
    ax.text(x1*1.1,y1,f'{alpha_spec*1e2:3.1f}%')

#}}}

def Compute_Averaged_Spectrum(infile,l0,l1,z0,z1,nperseg,alpha_spec): #{{{
    #------------------------------------------------------------------------------
    #                    Read In
    #------------------------------------------------------------------------------
    print("IN << "+infile)
    var = xr.open_dataset(infile)['u'].isel(lat=0,time=slice(l0,l1)).sel(lev=slice(z0,z1))

    lon  = np.array(var.lon)
    lev  = np.array(var.lev)
    time = np.array(var.time)

    var = np.array(var)

    nt, nz, nx = np.shape(var)

    #------------------------------------------------------------------------------
    #                    Block average in longitude
    #------------------------------------------------------------------------------
    dx = 5                       #  Average every 5 data in lon for var[{time},{lon}]
    nx1 = (len(lon) //dx) *dx

    tmp = np.array(var[:,:,:nx1]).reshape(len(time),len(lev),-1,dx)  #   First axis

    var = np.nanmean(tmp,axis=3)

    lon = lon[:nx1].reshape(-1, dx).mean(axis=1)

    nx = len(lon)

    #------------------------------------------------------------------------------
    #                    Compute Spectrum
    #------------------------------------------------------------------------------
    freq, dummy = signal.welch(var[:,0,0], nperseg=nperseg)

    positive = np.argwhere(freq > 0)[:,0]
    Pxx = np.zeros(len(positive))
    icount = 0.

    for i in range(nx):
        for k in range(nz):
            if np.sum(np.isnan(var[:,k,i])) == 0:
                freq, tmp_Pxx = signal.welch(var[:,k,i], nperseg=nperseg, window='hann')

                Pxx    += tmp_Pxx[positive]
                icount += 1.

    freq = freq[positive]

    Pxx /= icount

    #
    edof = 2. *(np.floor( nt /int(nperseg/2) )-1) *(8/3)
         #  "*(8/3)" is owing to the use of the Hanning window (Thomson and Emery 2014, Table 5.5, p. 479)
    relerr_low, relerr_high = spec_err(edof, alpha_spec)

    #------------------------------------------------------------------------------
    return Pxx,freq, edof, relerr_low, relerr_high
#}}}

def k_m_Spectrum(ForcingPeriod, input_params, ax, lvls):  #{{{
    #--------------------------------------
    #              Read In
    #--------------------------------------
    psd, coords = read_variable_from_netcdf(
            input_params.DirectoryName+f"/tmp_psd_{int(ForcingPeriod):d}dy.nc", "psd",
            return_coords=True
    )
    ft_p,fz_p,fx_p = coords["ft"], coords["fz"], coords["fx"]

    #--------------------------------------
    #            Get Parameters
    #--------------------------------------
    T_target, exp_name, N2 = \
         input_params.T_target, input_params.exp_name, input_params.N2

    ω_lowcut, ω_highcut, ω_tap = get_filter_omega(T_target)
    ω_cent = 0.5 *(ω_lowcut + ω_highcut)

    k_lowcut, k_highcut, k_tap, m_lowcut, m_highcut, m_tap, l_cut = get_filter_km(T_target)

    #--------------------------------------
    #              Draw Figure
    #--------------------------------------
    var = psd
    krange, mrange = [-2.0e-5,2.0e-5], [-4e-2,4e-2]  #  range in k and m
    ft,fz,fx = ft_p,fz_p,fx_p
    ω_target =  ω_cent

    txt_psd_unit = r'[cm$^2$ s$^{-2}$ (day$^{-1}$)$^{-1}$ (m$^{-1}$)$^{-2}$]'

    #    Interpolate PSD onto the target frequency
    interp_fun = interp1d(ft,var, axis=0, kind='linear', fill_value="extrapolate")
    var_intp = interp_fun( ω_target )

    #     contour
    cs = ax.contourf(fx,fz,var_intp,lvls, cmap='jet',extend='both', locator=ticker.LogLocator())

    #     color bar
    cbar = plt.colorbar(cs,ax=ax,location='right',pad=0.10)
    cbar.ax.tick_params(labelsize=9)
    cbar.set_label(txt_psd_unit,fontsize=9)

    #    axes
    ax = define_figure_axes(ax, krange, mrange)

    #    Dispersion relations
    Overplot_DispersionRelation(np.abs(ω_cent), fx, ax)

    #   Mark the target ω and m
    dk,dm = 0.1e-5, 0.1e-2
    i = 0
    k0,k1,m0,m1 = k_lowcut[i],k_highcut[i], m_lowcut[i], m_highcut[i]
    rect = patches.Rectangle(
        (k0,m0),     # bottom-left corner
        (k1-k0),     # width
        (m1-m0),     # height
        linewidth=1,
        edgecolor='red',
        facecolor='none'
    )
    ax.add_patch(rect)
    #ax.text(k0+dk,m0+dm,f'{i+1:d}',color='r',weight='bold')

    return
#}}}

def x_z_Variance(ForcingPeriod, input_params, ax, lvls, n_reflect=0): #{{{
    #--------------------------------------
    #              Read In
    #--------------------------------------
    var_filtered_variance, coords = read_variable_from_netcdf(
        input_params.DirectoryName+f"/tmp_variance_{int(ForcingPeriod):d}dy.nc",
        "var_filtered_variance", return_coords=True
    )
    lev, lon = coords["lev"], coords["lon"]

    #--------------------------------------
    #            Get Parameters
    #--------------------------------------
    N2 = input_params.N2

    #--------------------------------------
    #              Draw Figure
    #--------------------------------------
    lon_Ticks = np.arange(0,61,20)
    lon_TickNames = [f"{int(x)}\N{degree sign}" for x in lon_Ticks]

    iplot = 0

    #ax.set_title(txt_labels[iplot+1]+f"  Box {iplot+1:d}")
    cs=ax.contourf(lon,lev,var_filtered_variance[iplot,:,:],lvls,
                   cmap='jet', extend='both')
    cbar = plt.colorbar(cs,ax=ax)
    cbar.set_label(r'Variance (cm$^2$ s$^{-2}$)')
#    cbar.set_ticks(cbar_ticks)
    ax.tick_params(right=True,top=True,which='both', labelsize=9)

    ax.set_xlabel('Longitude', fontsize=10)
    ax.set_ylabel('Depth (m)', fontsize=10)
    ax.set_xticks(lon_Ticks)
    ax.set_xticklabels(lon_TickNames)
    ax.xaxis.set_minor_locator(AutoMinorLocator(4))

    ax.set_ylim([5000.,0.])
    ax.yaxis.set_ticks(np.arange(0., 5001., 1000.))
    ax.yaxis.set_minor_locator(AutoMinorLocator(5))

    xend, zend = [lon[0],lon[-1]], [0.,lev[-1]]
    x0 = 0.5*(lon[0] + lon[-1])
    xs1,zs1 = ray_K_LongRossby(ForcingPeriod,N2,xend,zend, 'K', x0=x0, nrepeat=n_reflect)
    ax = draw_arrows_mid_head(ax,xs1,zs1, clr='k')

    return
#}}}

def f2T(x): #{{{
    return 1./x
def T2f(x):
    return 1./x
#}}}

def define_figure_axes(ax, krange, mrange): #{{{
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

    return ax
#}}}

#=================================================================
#=================================================================
if __name__ == "__main__":
    #-------------------------------------------------------------
    #                       Input Parameters
    #-------------------------------------------------------------
    alpha_spec = 0.95

    MY_WORK_ROOT = str( pathlib.Path(os.environ["MY_WORK_ROOT"]) )
    indir = MY_WORK_ROOT+"/data/ofes_exp/results/box_bounded/compiled/"

    current_dir = os.getcwd()                        #  Get current directory path
    IN_DIR_TMP = current_dir.replace("code", "data") #  Replace "code" with "data"

    filename_pdf = 'fig.pdf'

    #-----
    grid_name   = "box_bounded"
    period_name = "2000-2029"
    time_width  = 1024           #  # of samples in the time direction
    zstrt, zend = 500., 5000.    #  Depth range
    lstrt, lend = -4096, None    #  Last 4096 days
    N2          = 5e-6           #  Buoyancy frequency squared [s^-2]; for the dispersion relation

    list_ForcingPeriod = [140, 220]

    #-------------------------------------------------------------
    fig, axes, _ = new_page(nrows=3, ncols=2)
    #-------------------------------------------------------------

    iloop, ForcingPeriod = 0, list_ForcingPeriod[0]

    exp_name    = f'{int(ForcingPeriod):d}dy'
    input_params = \
            InputParams( grid_name, exp_name, period_name, time_width,
                         zstrt,zend, lstrt,lend, N2,
                         ForcingPeriod/3.,            #  T_target (2π/ω3 [day])
                         IN_DIR_TMP
                        )

    Averaged_Spectrum(
            indir+"u_eq_box_bounded_"+exp_name+"_2000-2029.nc", exp_name,
            [ForcingPeriod, ForcingPeriod/2., ForcingPeriod/3.],alpha_spec,
            axes[iloop], ybot=1e-1, yerr=3e1)

    k_m_Spectrum(ForcingPeriod, input_params, axes[iloop+2], 10**np.arange(6.0, 10.6, 0.5))
    x_z_Variance(ForcingPeriod, input_params, axes[iloop+4],     np.arange(0.0,  0.11,0.01),
                 n_reflect=5)

    #-------------------------------------------------------------

    iloop, ForcingPeriod = 1, list_ForcingPeriod[1]

    exp_name    = f'{int(ForcingPeriod):d}dy'
    input_params = \
            InputParams( grid_name, exp_name, period_name, time_width,
                         zstrt,zend, lstrt,lend, N2,
                         ForcingPeriod/3.,   #  T_target (2π/ω3 [day])
                         IN_DIR_TMP,
                        )

    Averaged_Spectrum(
            indir+"u_eq_box_bounded_"+exp_name+"_2000-2029.nc", exp_name,
            [ForcingPeriod, ForcingPeriod/2., ForcingPeriod/3.],alpha_spec,
            axes[iloop], ybot=1e1, yerr=3e2)

    k_m_Spectrum(ForcingPeriod, input_params, axes[iloop+2], 10**np.arange(6.5, 11.1, 0.50) )
    x_z_Variance(ForcingPeriod, input_params, axes[iloop+4],     np.arange(0.0,  1.01,0.10),
                 n_reflect=7)

    #-------------------------------------------------------------
    txt_labels = [f"({c})" for c in string.ascii_lowercase]

    for i,(x0,y0) in enumerate(zip((0.10,  0.60,  0.10, 0.60, 0.10,  0.60),
                                   (0.945, 0.945, 0.63, 0.63, 0.315, 0.315))):
        plt.figtext(x0,y0,txt_labels[i], fontsize=12)

    #-------------------------------------------------------------
    plt.tight_layout()
    plt.savefig(filename_pdf)
    #-------------------------------------------------------------


#  ======>>> No longer in use
#def Averaged_Spectrum(infile,exp_name,MarkPeriods,alpha_spec, ax,ytop=10.,yerr=1.): #{{{
#    #===================================================================================
#    #                    Parameters
#    #===================================================================================
#    z0, z1 = 500., 5000.              #   Depth range of the analysis
#
#    l0, l1 = -4096, None              #   Last 4096 days
#    #l0, l1 = 6000, 10096               #   Day 6000 ~ 10095
#    #l0, l1 = -3650, None              #   Read the last "l0" timesteps
#
#    nperseg = 1024                    #  segment length
#
#    #===================================================================================
#    #                     Compute
#    #===================================================================================
#    Pxx,freq, edof, relerr_low, relerr_high = \
#            Compute_Averaged_Spectrum(infile,l0,l1,z0,z1,nperseg,alpha_spec)
#
#    #===================================================================================
#    #                    Plots
#    #===================================================================================
#    #ax.set_title(panel_title)
#    #ax.set_title(f'Average over longitudes and\n'+\
#    #                     f'{int(z0):4d}-{int(z1):4d} m depths', fontsize=10)
#    ax.plot(freq, Pxx *freq, 'k-', label=exp_name)
#
#    ax.set(xlabel='Frequency (cpd)', \
#                   ylabel=r'Variance (cm$^2$ s$^{-2}$)')
#
#    ax.set_xlim([1./500.,1./10.])
#    ax.set_xscale('log')
#    xax2 = ax.secondary_xaxis('top', functions=(f2T,T2f))
#    xax2.set_xlabel('Period (day)')
#    #
#    xax2.xaxis.set_ticks(MarkPeriods)
#    xax2.xaxis.set_ticklabels([f"{int(x):d}" for x in MarkPeriods])
#
#    for T in MarkPeriods:
#        ax.axvline(x=1./T, c='k', ls=':', lw=0.5)
#
#    ax.set_ylim(bottom=0.,top=ytop)
#
#    #   Error bar
#    x1, y1 = 3e-2, yerr
#    err_low, err_hig = y1 *relerr_low, y1 *relerr_high
#    ax.plot(x1,y1,'ko')
#    ax.plot([x1,x1],[err_low,err_hig],c='k')
#    ax.text(x1*1.1,y1,f'{alpha_spec*1e2:3.1f}%')
#
#    #   Axis
#    ax.tick_params(right=True,which='both')
#    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
#}}}
#  <<<====== No longer in use

