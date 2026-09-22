"""
  Plot the variance of zonal velocity along the equator
"""
from rich import traceback
traceback.install()

import numpy as np
import xarray as xr
from scipy import signal
import sys, pathlib, os

import matplotlib.pyplot as plt
from matplotlib import ticker
from matplotlib.ticker import AutoMinorLocator

from common.new_page import ( new_page, )
from themes.MidDepth.OFES_IdealExp.Tools.common.Bandpass_Frequencies import (
        Bandpass_Frequencies, )
from themes.MidDepth.OFES_IdealExp.Tools.common.EqWaves_Rays import (
        ray_K_LongRossby, ray_Rossby_full, draw_arrows_mid_head, )
from themes.MidDepth.OFES_IdealExp.Tools.common.Variance_from_PSD import (
        Variance_from_PSD, )

MY_WORK_ROOT = str( pathlib.Path(os.environ["MY_WORK_ROOT"]) )

#===================================================================================
#                  User Input
#===================================================================================
GRID_NAME = "box_bounded"

EXP_NAME = "180dy+90dy"

TargetPeriod = 180.
T_target = TargetPeriod / np.array([1.,2.,3.])  #   Target period [day] (e.g., [180.,90.,60.])

DIR_IN = MY_WORK_ROOT+"/data/ofes_exp/results/"+GRID_NAME+"/compiled/"
INPUT_FILENAME = DIR_IN+"u_eq_"+GRID_NAME+"_"+EXP_NAME+"_2000-2029.nc"

#-------
#GRID_NAME = "box_bounded", "box_cyclic"

#EXP_NAME = "180dy+90dy", "180dy", "90dy", "180dy+90dy_N2_1e-5", "180dy+90dy",
#           "BottomDamp_180dy", "BottomDamp_180dy+90dy", "ReadIC_triad12_sol1"
#EXP_NAME, TargetPeriod = "220dy", 220.   #  140, 160, 200, 220

#DIR_IN = "/S/data01/G4006/y0330/data/ofes_exp/results/"+GRID_NAME+"/compiled/"+EXP_NAME[:10]+"/"

#--------
nperseg = 1024

l0, l1 = -4096, None                          #   Last 4096 days
#l0, l1 = 6000, 10096                          #   Day 6000 ~ 10095
#l0, l1 = -3650, None                          #   Read the last "l0" timesteps

df_log10 = {"default" : 0.125,   #    Separation in log_10 frequency
            "wide" :    0.175,   #     Necessary to compute variance
            "narrow" :  0.100,} ["default"]

#--------
if df_log10 == 0.125:
    OUTFILE_FIG = 'fig_Exp'+EXP_NAME+'.pdf'
else:
    OUTFILE_FIG = \
        'fig_Exp'+EXP_NAME+f"_dflog10={df_log10:.3f}".replace('.','p')+'.pdf'

#===================================================================================
#                Constants
#===================================================================================
N2 = 5.0e-6             #  Prescribe the background buoyancy frequency [s^-2]

if EXP_NAME == '180dy+90dy_N2_1e-5':
    N2 = 1e-5            #  N^2 [s^-2], for sensitivity experiments

#===================================================================================
#                Functions
#===================================================================================
def f2T(x):
    return 1./x
def T2f(x):
    return 1./x

#------------------------------------------------------------------------
def main():
    #===================================================================================
    #                    Read In
    #===================================================================================
    var = xr.open_dataset(INPUT_FILENAME)['u'].isel(lat=0,time=slice(l0,l1))

    lon  = np.array(var.lon)
    lev  = np.array(var.lev)
    #time = np.array(var.time)
    var  = np.array(var)

    nt, nz, nx = np.shape(var)
    print('nt = ',nt)

    #===================================================================================
    #                    Block average in longitude
    #===================================================================================
    dx = 5                       #  Average every 5 data in lon for var[{time},{lon}]
    nx1 = (len(lon) //dx) *dx
    tmp1 = np.array(var[:,:,:nx1]).reshape(nt,nz,-1,dx)  #   First axis
    var = np.nanmean(tmp1,axis=3)

    lon = lon[:nx1].reshape(-1, dx).mean(axis=1)

    nx = len(lon)

    #===================================================================================
    #                    Compute Spectrum
    #===================================================================================
    vari = np.zeros((len(T_target),nz,nx)) *np.nan
    #Pxx_180dy, Pxx_90dy, Pxx_60dy = np.zeros((3,nz,nx)) *np.nan

    for i in range(nx):
        for k in range(nz):
            if np.sum(np.isnan(var[:,k,i])) == 0:
                freq, tmp_Pxx = signal.welch(var[:,k,i], nperseg=nperseg, window='hann')
                                                                        #  Compute spectrum

                for p in range(len(T_target)):
                    f_low, f_high = Bandpass_Frequencies(T_target[p])   #  Frequency range

                    if np.logical_and(i==1, k==1):                      #  Print
                        print('1./f_high, T_target[p], 1./f_low = ',
                             f'{int(np.round(1./f_high)):d}, '+\
                             f'{int(np.round(T_target[p])):d}, '+\
                             f'{int(np.round(1./f_low)):d}')

                    vari[p,k,i] = Variance_from_PSD(freq, tmp_Pxx, f_low, f_high)

    #            ifreq = np.argmin(np.abs(freq - 1./T_180dy))
    #            Pxx_180dy[k,i] = tmp_Pxx[ifreq]
    #            T_act_180dy = 1./freq[ifreq]
    #
    #            ifreq = np.argmin(np.abs(freq - 1./T_90dy))
    #            Pxx_90dy[k,i] = tmp_Pxx[ifreq]
    #            T_act_90dy = 1./freq[ifreq]
    #
    #            ifreq = np.argmin(np.abs(freq - 1./T_60dy))
    #            Pxx_60dy[k,i] = tmp_Pxx[ifreq]
    #            T_act_60dy = 1./freq[ifreq]

    #===================================================================================
    #                    Plots
    #===================================================================================
    print("OUT >> "+OUTFILE_FIG)

    lvls_180dy = 10**np.arange( 0.0, 2.1, 0.1)   #   Variance
    lvls_90dy  = 10**np.arange(-1.0, 1.1, 0.1)
    lvls_60dy  = 10**np.arange(-1.5, 0.4, 0.1)

    cbart_180dy = 10**np.arange( 0.0, 2.1, 1.0)
    cbart_90dy  = 10**np.arange(-1.0, 1.1, 1.0)
    cbart_60dy  = 10**np.arange(-1.0, 0.1, 1.0)

    #lvls_180dy = 10**np.arange( 2.0, 5.1, 0.25)  #   Spectral power
    #lvls_90dy  = 10**np.arange( 1.5, 3.6, 0.25)
    #lvls_60dy  = 10**np.arange( 1.0, 3.1, 0.2)
    #
    #cbart_180dy = 10**np.arange(2.0, 5.1, 1.0)
    #cbart_90dy  = 10**np.arange(2.0, 3.1, 1.0)
    #cbart_60dy  = 10**np.arange(1.0, 3.1, 1.0)

    nrows, ncols = 2, 2
    fig, subplots, iplot = new_page(nrows=nrows, ncols=ncols)
    txt_labels = ['(a)','(b)','(c)']

    LonTicks = np.arange(0,61,20)
    LonTickNames = np.char.mod('%d', LonTicks)
    LonTickNames = [s + '\N{degree sign}' for s in LonTickNames]

    #
    for vari_1, lvls, cbar_ticks, Tper in \
            zip(vari, \
                [ lvls_180dy,  lvls_90dy,  lvls_60dy], \
                [cbart_180dy, cbart_90dy, cbart_60dy], \
                T_target):

        ax = subplots[iplot]

        ax.set_title(txt_labels[iplot])
        #ax.set_title(txt_labels[iplot]+"  "+fr'T={int(Tper):d} days')

        cs=ax.contourf(lon,lev, vari_1, lvls,\
                       locator=ticker.LogLocator(),cmap='jet',extend='both')
        ax.set(xlabel='Longitude',ylabel='Depth (m)')
        cbar = plt.colorbar(cs,ax=ax, location='bottom')
        cbar.set_label(r'Variance (cm$^2$ s$^{-2}$)')
        cbar.set_ticks(cbar_ticks)
        #
        ax.set_xlim([0.,62.])
        ax.set_xticks(LonTicks)
        ax.set_xticklabels(LonTickNames)
        ax.xaxis.set_minor_locator(AutoMinorLocator(4))
        #
        ax.set_ylim([5000.,0.])
        #ax.set_ylim([250.,0.])
        ax.yaxis.set_minor_locator(AutoMinorLocator(5))
        #
        ax.tick_params(axis='both', which='both', direction='out',
                       top=True, bottom=True, left=True, right=True)
        #     Ray

        xend, zend = [lon[0],lon[-1]], [lev[0],lev[-1]]
        if   iplot == 0:
            x0 = 0.5*(lon[0] + lon[-1])
            xs1,zs1 = ray_K_LongRossby(Tper,N2,xend,zend, 'K', x0=x0, nrepeat=6)
            ax = draw_arrows_mid_head(ax,xs1,zs1, clr='k')

        elif iplot == 1:
            x0 = 0.75*(lon[0] + lon[-1])
            xs1,zs1 = ray_K_LongRossby(Tper,N2,xend,zend, 'K', x0=x0)
            ax = draw_arrows_mid_head(ax,xs1,zs1, clr='k')

        #elif Tper == 60.:
        #    x0, xend2, zend2 = 40., [20.,40.], [1500.,3500.]
        #    xs2,zs2 = ray_K_LongRossby(Tper,N2,xend2,zend2, 'R', x0=x0, nrepeat=1)
        #    ax = draw_arrows_mid_head(ax,xs2,zs2, clr='k')
        #    #
        #    Lx = 750e3                      #   zonal wavelength [m]
        #    k = 2.*π/(-1.*Lx)
        #    x0, x1, z0 = 40., 35., 2000.    #   points of ray
        #    z1 = ray_Rossby_full(Tper,N2, k, x0=x0,x1=x1,z0=z0)
        #    xs, zs = np.array([x0,x1]), np.array([z0,z1])
        #    ax = draw_arrows_mid_head(ax,xs,zs, clr='r')

        #
        iplot += 1

    fig.delaxes(subplots[3])
    #
    plt.tight_layout()
    plt.savefig(OUTFILE_FIG)


if __name__ == "__main__":
    sys.exit(main())
