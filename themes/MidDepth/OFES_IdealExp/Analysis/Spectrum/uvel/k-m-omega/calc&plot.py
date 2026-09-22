"""
   Plot the spectral power of zonal velocity along the equator in k-m space
     at three periods and superpose dispersion curves
"""

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                     imports
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
from rich import traceback
traceback.install()

import numpy as np
import xarray as xr
from scipy.signal import detrend
from scipy.interpolate import interp1d
import sys, os, pathlib
from dataclasses import dataclass, field

import matplotlib.pyplot as plt
from matplotlib import ticker
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import AutoMinorLocator
import matplotlib.patches as patches
import time as time_measure

from common.new_page import new_page
from common.fft_spectrum_3d import fft_spectrum_3d
from themes.MidDepth.OFES_IdealExp.Tools.common.Overplot_DispersionRelation import (
        Overplot_DispersionRelation,
        )
from themes.MidDepth.OFES_IdealExp.Tools.common.Fixed_Parameters import (
        degrees_to_meters,
        )
from themes.MidDepth.OFES_IdealExp.Analysis.Filtering.filtering.Below500m.main_filtering import (
        get_filter_km_60days, )

MY_WORK_ROOT = str( pathlib.Path(os.environ["MY_WORK_ROOT"]) )

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                   constants
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
π       = np.pi
deg2met = degrees_to_meters() #  1 deg in meters

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                   functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def f2T(x):
    return 1./x
def T2f(x):
    return 1./x

@dataclass(slots=True)
class PassBands:
    k_lowcut:  np.ndarray = field(default_factory=lambda: np.array([]))
    k_highcut: np.ndarray = field(default_factory=lambda: np.array([]))
    m_lowcut:  np.ndarray = field(default_factory=lambda: np.array([]))
    m_highcut: np.ndarray = field(default_factory=lambda: np.array([]))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                   main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def main():
    #============================================================================================#
    #           Input Parameters
    elapsed_time_start = time_measure.perf_counter()

    GRID_NAME = "box_bounded"
    #GRID_NAME = "box_cyclic"

    #     Set experiment name
    #EXP_NAME = "180dy+90dy", "BottomDamp_180dy+90dy", "180dy", "BottomDamp_180dy",
    #           "90dy", "BottomDamp_180dy", "ReadIC_triad1_sol1", "ReadIC_Case3_triad12_sol1",
    #           "180dy+90dy_N2_1e-5", ...
    EXP_NAME = "180dy+90dy"
    Target_Periods = [180., 90., 60.]

    ##    For single frequeny experiment (140, 160, 200, 220)
    #ForcingPeriod = 220.
    #EXP_NAME = f"{int(ForcingPeriod):d}dy"
    #Target_Periods = [ForcingPeriod, ForcingPeriod/2., ForcingPeriod/3.]

    #    Analysis domain
    PERIOD_NAME, time_width = "2000-2029", 1024
    #            time_width = # of samples in the time direction

    z0, z1 = 500., 5000.        #  Depth range of the analysis [m]

    l0,l1  = -4096, None        #  Last 4096 days

    #    Buoyancy frequency
    N2 = 5e-6                   #  Buoyancy frequency [s^-2]
    if EXP_NAME == "180dy+90dy_N2_1e-5": #     sensitivity experiment
        N2 = 1e-5

    #    Output file name
    OUTPUT_FILENAME = 'fig_Exp'+EXP_NAME+'.pdf'

    #============================================================================================#
    #                                Read In Model Zonal Velocity
    #============================================================================================#
    print('Reading...')
    infile = MY_WORK_ROOT+"/data/ofes_exp/results/"+GRID_NAME+"/compiled/"+\
             "u_eq_"+GRID_NAME+"_"+EXP_NAME+"_"+PERIOD_NAME+".nc"
    print("IN << "+infile)

    var = xr.open_dataset(infile)['u'].isel(lat=0,time=slice(l0,l1)).sel(lev=slice(z0,z1))

    lon  = np.array(var.lon)
    lev  = np.array(var.lev)
    var  = np.array(var)

    time = np.arange(0., len(var[:,0,0]))   #  time in days

    #   Remove NaN in the x direction
    i = np.argwhere(~np.isnan(var[0,0,:]))[:,0]

    var = var[:,:,i]
    lon = lon[i]

    #============================================================================================#
    #                                  Compute Spectrum
    #============================================================================================#
    print('Computing spectrum...')

    nt,nz,nx = np.shape(var)     #   # of time steps in the data
    print('nt = ',nt)

    #-------
    lon_met = lon *deg2met                 #  X coordinate, in meters

    #-------
    lev_flipped = np.flip(lev *(-1.))      #  Z coordinate is positive upward
                                           #  Flipped in the vertical direction such that
                                           #     array is from -5000 m to 0 m (say)

    #-------
    hanning_1d = np.hanning(time_width)
    hanning_window = np.tile(hanning_1d[:, None, None], (1, nz,nx))

    #-------
    l0, l1, icnt = 0, time_width, 0.   #  initial l0, l1, and a counter

    psd = np.zeros((time_width,nz,nx))

    while l1 < nt:

        x_in   = lon_met
        z_in   = lev_flipped

        var_in = var[l0:l1,:,:]

        var_in = detrend(var_in, type='constant', axis=0)  #  subtract the temporal average
        var_in = var_in * hanning_window
        var_in = np.flip(var_in, axis=1)    #  Flipped in the vertical direction such that
                                            #     array is from -5000 m to 0 m (say)

        tmp_psd, fx,fz,ft = fft_spectrum_3d(var_in, x_in, z_in, time[l0:l1])
        tmp_psd = tmp_psd *(8./3.)    #  Correction for the Hanning window in time

        psd = psd + tmp_psd

        icnt = icnt + 1.

        l0 += int(time_width /2)   #      advance by time_width/2
        l1 += int(time_width /2)

    psd = psd /icnt
    edof = icnt *2.

    #-------
    fx = 2.*π *fx          #  Zonal wavenumber [cycle per meter] to [rad meter^-1]
    fz = 2.*π *fz          #  Vertical wavenumber [cycle per meter] to [rad meter^-1]

    #============================================================================================#
    #                                         Plots
    #============================================================================================#
    print('Plotting...')

    #--------------------------------------------------------------------

    xrange, yrange = [-2e-5,2e-5], [-4e-2,4e-2]  #  range in k and m

    nrows, ncols = 2, 2

    PassBands.k_lowcut, PassBands.k_highcut, _, PassBands.m_lowcut, PassBands.m_highcut, \
            _, _ = get_filter_km_60days()

    with PdfPages(OUTPUT_FILENAME) as pdf:
        #---------------------------------------------------------------------
        #            Spectrum

        fig, axes, iplot = new_page(nrows,ncols)
        iplot = 0

        for T in Target_Periods:     #  Loop in period

            ax = axes[iplot]

            #    Contour levels

            if   iplot == 0:
                lvls = 10**np.arange(7.0, 13.1, 0.5)   #  Lowest frequency
            elif iplot == 1:
                lvls = 10**np.arange(6.5, 11.1, 0.5)   #  Middle frequency
            elif iplot == 2:
                lvls = 10**np.arange(6.0, 10.6, 0.5)   #  Highest frequency
                #lvls = 10**np.arange(9.5, 10.7, 0.1)   #  (emphasize high power)

            #    Interpolate PSD onto the target frequency
            interp_fun = interp1d(ft,psd, axis=0, kind='linear', fill_value="extrapolate")
            psd_1 = interp_fun(1./T)
            title = f'{int(T):d} days'
            f_cpd = 1./T

            #ift = np.argmin(np.abs(ft - 1./T))      #  Or, pick the selected period
            #psd_1 = psd[ift,:,:]
            #title = f'{1./ft[ift]:3.1f} days'
            #f_cpd = ft[ift]

            #     title
            ax.text( 0.04, 0.96, title,
                        transform=ax.transAxes, ha="left", va="top",
                        bbox=dict(facecolor="white", edgecolor="black")
            )

            #     contour
            cs = ax.contourf(fx,fz,psd_1, levels=lvls, locator=ticker.LogLocator(),
                             cmap='jet',extend='both')

            #     color bar
            cbar = plt.colorbar(cs,ax=ax,location='bottom',pad=0.17)
            cbar.ax.tick_params(labelsize=9)
            cbar.set_label(r'[cm$^2$ s$^{-2}$ (day$^{-1}$)$^{-1}$ (m$^{-1}$)$^{-2}$]',\
                           fontsize=9)

            #     X and Y axes
            ax.axhline(y=0., c='k',ls='-',lw=0.5)
            ax.axvline(x=0., c='k',ls='-',lw=0.5)

            ax.set_xlim(xrange)
            ax.set_ylim(yrange)

            ax.set_xlabel(r'$k$ (m$^{-1}$)',fontsize=10)
            ax.set_ylabel(r'$m$ (m$^{-1}$)',fontsize=10)

            ax.tick_params(axis='both', which='major', labelsize=9)
            ax.yaxis.set_ticks(np.arange(-0.04, 0.041, 0.02))
            ax.xaxis.set_minor_locator(AutoMinorLocator(5))
            ax.yaxis.set_minor_locator(AutoMinorLocator(4))

            #     Secondary axes
            xax2 = ax.secondary_xaxis('top', functions=(f2T,T2f))    #  Secondary axis on X
            xax2.set_xlabel('Zonal Wavelength (km)',fontsize=10)
            xax2.tick_params(axis='both', which='major', labelsize=9)

            yax2 = ax.secondary_yaxis('right', functions=(f2T,T2f))  #  Secondary axis on Y
            yax2.set_ylabel('Vertical Wavelength (km)',fontsize=10)
            yax2.tick_params(axis='both', which='major', labelsize=9)

            tmp = [-0.4, -0.75, 0.75, 0.4]
            tick_labels = ['-400', '-750', '750', '400']
            tmp1 = np.array(tmp)*1e6 /(2.*π)
            xax2.xaxis.set_ticks(tmp1)
            xax2.xaxis.set_ticklabels(tick_labels)

            tmp = [-250, -1000, 1000, 250]
            tick_labels = ['-0.25', '-1', '1', '0.25']
            tmp1 = np.array(tmp) /(2.*π)
            yax2.yaxis.set_ticks(tmp1)
            yax2.yaxis.set_ticklabels(tick_labels)

            #     Additional marks
            Overplot_DispersionRelation(2.*π *f_cpd /86400., fx, ax, N2=N2)

            if EXP_NAME == "180dy+90dy":
                if   iplot == 0:
                    k_HE_180, m_KE_180 = MarkHighestEnergy_180day(psd_1,fx,fz,iplot,axes,EXP_NAME)
                elif iplot == 2:
                    ax = Mark_PassBands(ax,PassBands)    #  Mark passbands for 60-day variability

                MarkHarmonics(iplot,axes,k_HE_180, m_KE_180)

            iplot += 1

        fig.delaxes(axes[3])

        plt.figtext(0.09, 0.96, '(a)', fontsize=13)
        plt.figtext(0.58, 0.96, '(b)', fontsize=13)
        plt.figtext(0.09, 0.47, '(c)', fontsize=13)

        plt.tight_layout()
        pdf.savefig()
        plt.close(fig)


    #======================================================================
    #======================================================================

    elapsed_time_end = time_measure.perf_counter()
    elapsed = elapsed_time_end - elapsed_time_start
    print(f"Elapsed time: {elapsed:.6f} seconds")

def MarkHighestEnergy_180day(psd_1,fx,fz,iplot,axes,EXP_NAME):
    #   Mark the highest energy at 180-day period
    ifz, ifx = np.unravel_index(np.argmax(psd_1), psd_1.shape) #  find the indices for max
    kmark, mmark = fx[ifx], fz[ifz]
    axes[iplot].plot(kmark, mmark, 'wx')

    #  Write out parameters to a text file
    with open('out_'+EXP_NAME+'.txt','w') as fa:
        fa.write(f"Highest energy at 180-day period is at k={kmark:.3e} and m={mmark:.3e}")
        fa.write('\n')

    return kmark, mmark

def MarkHarmonics(iplot,axes,kmark,mmark,):
    if   iplot == 1:
        axes[iplot].plot(kmark*2, mmark*2, 'wx')
    elif iplot == 2:
        axes[iplot].plot(kmark*3, mmark*3, 'wx')

def Mark_PassBands(ax,PassBands):
    #   Mark the target ω and m
    dk,dm = 0.1e-5, 0.1e-2
    for i in range(len(PassBands.k_lowcut)):
        k0,k1 = PassBands.k_lowcut[i], PassBands.k_highcut[i]
        m0,m1 = PassBands.m_lowcut[i], PassBands.m_highcut[i]
        rect = patches.Rectangle(
            (k0,m0),     # bottom-left corner
            (k1-k0),     # width
            (m1-m0),     # height
            linewidth=1,
            edgecolor='red',
            facecolor='none'
        )
        ax.add_patch(rect)
        ax.text(k0+dk,m0+dm,f'{i+1:d}',color='r',weight='bold')

    return ax

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                Script entry point
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if __name__ == "__main__":
    sys.exit(main())
