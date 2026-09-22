"""
   Plot bicoherence corresponding to the maximum of the real part of bispectrum
     at each (k3, m3) grid
"""
from rich import traceback
traceback.install()

import numpy as np
import sys, os
from scipy.stats import t
import string

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import AutoMinorLocator
import matplotlib.colors as colors

from common.new_page import (new_page,)
from common.sub_NetCDF_IO import read_variable_from_netcdf as read_in
from themes.MidDepth.OFES_IdealExp.Tools.common.Overplot_DispersionRelation import (
        Overplot_DispersionRelation,)

π = np.pi
day2sec = 86400.   #  1 day in seconds
ω_90dy, ω_60dy = 2.*π/(90.*day2sec), 2.*π/(60.*day2sec)   #  [rad second^-1]


#==========================================================================================
#     Input Parameters
#==========================================================================================

#  Filter

FILTER = { "NoSmth" : "NoSmth",
           "3ptBox" : "3ptBox",}["3ptBox"]

#  Component

COMPONENT = {"xyz": "xyz",     #  u ∂u/dx + v ∂u/∂y + w ∂u/∂z
             "xz":  "xz",      #  u ∂u/dx + w ∂u/∂z
             "x":   "x",}["x"] #  u ∂u/dx

#  Statistical significance level

STAT_Sig_Lev = "99%"

#  DOF for bicoherence
#    By Elgar and Guza (1988, IEEE Transactions of Acoustics, Speech, and Signal
#        Processing, Vol. 36, No. 10)
#    ``For the general case of averaging over m ensembles, and merging over
#      n x n bifrequency squares, ... Each bispectral estimate has 2nm degrees of
#      freedom.'

m = 7   #  1024-day half-overlapping segments for 4096 days long 4096/512-1 = 7
n = 3   #  3 x 3 box filter in (k1,m1) and (k2,m2) space
DOF_bicoherence = 2*n*m  #  = 42

BICOH_SIG = {"99%": np.sqrt(9.2 /DOF_bicoherence),
             "95%": np.sqrt(6.0 /DOF_bicoherence),
             "90%": np.sqrt(4.6 /DOF_bicoherence),
             "80%": np.sqrt(3.2 /DOF_bicoherence),}[STAT_Sig_Lev]

#  DOF for synthetic data test

DOF_synthetic = 500   #  500 samples. See Analysis/Bispectrum/calc_Synth/main.py

#==========================================================================================
#     Derived Parameters
#==========================================================================================

EXP_NAMES  = ["180dy+90dy",   "180dy",]             #  Model run
FREQ_NAMES = ["180-90-60day", "180-180-90day", ]    #  Frequency combination

#  Input directory

current_dir = os.getcwd()                             #  Get current directory path
INDIR_DATA = current_dir.replace("code", "data")      #  Replace "code" with "data"
INDIR_DATA = INDIR_DATA.replace("/significance", "/") #  Remove "significance/"

#  Output

OUTFILE_FIG = "./fig_"+FILTER+".pdf"


#==========================================================================================
def Read_and_Process(INFILE_NETCDF, STAT_Sig_Lev, DOF_synthetic,):
    #==========================================================================================
    #  Read In
    #==========================================================================================
    BSR, coords   = read_in( INFILE_NETCDF, "BSR_k3m3", return_coords=True, )
    BC            = read_in( INFILE_NETCDF, "BC_k3m3", )
    BSR_mean_sync = read_in( INFILE_NETCDF, "BSR_mean_sync_k3m3", )
    BC_sync       = read_in( INFILE_NETCDF, "BC_sync_k3m3", )
    BSR_vari_sync = read_in( INFILE_NETCDF, "BSR_vari_sync_k3m3", )

    m_BSR = coords["m_BSR"]
    k_BSR = coords["k_BSR"]

    #==========================================================================================
    #  Maximum value on each (m3,k3) grid
    #==========================================================================================
    idx = np.nanargmax(BSR, axis=0)

    BSR_max              = np.nanmax(BSR, axis=0)
    BC_at_max            = np.take_along_axis(BC,            idx[None, :, :], axis=0)[0]
    BSR_mean_sync_at_max = np.take_along_axis(BSR_mean_sync, idx[None, :, :], axis=0)[0]
    BC_sync_at_max       = np.take_along_axis(BC_sync      , idx[None, :, :], axis=0)[0]
    BSR_vari_sync_at_max = np.take_along_axis(BSR_vari_sync, idx[None, :, :], axis=0)[0]

    scale = np.nanmax(BSR)
    print("scale = ", scale)

    #==========================================================================================
    #  Confidence interval for real bispectrum
    #==========================================================================================
    stddev = np.sqrt(BSR_vari_sync_at_max)

    alpha = 1. - int(STAT_Sig_Lev[0:2]) *0.01

    if not np.isclose(100.*(1-alpha), float(STAT_Sig_Lev[0:2]), rtol=1e-3):
        raise ValueError("100 x (1-α) /= STAT_Sig_Lev. The former = "
                         f"{100.*(1-alpha):.2f}, "
                         f"the latter = {float(STAT_Sig_Lev[0:2]):.2f}")

    tval = t.ppf(1.-alpha*0.5, DOF_synthetic-1)

    ##  Check t.ppf. Compare with Table A4.3A in Thomson and Emery (2014)
    #print("t.cdf(1.-0.05/2=0.975, 18) = ",t.ppf(1.-0.05/2., 18))  #  = 2.100922...
    #print("t.cdf(1.-0.10/2=0.950, 26) = ",t.ppf(1.-0.10/2., 26))  #  = 1.705617...

    conf_int = tval *stddev /np.sqrt(DOF_synthetic) #  Confidence interval for mean
                                                    #  Eq.(3.40) in Thomson and Emery (2014)

    #==========================================================================================
    return ( m_BSR,k_BSR,
            BSR_max, BC_at_max,
            BSR_mean_sync_at_max, BC_sync_at_max, conf_int,
            scale,
            )


def main():
    EXP_NAME, FREQ_NAME = "180dy+90dy", "180-90-60day"
    DIR_SUB             = "/results/Exp_"+EXP_NAME+"/"+FREQ_NAME+"/"+COMPONENT+"/"
    (m_BSR,k_BSR, BSR_max_1, BC_at_max_1,
      BSR_mean_sync_at_max_1, BC_sync_at_max_1, conf_int_1,
      scale_1,) = \
              Read_and_Process(INDIR_DATA+DIR_SUB+"out_"+FILTER+".nc",
                               STAT_Sig_Lev, DOF_synthetic,)

    EXP_NAME, FREQ_NAME = "180dy", "180-180-90day"
    DIR_SUB             = "/results/Exp_"+EXP_NAME+"/"+FREQ_NAME+"/"+COMPONENT+"/"
    (m_BSR,k_BSR, BSR_max_2, BC_at_max_2,
      BSR_mean_sync_at_max_2, BC_sync_at_max_2, conf_int_2,
      scale_2,) = Read_and_Process(INDIR_DATA+DIR_SUB+"out_"+FILTER+".nc",
                                   STAT_Sig_Lev, DOF_synthetic,)

    #==========================================================================================
    #                    Figure
    #==========================================================================================

    xlab, ylab = "$k_3$ (m$^{-1}$)", "$m_3$ (m$^{-1}$)"
    xax_range,yax_range = [-2e-5,2e-5], [-4e-2,4e-2]
    xcd,ycd = k_BSR,m_BSR

    xcd_ext = np.append(xcd, xcd[-1]+(xcd[-1]-xcd[-2]))

    nrows, ncols = 2, 2

    #----------

    print('OUT >>> '+OUTFILE_FIG)

    with PdfPages(OUTFILE_FIG) as pdf:

        #-------------------------------------------------------------------------------------
        #-------------------------------------------------------------------------------------
        fig, axes, iplot = new_page(nrows, ncols)

        titles = [f"({c})" for c in string.ascii_lowercase]

        #-------------------------------------------------------------------------------------
        #   Bicoherence at BSR = max

        clab, cmap, lvls = 'Bicoherence', 'seismic', np.arange(0.,1.1,0.1)

        for arr, ω, in zip((BC_at_max_1, BC_at_max_2), (ω_60dy, ω_90dy,),):
            ax = axes[iplot]
            cs = ax.contourf(xcd,ycd,arr,lvls, cmap=cmap, extend='both')
            fig.colorbar(cs, ax=ax, label=clab, location='bottom')
            ax.contour(xcd,ycd,arr,levels=[BICOH_SIG],)
            ax = set_ax(ax, titles[iplot], xax_range, yax_range, xlab,ylab)
            Overplot_DispersionRelation(np.abs(ω), xcd_ext, ax)
            iplot += 1

        #-------------------------------------------------------------------------------------
        #   Confidence interval from synthetic data test

        clab, cmap, lvls, cbar_ticks = \
                ('$B$ / '+STAT_Sig_Lev+' Conf. Int.', 'seismic', np.arange(0.,2.1,0.1),
                 np.arange(0,2.1,0.5),)

        for arr, ω, in zip(
                (BSR_max_1/conf_int_1, BSR_max_2/conf_int_2,), (ω_60dy, ω_90dy,),):
            ax = axes[iplot]
            cs = ax.contourf(xcd,ycd,arr,lvls, cmap=cmap, extend='both')
            fig.colorbar(cs, ax=ax, label=clab, location='bottom', ticks=cbar_ticks)
            ax = set_ax(ax, titles[iplot], xax_range, yax_range, xlab,ylab)
            Overplot_DispersionRelation(np.abs(ω), xcd_ext, ax)
            iplot += 1

        #-------------------------------------------------------------------------------------
        plt.tight_layout()
        pdf.savefig()
        plt.close(fig)


        #-------------------------------------------------------------------------------------
        #-------------------------------------------------------------------------------------
        for arr1, arr2, arr3, arr4, ω in zip(
                (              BC_sync_at_max_1,               BC_sync_at_max_2, ),
                (            conf_int_1/scale_1,             conf_int_2/scale_2, ),
                (BSR_mean_sync_at_max_1/scale_1, BSR_mean_sync_at_max_2/scale_2, ), 
                (             BSR_max_1/scale_1,              BSR_max_2/scale_2, ),
                (                        ω_60dy,                         ω_90dy, ),
                ):

            fig, axes, iplot = new_page(nrows, ncols)

            #-------------------------------------------------------------------------------------
            #   Bicoherence from synthetic data test

            ttl, clab, cmap = \
                    'Bicoherence, synthetic', 'bicoherence', 'seismic'

            ax = axes[iplot]
            cs = ax.contourf(xcd,ycd,arr1,np.arange(0.,1.1,0.1), cmap=cmap, extend='both')
            fig.colorbar(cs, ax=ax, label=clab, location='bottom')
            ax.contour(xcd,ycd,arr1,levels=[BICOH_SIG],)
            ax = set_ax(ax, ttl, xax_range, yax_range, xlab,ylab)
            Overplot_DispersionRelation(np.abs(ω), xcd_ext, ax)
            iplot += 1

            #-------------------------------------------------------------------------------------
            #   Confidence interval synthetic data test

            ttl, clab, cmap, v0 = \
                    STAT_Sig_Lev+' Confidence Interval, Scaled', 'BSR', 'RdBu_r', 1.

            norm = colors.SymLogNorm( linthresh=1e0, linscale=1e0, vmin=-v0, vmax=v0)

            ax = axes[iplot]
            pcm = ax.pcolormesh(xcd,ycd,arr2, cmap=cmap, shading='auto', norm=norm)
            fig.colorbar(pcm, ax=ax, label=clab, location='bottom')
            ax = set_ax(ax, ttl, xax_range, yax_range, xlab,ylab)
            Overplot_DispersionRelation(np.abs(ω), xcd_ext, ax)
            iplot += 1

            #-------------------------------------------------------------------------------------
            #   Mean of real bispectrum from synthetic data test

            ttl, clab, cmap, v0 = \
                    'Mean of BSR from test, Scaled', 'BSR', 'RdBu_r', 1.

            norm = colors.SymLogNorm( linthresh=1e0, linscale=1e0, vmin=-v0, vmax=v0)

            ax = axes[iplot]
            pcm = ax.pcolormesh(xcd,ycd,arr3, cmap=cmap, shading='auto', norm=norm)
            fig.colorbar(pcm, ax=ax, label=clab, location='bottom')
            ax = set_ax(ax, ttl, xax_range, yax_range, xlab,ylab)
            Overplot_DispersionRelation(np.abs(ω), xcd_ext, ax)
            iplot += 1

            #-------------------------------------------------------------------------------------
            #   Maximum of real bispectrum at each (k3,m3) grid

            ttl, clab, cmap, v0 = \
                    'Real Bispectrum, Scaled', 'BSR', 'RdBu_r', 1.

            norm = colors.SymLogNorm( linthresh=1e0, linscale=1e0, vmin=-v0, vmax=v0)

            ax = axes[iplot]
            pcm = ax.pcolormesh(xcd,ycd,arr4, cmap=cmap, shading='auto', norm=norm)
            fig.colorbar(pcm, ax=ax, label=clab, location='bottom')
            ax = set_ax(ax, ttl, xax_range, yax_range, xlab,ylab)
            Overplot_DispersionRelation(np.abs(ω), xcd_ext, ax)
            iplot += 1


            plt.tight_layout()
            pdf.savefig()
            plt.close(fig)


def set_ax(ax, ttl, xax_range, yax_range, xlab,ylab):
    ax.set_title(ttl)
    ax.set_xlim(xax_range)
    ax.set_ylim(yax_range)
    ax.axvline(x=0.,lw=0.5,c='k')
    ax.axhline(y=0.,lw=0.5,c='k')
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.xaxis.set_ticks(np.arange(-2.0e-5, 2.1e-5, 1.0e-5))
    ax.yaxis.set_ticks(np.arange(-0.04, 0.041, 0.02))
    ax.yaxis.set_minor_locator(AutoMinorLocator(4))
    ax.set(xlabel=xlab, ylabel=ylab)
    return ax


#  Script entry point
if __name__ == "__main__":
    sys.exit(main())
