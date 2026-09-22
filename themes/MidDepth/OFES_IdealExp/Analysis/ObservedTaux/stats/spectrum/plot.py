"""
         Plot the spectrum of zonal wind stress along the equator
"""
from rich import traceback
traceback.install()

import numpy as np
import os, pathlib, sys

import matplotlib.pyplot as plt
from matplotlib import ticker
import cartopy.mpl.ticker as cticker
from matplotlib.ticker import AutoMinorLocator

from common.new_page import new_page
from themes.MidDepth.OFES_IdealExp.Tools.Taux.sub_read_in_taux import (
        Read_In_taux,)

MY_WORK_ROOT = str( pathlib.Path(os.environ["MY_WORK_ROOT"]) )

#======================================================================
#                          Parameters
#======================================================================

infile = MY_WORK_ROOT+'/data/ofes_exp/results/io0p1/compiled/'+\
        'taux_control_io0p1_2000-2023.nc'
dx, dy = 10, 1    #  Average every dx (or dy) data of invar[{time},{lev},{lon}]

nperseg_spectrum = 1024     #  Length of segments in spectral analysis

#======================================================================
#======================================================================
def f2T(x):
    return 1./x
def T2f(x):
    return 1./x


def Compute_spectrum(taux,lon,lat, nperseg=1024):
    import numpy as np
    from scipy import signal

    print(" ====> Compute_spectrum <====")

    freq, taux_den = signal.welch(taux, nperseg=nperseg, axis=0)

    positive = np.argwhere(freq > 0.)[:,0]
    freq     =     freq[positive]
    taux_den = taux_den[positive,:,:]

    taux_vps = taux_den * freq[:,np.newaxis,np.newaxis]

    return taux_den, taux_vps, freq


def main():
    taux,lon,lat,_ = Read_In_taux(infile, dx, dy, xstrt=40., xend=105.)
    _, vp_spec, f_spec = Compute_spectrum(taux,lon,lat, nperseg=nperseg_spectrum)

    #======================================================================
    #                               Plots
    #======================================================================
    lon_formatter = cticker.LongitudeFormatter()

    jeq = np.argmin( np.abs(lat) )

    #----------------------------------------------------------------------------------
    #         Variance Preserving Spectrum in longitude-frequency domain
    #----------------------------------------------------------------------------------
    plt.rcParams['font.size'] = 10
    fig, axes, iplot = new_page(1,2, figsize=(7,4))

    ax = axes[0]

    lvls = 10**np.arange(-2.5, -1.2, 0.1)
    cs0 = ax.contourf(lon,f_spec,vp_spec[:,jeq,:],levels=lvls,locator=ticker.LogLocator(),\
                      cmap='jet',extend='both')

    cbar0 = plt.colorbar(cs0,ax=ax, location='bottom')

    cbar0.set_label(r'Variance (dyn$^2$ cm$^{-4}$)')

    cbar0.set_ticks([10**-2.5, 1e-2, 10**-1.5])
    cbar0.set_ticklabels([r'10$^{-2.5}$', r'10$^{-2}$', r'10$^{-1.5}$'])

    ax.set(xlabel='Longitude', ylabel='Frequency (cpd)')
    #
    ax.set_xlim([40.,100.])
    ax.set_xticks(np.arange(40,101,20))
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    #
    ax.set_yscale('log')
    ax.set_ylim([1./300.,1./10.])
    #
    yax2 = ax.secondary_yaxis('right', functions=(f2T,T2f))
    yax2.set_ylabel('Period (day)')
    yax2.yaxis.set_ticks([180., 90., 50., 30., 10.])
    yax2.yaxis.set_ticklabels(['180', '90', '50', '30', '10'])

    #ax.axhline(y=1./ 50., c='r', ls=':')
    #ax.axhline(y=1./150., c='r', ls=':')

    fig.delaxes(axes[1])

    plt.tight_layout()
    plt.savefig('fig.pdf')


if __name__ == "__main__":
    sys.exit(main())
