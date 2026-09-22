"""
   Plot the spectra of zonal velocity along the equator averaged over
      longitudes and depths and compare the results between 180dy+90dy and
      180dy.
"""
import numpy as np
import sys,pathlib,os
import matplotlib.pyplot as plt
from common.new_page import new_page

sys.path.append("..")
from Read_and_Compute_Time_Spectrum import (Read_and_Compute_Time_Spectrum, )


#--------------------

EXP_NAME_1 = '180dy+90dy'
EXP_NAME_2 = '180dy'

MY_WORK_ROOT = str( pathlib.Path(os.environ["MY_WORK_ROOT"]) )
INDIR = MY_WORK_ROOT+"/data/ofes_exp/results/box_bounded/compiled/"
INPUT_FILENAME_1 = INDIR+"u_eq_box_bounded_"+EXP_NAME_1+"_2000-2029.nc"
INPUT_FILENAME_2 = INDIR+"u_eq_box_bounded_"+EXP_NAME_2+"_2000-2029.nc"

z0, z1 = 500., 5000.              #   Depth range of the analysis
#z0, z1 = 0., 500.                  #   Near the surface

l0, l1 = -4096, None              #   Last 4096 days

nperseg = 1024                    #  Segment length for Welch method

alpha_spec = 0.95                 #  95% for spectrum error bar

filename_pdf = 'fig.pdf'

if z0 == 0. and z1 == 500.:
    filename_pdf = 'fig_surf.pdf'

#===================================================================================
#===================================================================================
def f2T(x):
    return 1./x
def T2f(x):
    return 1./x


def main():
    #------------------------------------
    #        Read and Compute
    #------------------------------------
    freq, Pxx_1, relerr_low, relerr_high,  = \
            Read_and_Compute_Time_Spectrum(
                    INPUT_FILENAME_1,l0=l0,l1=l1,z0=z0,z1=z1,
                    nperseg=nperseg, alpha_spec=alpha_spec,
                    )

    freq, Pxx_2, relerr_low, relerr_high,  = \
            Read_and_Compute_Time_Spectrum(
                    INPUT_FILENAME_2,l0=l0,l1=l1,z0=z0,z1=z1,
                    nperseg=nperseg, alpha_spec=alpha_spec,
                    )

    #------------------------------------
    #         Draw figures
    #------------------------------------

    fig, axes, _ = new_page(nrows=1, ncols=2, figsize=(8,3))

    #  power spectrum density
    var_1, var_2, ylabel, = (Pxx_1, Pxx_2, r'PSD (cm$^2$ s$^{-2}$ cpd$^{-1}$)',)

    # axis, error bar location
    ylimit, yscale, xerr, yerr, = ([1e0, 10**(4.5)], 'log', 3e-2, 1e2,)

    #--------------------------------
    ax = axes[0]

    ax.plot(freq, var_1, 'k-',  label=EXP_NAME_1)
    ax.plot(freq, var_2, 'r--', label=EXP_NAME_2)

    #   Axis
    ax.set(xlabel='Frequency (cpd)', ylabel=ylabel)
    #ax.legend()

    ax.set_xlim([1./500.,1./10.])
    ax.set_xscale('log')
    xax2 = ax.secondary_xaxis('top', functions=(f2T,T2f))
    xax2.set_xlabel('Period (day)')
    xax2.xaxis.set_ticks([180., 90., 60.])
    xax2.xaxis.set_ticklabels(['180', '90', '60'])

    ax.set_ylim(ylimit)
    ax.set_yscale(yscale)

    #   Mark periods
    for T in [180.,90.,60.]:
        ax.axvline(x=1./T, c='k', ls=':', lw=0.5)

    #   Error bar
    err_low, err_hig = yerr *relerr_low, yerr *relerr_high
    ax.plot(xerr,yerr,'ko',ms=4)
    ax.plot([xerr,xerr],[err_low,err_hig],c='k')
    ax.text(xerr*1.1,yerr,f'{alpha_spec*1e2:3.1f}%')

    ax.tick_params(right=True,which='both')

    fig.delaxes(axes[1])

    plt.tight_layout()
    plt.savefig(filename_pdf)

if __name__ == "__main__":
    sys.exit(main())
