"""
   Plot the spectra of zonal velocity along the equator averaged over
      longitudes and depths obtained from Exp180dy+90dy
"""
import numpy as np
import sys, pathlib, os

import matplotlib.pyplot as plt

from common.new_page import new_page

from themes.MidDepth.OFES_IdealExp.Tools.common.Variance_from_PSD import (
    Variance_from_PSD, )
from themes.MidDepth.OFES_IdealExp.Tools.common.Bandpass_Frequencies import (
    Bandpass_Frequencies, )

sys.path.append("..")
from Read_and_Compute_Time_Spectrum import (Read_and_Compute_Time_Spectrum, )


#--------------------

EXP_NAME = '180dy+90dy'

MY_WORK_ROOT = str( pathlib.Path(os.environ["MY_WORK_ROOT"]) )
INPUT_FILENAME = MY_WORK_ROOT+"/data/ofes_exp/results/box_bounded/compiled/"+\
                 "u_eq_box_bounded_"+EXP_NAME+"_2000-2029.nc"

z0, z1 = 500., 5000.              #   Depth range of the analysis

l0, l1 = -4096, None              #   Last 4096 days
#l0, l1 = 6000, 10096              #   Day 6000 ~ 10095
#l0, l1 = -3650, None              #   Read the last "l0" timesteps

nperseg = 1024                    #  Segment length for Welch method
#nperseg = 365

alpha_spec = 0.95                 #  95% for spectrum error bar

df_log10 = {"default" : 0.125,    #  Separation in log_10 frequency
            "wide" :    0.175,    #    necessary to compute variance
            "narrow" :  0.100,} ["default"]

OUTFILE_FIG = 'fig.pdf'

if df_log10 == 0.125:
    OUTFILE_TXT = 'out.txt'
else:
    OUTFILE_TXT = f"out_dflog10={df_log10:.3f}".replace('.','p')+'.txt'


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
    freq, Pxx, relerr_low, relerr_high,  = \
            Read_and_Compute_Time_Spectrum(
                    INPUT_FILENAME,l0=l0,l1=l1,z0=z0,z1=z1,
                    nperseg=nperseg, alpha_spec=alpha_spec,
                    )

    #------------------------------------
    #        Compute variance
    #------------------------------------
    T_list = [180., 90.]
    vari = np.zeros(len(T_list))

    for i, T in enumerate(T_list):
        low_cut, high_cut = Bandpass_Frequencies(T, df_log10=df_log10)
        vari[i] = Variance_from_PSD(freq, Pxx, low_cut, high_cut)

    #------------------------------------
    #        Text output
    #------------------------------------
    tmp, text_out = [], []

    #  VPS at fixed frequencies

    for i, T in enumerate(T_list):
        text_out.append(f"T = {T:.2f},  variance = {vari[i]:.2f} cm s^-1")
        tmp.append(vari[i])

    text_out.append(f'ratio = {tmp[1]/tmp[0]*1e2:.1f} %')

    #  Periods in resonance with 180-day variability

    T = np.zeros(20)
    T[0] = 180.
    for i in range(1,len(T)):
        if i == 1:
            T[i] = 1./(1./T[0] + 1./T[0])
        else:
            T[i] = 1./(1./T[0] + 1./T[i-1])
    text_out.append("Periods in resonance with 180-day variability = "+
          ", ".join(f"{v:.2f}" for v in T[:])+" days"
    )

    T_res_180dy = T

    with open(OUTFILE_TXT,'w') as fa:
        for text in text_out:
            fa.write(text)
            fa.write('\n')
            print(text)

    #------------------------------------
    #         Draw figures
    #------------------------------------
    fig, axes, _ = new_page(nrows=1, ncols=2, figsize=(8,3))

    #>>>>>>>>>>>>
    #  power spectrum density
    var, ylabel, = (Pxx, r'PSD (cm$^2$ s$^{-2}$ cpd$^{-1}$)',)

    # axis, error bar location
    ylimit, yscale, xerr, yerr, = ([1e0, 10**(4.5)], 'log', 3e-2, 1e2,)

    #--------------------------------
    ax = axes[0]

    ax.plot(freq, var, 'k-')

    #  Axis settings
    ax.set(xlabel='Frequency (cpd)', ylabel=ylabel)
    ax.tick_params(right=True,which='both')

    ax.set_xlim([1./500.,1./10.])
    ax.set_xscale('log')
    xax2 = ax.secondary_xaxis('top', functions=(f2T,T2f))
    xax2.set_xlabel('Period (day)')
    xax2.xaxis.set_ticks([180., 90., 60.])
    xax2.xaxis.set_ticklabels(['180', '90', '60'])

    ax.set_ylim(ylimit)
    ax.set_yscale(yscale)

    #   Mark periods
    for T in T_res_180dy[0:3]:
#    for T in T_res_180dy:
        ax.axvline(x=1./T, c='k', ls=':', lw=0.5)

    #   Error bar
    err_low, err_hig = yerr *relerr_low, yerr *relerr_high
    ax.plot(xerr,yerr,'ko',ms=4)
    ax.plot([xerr,xerr],[err_low,err_hig],c='k')
    ax.text(xerr*1.1,yerr,f'{alpha_spec*1e2:3.1f}%')

    fig.delaxes(axes[1])
    #<<<<<<<<<<<<<

    plt.tight_layout()
    plt.savefig(OUTFILE_FIG)

if __name__ == "__main__":
    sys.exit(main())
