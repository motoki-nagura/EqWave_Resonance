"""
  Plot the spectra of zonal velocity on the equator averaged over
     longitudes and depths and compare the results between EXP_180+90dy
     and EXP_180+90dy_Damp
"""
import numpy as np
import sys,pathlib,os

import matplotlib.pyplot as plt

from common.new_page import new_page

from themes.MidDepth.OFES_IdealExp.Tools.common.Variance_from_PSD import (
    Variance_from_PSD, )
from themes.MidDepth.OFES_IdealExp.Tools.common.Bandpass_Frequencies import (
    Bandpass_Frequencies, )

sys.path.append("..")
from Read_and_Compute_Time_Spectrum import (Read_and_Compute_Time_Spectrum, )


#--------------------

EXP_NAME1 = '180dy+90dy'
EXP_NAME2 = 'BottomDamp_180dy+90dy'

MY_WORK_ROOT = str( pathlib.Path(os.environ["MY_WORK_ROOT"]) )
INDIR = MY_WORK_ROOT+"/data/ofes_exp/results/box_bounded/compiled/"
INPUT_FILENAME_1 = INDIR+"u_eq_box_bounded_"+EXP_NAME1+"_2000-2029.nc"
INPUT_FILENAME_2 = INDIR+"u_eq_box_bounded_"+EXP_NAME2+"_2000-2029.nc"

z0, z1 = 500., 5000.              #   Depth range of the analysis

l0, l1 = -4096, None              #   Last 4096 days

nperseg = 1024                    #  Segment length for Welch method

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
    #        Compute variance
    #------------------------------------
    T_list = [180., 90., 60.]
    vari_1, vari_2 = np.zeros((2,len(T_list)))

    for i, T in enumerate(T_list):
        low_cut, high_cut = Bandpass_Frequencies(T, df_log10=df_log10)
        vari_1[i] = Variance_from_PSD(freq, Pxx_1, low_cut, high_cut)
        vari_2[i] = Variance_from_PSD(freq, Pxx_2, low_cut, high_cut)

    #amp = np.sqrt(vari)

    #------------------------------------
    #        Text output
    #------------------------------------
    txt_out = []

    for i in range(len(T_list)):
        v1, v2 = vari_1[i], vari_2[i]
        ratio = v2 / v1 *1e2
        reduction_rate = 100. - ratio
        txt_out.append(f"Variance at {int(T_list[i]):d} days = "+\
                          f"{v1:6.2f} ("+EXP_NAME1+"), "+\
                          f"{v2:6.2f} ("+EXP_NAME2+"), "+\
                          f" ratio = {ratio:6.2f}%, "+\
                          f" reduction rate = {reduction_rate:6.2f} %")

    #for T in [180., 90., 60.]:
    #    i = np.argmin(np.abs(freq - 1./T))
    #    v1, v2 = Pxx_1[i]*freq[i], Pxx_2[i]*freq[i]
    #    ratio = v2/v1 *1e2
    #    reduction_rate = 100. - ratio
    #    txt_out.append(f"VPS at {int(1./freq[i]):d} days = {v1:6.2f} ("+EXP_NAME1+"), "+\
    #                                                     f"{v2:6.2f} ("+EXP_NAME2+"), "+\
    #                                                     f" ratio = {ratio:6.2f}%, "+\
    #                                                     f" reduction rate = {reduction_rate:6.2f}%")

    with open(OUTFILE_TXT,'w') as fa:
        for text in txt_out:
            fa.write(text)
            fa.write('\n')
            print(text)

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

    ax.plot(freq, var_1, 'k-', label=EXP_NAME1)
    ax.plot(freq, var_2, 'r-', label=EXP_NAME2)

    #   Axis
    ax.set(xlabel='Frequency (cpd)', ylabel=ylabel)

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
    plt.savefig(OUTFILE_FIG)

if __name__ == "__main__":
    sys.exit(main())
