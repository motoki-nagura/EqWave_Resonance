"""
    Plot the time series of wavelet power for the three periods
"""
import netCDF4
import numpy as np
import sys,os

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import ticker
from matplotlib.ticker import AutoMinorLocator

from common.new_page import (new_page, )
from themes.MidDepth.OFES_IdealExp.Tools.common.Bandpass_Frequencies import (
        Bandpass_Frequencies, )

exp_name = '180dy+90dy'
#exp_name = '180dy'
#exp_name = '90dy'
#exp_name = 'BottomDamp_180dy'

df_log10 = {"default" : 0.125,   #    Separation in log_10 frequency
            "wide" :    0.175,
            "narrow" :  0.100,} ["default"]

#---
current_dir = os.getcwd()                      #  Get current directory path
IN_DIR = current_dir.replace("code", "data")  #  Replace "code" with "data"

ext = ""
if df_log10 != 0.125:
    ext = f"_df{df_log10:.3f}".replace(".","p")

INPUT_FILE  = IN_DIR+'/out_'+exp_name+ext+'.nc'
OUTPUT_FILE =        './fig_'+exp_name+ext+'.pdf'


def main():
    #==============================================================================
    #             Read In
    #==============================================================================
    print("IN << "+INPUT_FILE)
    f1 = netCDF4.Dataset(INPUT_FILE, 'r')
    lon       = np.copy(f1.variables['lon'])
    lev       = np.copy(f1.variables['lev'])
    time      = np.copy(f1.variables['time'])
    power_avg = np.copy(f1.variables['power_avg'])
    scale_avg = np.copy(f1.variables['scale_avg'])
    period    = np.copy(f1.variables['period'])
    coi       = np.copy(f1.variables['coi'])
    sig95     = np.copy(f1.variables['sig95'])
    T_target  = np.copy(f1.variables['T_target'])
    x0        = np.copy(f1.variables['x0'][0])
    x1        = np.copy(f1.variables['x1'][0])
    z0        = np.copy(f1.variables['z0'][0])
    z1        = np.copy(f1.variables['z1'][0])
    f1.close()

    #        float lon(lon) ;
    #        float lev(lev) ;
    #        double time(time) ;
    #                time:units = "Days since 2000-01-01" ;
    #        float power_avg(power_avg(period, time) ;
    #        float scale_avg(Tband, time, lev, lon) ;
    #        float period(period) ;
    #        float coi(time) ;
    #        float sig95(period, time) ;
    #        float T_target(Tband) ;

    #==============================================================================
    #                      Computation
    #==============================================================================

    scale_avg_xavg  = np.nanmean(scale_avg,      axis=3)
    scale_avg_xzavg = np.nanmean(scale_avg_xavg, axis=2)   #  (Tband, time)

    z0 = 3000.
    k0 = np.argmin( np.abs(lev-z0) )
    scale_avg_zfix = np.squeeze( scale_avg[:,:,k0,:] )

    #-------------------------------------------------------------------------------
    #   Make a version of scale_avg which is NaN if it is in the cone of influence
    #-------------------------------------------------------------------------------
    Tmax = np.zeros(len(T_target))
    #Tmin, Tmax = np.zeros((2,len(T_target)))

    for i,T in enumerate(T_target):
        lc, hc = Bandpass_Frequencies(T)
    #    Tmin[i] = 1./hc
        Tmax[i] = 1./lc

    # 
    scale_avg_xzavg_copy = np.copy(scale_avg_xzavg)

    for p in range(len(T_target)):
        in_coi = np.argwhere(coi < Tmax[p])[:,0]
        scale_avg_xzavg_copy[p,in_coi] = np.nan

    p,l0,l1 = 0,280,290
    #print(f'scale_avg_xzavg_copy[{p:d},{l0:d}:{l1:d}] = ',scale_avg_xzavg_copy[p,l0:l1])
    #print(f'coi[{l0:d}:{l1:d}] = ',coi[l0:l1])
    #print(f'Tmax[{p:d}] = ',Tmax[p])

    #==============================================================================
    #             Plots
    #==============================================================================
    print("OUT << "+OUTPUT_FILE)

    with PdfPages(OUTPUT_FILE) as pdf:

        #------------------------------------------------------------------
        # -----  Scale Averaged Power

        plt.rcParams['font.size'] = 8                    #  Change global font size
        fig, axes, iplot = new_page(1,1, figsize=(7,2))

        #------------------------------------------------------------------

        ax = axes
        #ax = axes[iplot]

        #ax.set_title('(a)')

        TimeSeries_coi = scale_avg_xzavg_copy
        TimeSeries     = scale_avg_xzavg

        # Plot first time series (left y-axis)
        i = 0
        ax.plot(time, TimeSeries_coi[i,:], 'k-', lw=1.0, label=f'{int(T_target[i]):d} days')
        ax.plot(time, TimeSeries[    i,:], 'k:', lw=1.0)
        ax.set_ylabel("Wavelet Power")
        ax.set_xlabel("Time (days)")
        ax.set_xlim([0., time.max()])
        ax.set_ylim(bottom=0.)
        ax.xaxis.set_minor_locator(AutoMinorLocator(4))

        # Create second y-axis sharing the same x-axis
        ax_right = ax.twinx()
        ax_right.set_ylabel("Wavelet Power")

        # Plot second time series (right y-axis)
        i = 1
        ax_right.plot(time, TimeSeries_coi[i,:], 'r-', lw=1.0, label=f'{int(T_target[i]):d} days')
        ax_right.plot(time, TimeSeries[    i,:], 'r:', lw=1.0)

        i = 2
        ax_right.plot(time, TimeSeries_coi[i,:], 'b-', lw=1.0, label=f'{int(T_target[i]):d} days')
        ax_right.plot(time, TimeSeries[    i,:], 'b:', lw=1.0)

        ax_right.set_ylim(bottom=0.)

        #  Add legend
        lines_left,  labels_left  =       ax.get_legend_handles_labels()
        lines_right, labels_right = ax_right.get_legend_handles_labels()

        ax_right.legend(  lines_left  + lines_right,
                         labels_left + labels_right,
                        loc='upper right', bbox_to_anchor=(0.98, 0.15),
                        framealpha=1.0, frameon=False, ncol=3,
                        fontsize='small'
        )

        #---------------------
        plt.tight_layout()
        pdf.savefig()
        plt.close(fig)

        #------------------------------------------------------------------
        # -----  Wavelet Power in time-period space

        plt.rcParams['font.size'] = 8                    #  Change global font size
        fig, axes, iplot = new_page(1,1, figsize=(8,4))

        #------------------------------------------------------------------

        ax = axes
        #ax = axes[iplot]

        dt = 1              #  sampling interval [days]

        lvls = 10**np.arange(-3., 3.6, 0.5)

        #ax.set_title("(b)")

        #           Power
        cs=ax.contourf(time, np.log2(period), power_avg, lvls, extend='both', cmap='jet', \
                               locator=ticker.LogLocator())
        cbar = plt.colorbar(cs,ax=ax, pad=0.03)
        cbar.set_label('Wavelet Power')

        #           95% confidence level
        extent = [time.min(), time.max(), 0, max(period)]
        ax.contour(time, np.log2(period), sig95, [-99, 1],
                           colors='k', linewidths=1, extent=extent)

        #           Cone of influence
        ax.fill(np.concatenate([time, time[-1:]+dt, time[-1:]+dt, time[:1]-dt, time[:1]-dt]),
                np.concatenate([np.log2(coi), [1e-9], np.log2(period[-1:]),
                                np.log2(period[-1:]), [1e-9]]),
                'k', alpha=0.3, hatch='x')

        #           Axis name, range, minor ticks
        ax.set_xlabel('Time (days)')
        ax.set_ylabel('Period (days)')

        ax.set_xlim([time.min(), time.max()])
        ax.set_ylim([np.log2(period.min()), np.log2(period.max())])
        #
        ax.xaxis.set_minor_locator(AutoMinorLocator(4))

        #           Add Y ticks
        Yticks = 2 ** np.arange( np.ceil(np.log2(period.min())),
                                 np.ceil(np.log2(period.max())) )
        ax.set_yticks(np.log2(Yticks))
        ax.set_yticklabels(Yticks)
        #

        #---------------------
        plt.tight_layout()
        pdf.savefig()
        plt.close(fig)

    #    #------------------------------------------------------------------
    #    #------------------------------------------------------------------
    #    #    T - Z  and T-X sections,  scale_avg
    #
    #    fig, axes = new_page(3,1)
    #
    #    trange = [np.min(time), np.max(time)]
    #
    #    for p in range(len(T_target)):
    #        if   p == 0:
    #            lvls = np.arange(0., 1e3+1., 1e2)
    #            #lvls = np.arange(0.,401.,20.)
    #        elif p == 1:
    #            lvls = np.arange(0., 1e2+1., 1e1)
    #            #lvls = np.arange(0., 0.51, 0.05)
    #        elif p == 2:
    #            lvls = np.arange(0., 5e1+1., 5.)
    #            #lvls = np.arange(0., 21., 2.)
    #
    #        axes[iplot].set_title('Zonal average of wavelet power '+\
    #                            f'({int(x0):d}-{int(x1):d}E; '+\
    #                            f'{int(T_target[p]):d} days)')
    #        cs=axes[iplot].contourf(time,lev,np.transpose(scale_avg_xavg[p,:,:]),lvls,\
    #                               cmap='jet',extend='both')
    #        cbar = plt.colorbar(cs,ax=axes[iplot])
    #        axes[iplot].set(xlabel='Time (days)', ylabel='Depth (m)')
    #        axes[iplot].invert_yaxis()
    #        axes[iplot].set_xlim(trange)
    #        axes[iplot].xaxis.set_minor_locator(AutoMinorLocator(5))
    #        iplot += 1
    #
    #    plt.tight_layout()
    #    pdf.savefig()
    #    plt.close(fig)
    #
    #    #------------------------------------------------------------------
    #    #    T - X section,  scale_avg
    #
    #    fig, axes = new_page(3,1)
    #
    #    trange = [np.min(time), np.max(time)]
    #
    #    for p in range(len(T_target)):
    #        if   p == 0:
    #            lvls = np.arange(0., 3e3+1., 3e2)
    #        elif p == 1:
    #            lvls = np.arange(0., 2e2+1., 2e1)
    #        elif p == 2:
    #            lvls = np.arange(0., 2e2+1., 2e1)
    #
    #        axes[iplot].set_title(f'Wavelet power at {int(lev[k0]):d} m for '+\
    #                             f'{int(T_target[p]):d} days')
    #        cs=axes[iplot].contourf(time,lon,np.transpose(scale_avg_zfix[p,:,:]),lvls,\
    #                               cmap='jet',extend='both')
    #        cbar = plt.colorbar(cs,ax=axes[iplot])
    #        axes[iplot].set(xlabel='Time (days)', ylabel='Longitude')
    #        axes[iplot].set_xlim(trange)
    #        axes[iplot].yaxis.set_minor_locator(AutoMinorLocator(2))
    #        iplot += 1
    #
    #    plt.tight_layout()
    #    pdf.savefig()
    #    plt.close(fig)

if __name__ == "__main__":
    sys.exit(main())
