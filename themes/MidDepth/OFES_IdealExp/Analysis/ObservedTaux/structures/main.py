"""
   * Draw regression coefficients of zonal wind stresses obtained from atmospheric reanalysis
   * Draw idealized wind forcing
"""
from rich import traceback
traceback.install()

import numpy as np
from scipy import signal
from scipy.interpolate import interp1d
import os, pathlib, sys

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import cartopy.mpl.ticker as cticker
from matplotlib.ticker import AutoMinorLocator
import matplotlib.pyplot as plt

from common.regr_wgt import (regr_wgt, )
from func_draw_fig import (
        lon2x_180dy, x2lon_180dy, lon2x_90dy, x2lon_90dy,
        generate_axes, make_smooth_seismic,
        )

from themes.MidDepth.OFES_IdealExp.Tools.Taux.sub_read_in_taux import (
        Read_In_taux,)
from themes.MidDepth.OFES_IdealExp.Tools.common.Bandpass_Frequencies import (
        Bandpass_Frequencies,)

MY_WORK_ROOT = str( pathlib.Path(os.environ["MY_WORK_ROOT"]) )

infile = MY_WORK_ROOT+'/data/ofes_exp/results/io0p1/compiled/'+\
        'taux_control_io0p1_2000-2023.nc'
dx, dy = 10, 1            #  Average every dx (or dy) data of invar[{time},{lev},{lon}]

T_target = [180., 90.]    #  Target periods [days]

df_log10 = {"default" : 0.125,   #    Separation in log_10 frequency
            "wide" :    0.175,
            "narrow" :  0.100,} ["default"]

order_filter  = 5         #  Filter order
sampling_freq = 1.        #  Sampling frequency (= 1 day^-1)


if df_log10 == 0.125:
    OUTFILE_FIG = "fig.pdf"
    OUTFILE_TXT = "out.txt"
else:
    OUTFILE_FIG = f"fig_df{df_log10:.3f}".replace(".","p")+".pdf"


#===========================================================================
def Time_Filtering(T_target, order_filter, sampling_freq, taux, lon, lat):
    print(" >=== Time filtering ===<")

    for p,Tper in enumerate(T_target):
        #    Compute low- and high-cut frequencies
        lc, hc = Bandpass_Frequencies(Tper, df_log10=df_log10, )
        print('p, Tper, 1/hc, 1/lc = ',p, Tper, 1/hc, 1/lc)

        #    Set up the filter
        sos = signal.butter(
                order_filter, [lc, hc], fs=sampling_freq, btype='band', output='sos')

        #    Apply the filter
        results = np.copy(taux) *np.nan

        for j in range(len(lat)):
            for i in range(len(lon)):
                if np.sum(~np.isnan(taux[:,j,i])) != 0:
                    results[:,j,i] = signal.sosfiltfilt(sos, taux[:,j,i])

        #    Save
        if   p == 0:
            taux_180dy = results
        elif p == 1:
            taux_90dy = results

    return taux_180dy, taux_90dy


def Regression_Analysis(taux_180dy, taux_90dy, lon, lat,):
    print(" >=== Regression Analysis ===<")

    regr = np.zeros((len(T_target),len(lat),len(lon))) *np.nan

    for p in range(len(T_target)):
        print('p = ',p)

        if   p == 0:
            var = taux_180dy
            x_index, y_index = 70., 0. #  Location of the index time series in coherence analysis
        elif p == 1:
            var = taux_90dy
            x_index, y_index = 85., 0. #  Location of the index time series in coherence analysis

        #   Index time series
        i_index = np.argmin( np.abs(lon - x_index) )
        j_index = np.argmin( np.abs(lat - y_index) )
        index = var[:,j_index,i_index]
        index = index /np.std(index)     #  normalize by standard deviation

        sig = np.zeros(len(index)) #  dummy
        alpha = 0.95

        for j in range(len(lat)):
            for i in range(len(lon)):
                tmpy = var[:,j,i]

                if np.sum( ~np.isnan(tmpy) ) != 0:
                    edof = len(index)   #   dummy
                    ints, slop, slop_err = \
                            regr_wgt(index, tmpy, sig,False, edof, alpha) #  Regression
                    regr[p,j,i] = slop

    return regr


def Generate_Idealized_Forcing(lon, lat,):
    print(" >=== Generate Idealized Forcing ===<")
    #==============================================================================
    #     Define structructure functions
    #==============================================================================
    Amp_180dy, lon_cent_180dy, lon_width_180dy, lat_width_180dy = 0.15, 70., 50., 15.
    Amp_90dy,  lon_cent_90dy,  lon_width_90dy,  lat_width_90dy  = 0.06, 85., 40., 15.

    #-----------
    x = 2.*np.pi *(lon-lon_cent_180dy) /lon_width_180dy
    xfunc_180dy = Amp_180dy *( 0.5*(np.cos(x) + 1.) )
    xfunc_180dy[ x < -np.pi ] = 0.
    xfunc_180dy[ x >  np.pi ] = 0.

    yfunc_180dy = Amp_180dy * np.exp( -0.5 *(lat**2)/lat_width_180dy )

    #-----
    x = 2.*np.pi *(lon-lon_cent_90dy) /lon_width_90dy
    xfunc_90dy = Amp_90dy *( 0.5*(np.cos(x) + 1.) )
    xfunc_90dy[ x < -np.pi ] = 0.
    xfunc_90dy[ x >  np.pi ] = 0.

    yfunc_90dy = Amp_90dy * np.exp( -0.5 *(lat**2)/lat_width_90dy )

    #==============================================================================
    #     Define structure functions on the grids of box_bounded
    #==============================================================================

    yu = np.arange(402) *0.1  -20.0   #  box_bounded
    xu = np.arange(622) *0.1  + 0.0

    Amp_180dy, lon_cent_180dy, lon_width_180dy, lat_width_180dy = 0.15, 30., 50., 15.
    Amp_90dy,  lon_cent_90dy,  lon_width_90dy,  lat_width_90dy  = 0.06, 45., 40., 15.

    #-----

    x = 2.*np.pi *(xu-lon_cent_180dy) /lon_width_180dy
    xfunc_model_180dy = 0.5*(np.cos(x) + 1.)
    xfunc_model_180dy[ x < -np.pi ] = 0.
    xfunc_model_180dy[ x >  np.pi ] = 0.

    yfunc_model_180dy = np.exp( -0.5 *(yu**2)/lat_width_180dy )

    func_model_180dy = Amp_180dy * yfunc_model_180dy[:,None] * xfunc_model_180dy[None,:]

    #-----

    x = 2.*np.pi *(xu-lon_cent_90dy) /lon_width_90dy
    xfunc_model_90dy = 0.5*(np.cos(x) + 1.)
    xfunc_model_90dy[ x < -np.pi ] = 0.
    xfunc_model_90dy[ x >  np.pi ] = 0.

    yfunc_model_90dy = np.exp( -0.5 *(yu**2)/lat_width_90dy )

    func_model_90dy = Amp_90dy * yfunc_model_90dy[:,None] * xfunc_model_90dy[None,:]

    #==============================================================================
    #     Compute the ratio and write out to a text file
    #==============================================================================
    ratio_vari = Amp_90dy**2 / Amp_180dy**2  #  Ratio of variance. The factor 1/2 is
                                             #    omitted because it does not affect the
                                             #    ratio
    text_out = f"Ratio of variance of taux at 90 day to that at 180 day is "+\
               f"{ratio_vari *1e2:.1f} %"

    if df_log10 == 0.125:
        with open(OUTFILE_TXT,'w') as fa:
            fa.write(text_out)
            fa.write('\n')
            print(text_out)


    return xfunc_180dy, yfunc_180dy, xfunc_90dy, yfunc_90dy, \
            func_model_180dy, func_model_90dy


def Draw_Figure(
        OUTFILE_FIG,
        regr,lon,lat,
        xfunc_180dy, yfunc_180dy, xfunc_90dy, yfunc_90dy,
        func_model_180dy, func_model_90dy,
        period, xu,yu,
        ):
    print(" >=== Draw Figure ===<")
    jeq  = np.argmin(np.abs(lat))
    i70e = np.argmin(np.abs(lon-70.))
    i85e = np.argmin(np.abs(lon-85.))

    for i in range(i70e,0,-1):
        if np.isnan(regr[0,jeq,i]):
            iwb = i
            break

    for i in range(i70e,len(lon),1):
        if np.isnan(regr[0,jeq,i]):
            ieb = i
            break

    print('Longitudes for the eastern and western boundaries are ',
          lon[ieb],' and ',lon[iwb],' E')

    #==============================================================================

    plt.rcParams['font.size'] = 9

    labels = ["(a)","(c)","(e)",
              "(b)","(d)","(f)"]

    lon_formatter = cticker.LongitudeFormatter()
    lat_formatter = cticker.LatitudeFormatter()

    XTicks     = np.arange(0,61,20)
    XTickNames = np.char.mod('%d', XTicks)
    XTickNames = [s + '\N{degree sign}' for s in XTickNames]

    print("OUT >> "+OUTFILE_FIG)

    with PdfPages(OUTFILE_FIG) as pdf:
        fig, axes = generate_axes()
        iplot = 0

        #----------

        for p in range(2):

            if p == 0:
                lvls = np.arange(-0.2,0.21,0.02)
                cbar_ticks = np.arange(-0.2,0.21,0.1)
                xindx, yindx = 70., 0.
                ilon0 = i70e
                xfunc, yfunc = xfunc_180dy, yfunc_180dy
                func_model = func_model_180dy
            elif p == 1:
                lvls = np.arange(-0.1,0.11,0.01)
                cbar_ticks = np.arange(-0.1,0.11,0.05)
                xindx, yindx = 85., 0.
                ilon0 = i85e
                xfunc, yfunc = xfunc_90dy, yfunc_90dy
                func_model = func_model_90dy

            cmap, norm = make_smooth_seismic(func_model, width=0.25)

            #      Observed, regressed taux on lon-lat map

            ax = axes[iplot]
            ax.set_title(labels[iplot]+"  "+f'{int(period[p]):d} days')
            ax.set_facecolor('0.7')              #  Background = gray
            cs = ax.contourf(lon,lat,regr[p,:,:],lvls, cmap=cmap,extend='both')
            cbar = plt.colorbar(cs, ax=ax, location='right', ticks=cbar_ticks)
            cbar.set_label(r'dyn cm$^{-2}$')
            ax.plot(xindx,yindx,'k.',ms=12)
            iplot += 1

            #      Idealized forcing in model grids

            ax = axes[iplot]
            ax.set_title(labels[iplot])
            cs = ax.contourf(xu,yu,func_model,lvls, cmap=cmap,norm=norm)
            cbar = plt.colorbar(cs,ax=ax, location='right', ticks=cbar_ticks)
            cbar.set_label(r'dyn cm$^{-2}$')
            iplot += 1

            #      Comparison by line plots

            f = interp1d(lat, regr[p,:,:], axis=0)
            regr_intp = f(0.)        #   interpolate onto the equator

            ax = axes[iplot]
            ax.set_title(labels[iplot])
            ax.plot(lon,regr_intp,'k-', label=r'Regressed $\tau^x$')
            ax.plot(lon,xfunc,    'k:', label='Idealized Structure')
            ax.set(ylabel='dyn cm$^{-2}$')
            ax.axhline(y=0,lw=1.0,c='k')
            iplot += 1

        #  Set axes
        for i in [0,3]:
            axes[i].set_xlim([40.,102.])
            axes[i].set_xticks(np.arange(40,101,20))
            axes[i].set_yticks(np.arange(-20,21,10))
            axes[i].xaxis.set_major_formatter(lon_formatter)
            axes[i].yaxis.set_major_formatter(lat_formatter)
            axes[i].xaxis.set_minor_locator(AutoMinorLocator(2))
            axes[i].yaxis.set_minor_locator(AutoMinorLocator(2))
            axes[i].tick_params(right=True,top=True,which='both')
        for i in [1,4]:
            axes[i].set_xticks(XTicks)
            axes[i].set_xticklabels(XTickNames)
            axes[i].yaxis.set_major_formatter(lat_formatter)
            axes[i].yaxis.set_minor_locator(AutoMinorLocator(2))
            axes[i].tick_params(right=True,top=True,which='both')
        for i in [2,5]:
            axes[i].set_xlim([40.,102.])
            axes[i].set_xticks(np.arange(40,101,20))
            axes[i].xaxis.set_major_formatter(lon_formatter)
            axes[i].xaxis.set_minor_locator(AutoMinorLocator(2))
            #
            if   i == 2:
                xax2 = axes[i].secondary_xaxis('top', functions=(lon2x_180dy,x2lon_180dy))
            elif i == 5:
                xax2 = axes[i].secondary_xaxis('top', functions=(lon2x_90dy,x2lon_90dy))
            xax2.xaxis.set_ticks(XTicks)
            xax2.xaxis.set_ticklabels(XTickNames)
            xax2.xaxis.set_minor_locator(AutoMinorLocator(2))

        pdf.savefig()
        plt.close(fig)


#===========================================================================
def main():
    taux,lon,lat,_ = \
        Read_In_taux(infile, dx, dy, xstrt=30., xend=130., ystrt=-20., yend=20.)

    #   Add one grid to the south (for presentation sake)
    taux = np.pad(taux, pad_width=((0, 0), (1, 0), (0, 0)), mode="edge")
    lat = np.insert(lat, 0, lat[0]-(lat[1]-lat[0]))

    taux_180dy, taux_90dy = \
        Time_Filtering(T_target, order_filter, sampling_freq, taux, lon, lat)

    regr = Regression_Analysis(taux_180dy, taux_90dy, lon, lat,)

    xfunc_180dy, yfunc_180dy, xfunc_90dy, yfunc_90dy, \
                func_model_180dy, func_model_90dy = \
        Generate_Idealized_Forcing(lon, lat,)

    yu = np.arange(402) *0.1  -20.0   #  box_bounded
    xu = np.arange(622) *0.1  + 0.0

    Draw_Figure(
        OUTFILE_FIG,
        regr,lon,lat,
        xfunc_180dy, yfunc_180dy, xfunc_90dy, yfunc_90dy,
        func_model_180dy, func_model_90dy,
        T_target, xu,yu,
        )


if __name__ == "__main__":
    sys.exit(main())
