"""
    Plot temperture and N^2 profiles in the initial conditions and last model year
"""
import numpy as np
import netCDF4, gsw, string, sys, pathlib, os

import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from matplotlib.backends.backend_pdf import PdfPages

from common.new_page import (new_page,)

exp_name = "180dy+90dy"
year_FC = "2029"

MY_WORK_ROOT = str( pathlib.Path(os.environ["MY_WORK_ROOT"]) )
infile_IC = MY_WORK_ROOT+"/code/themes/MidDepth/OFES_IdealExp/Model/"+\
            "InitialConditions/1D/TS/ConstantN/out_N2_5e-6.nc"

OUTPUT_FILENAME = 'fig.pdf'

def main():

    #  ---------------------  Read in T/S in i.c.  --------------------

    f1 = netCDF4.Dataset(infile_IC, 'r')

    f2 = netCDF4.Dataset(    "./temp_mean_"+year_FC+"_"+exp_name+".nc", 'r')
    f3 = netCDF4.Dataset("./salinity_mean_"+year_FC+"_"+exp_name+".nc", 'r')

    T_IC  = np.copy(f1.variables['temp'])
    S_IC  = np.copy(f1.variables['salt'])
    zt_IC = np.copy(f1.variables['lev'])

    T_FC  = np.copy(f2.variables['temp'])
    S_FC  = np.copy(f3.variables['salinity'])
    zt_FC = np.copy(f2.variables['lev'])

    f1.close()
    f2.close()
    f3.close()

    T_IC, T_FC = np.squeeze(T_IC), np.squeeze(T_FC)
    S_IC, S_FC = np.squeeze(S_IC), np.squeeze(S_FC)

    #  ---------------------  Compute N^2  --------------------

    lon0, lat0 = 0., 0.   #   dummy longitude & latitude

    p_IC  = gsw.conversions.p_from_z(zt_IC*(-1.),lat0)
    p_FC  = gsw.conversions.p_from_z(zt_FC*(-1.),lat0)
    sa_IC = gsw.conversions.SA_from_SP(S_IC,p_IC,lon0,lat0)
    sa_FC = gsw.conversions.SA_from_SP(S_FC,p_FC,lon0,lat0)
    ct_IC = gsw.conversions.CT_from_pt(sa_IC,T_IC)
    ct_FC = gsw.conversions.CT_from_pt(sa_FC,T_FC)
    N2_IC, pmid = gsw.Nsquared(sa_IC, ct_IC, p_IC, lat=None)
    N2_FC, pmid = gsw.Nsquared(sa_FC, ct_FC, p_FC, lat=None)

    #  ---------------------  Open PDF  --------------------

    txt_labels = [f"({c})" for c in string.ascii_lowercase]

    with PdfPages(OUTPUT_FILENAME) as pdf:

        fig, axes, iplot = new_page(nrows=1, ncols=2, figsize=(7,5))

        ax = axes[iplot]
        ax.set_title(txt_labels[iplot])
        ax.plot(T_IC, zt_IC, 'k-', label="Initial Condition")
        ax.plot(T_FC, zt_FC, 'r-', label="Year "+year_FC)
        ax.set(xlabel='Potential Temperature (\N{degree sign}C)', \
               ylabel='Depth (m)')
        iplot += 1

        ax = axes[iplot]
        ax.set_title(txt_labels[iplot])
        ax.plot(N2_IC, pmid, 'k-', label="Initial Condition")
        ax.plot(N2_FC, pmid, 'r-', label="Year "+year_FC)
        ax.set(xlabel=r'N$^2$ (s$^{-2}$)', ylabel='Depth (m)')
        ax.set_xlim(xmin=0.)
        iplot += 1

        for ax in axes:
            ax.set_ylim([5000.,0.])

            #ax.set_title('exp_name = '+exp_name)
            #ax.legend()
            ax.xaxis.set_minor_locator(AutoMinorLocator(4))
            ax.yaxis.set_minor_locator(AutoMinorLocator(4))
            ax.xaxis.set_ticks_position('both')
            ax.yaxis.set_ticks_position('both')

        plt.tight_layout()
        pdf.savefig()
        plt.close(fig)

    #---------

if __name__ == "__main__":
    sys.exit(main())
