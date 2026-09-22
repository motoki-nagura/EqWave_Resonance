"""
    Compute T and S corresponding to a constant N^2
"""
from rich import traceback
traceback.install()

import numpy as np
import netCDF4, gsw
import matplotlib.pyplot as plt
import sys

from themes.MidDepth.OFES_IdealExp.Tools.common.read_grid_dta import (
        read_grid_dta,)
#from tools.ofes_exp.ReadGrids.python.calc import read_grid_dta

#=================================================================
#             Parameters
#=================================================================

S0            = 34.5     #   The value of constant salinity [psu]

Ttop_init     = 25.      #   Temperature at the top level [Cdeg]
max_iteration = 100      #   Maximum number of iteration
tolerance     = 1e-9     #   Tolerance [s^-1]

mixed_layer_depth = -1.  #   Mixed layer depth where T = Ttop_init [m]
#mixed_layer_depth = 50.  #   Mixed layer depth where T = Ttop_init [m]
#mixed_layer_depth = 500. #   No mixed layer for negative "mixed_layer_depth"

grav          = 9.81     #  Gravitational acceleration rate [m s^-2]

N2_target     = 5.0e-6   #  Prescribe the background buoyancy frequency, N^2 [s^-2]
#N2_target     = 1.0e-5

#=================================================================
#             In / Out
#=================================================================
infile_grid="/S/data01/G4006/y0330/data/ofes_exp/settings/box_bounded/topog/grid.dta.out"

s = f"{N2_target:.0e}".replace('e-0', 'e-').replace('e+0', 'e+')
outfile_fig = f"fig_N2_"+s+".pdf"
outfile_nc  = f"out_N2_"+s+".nc"

#----------------------------------------------------------------
#      Function to compute N2

def calc_N2(temp1,salt1,z1,lon1,lat1):
    p  = gsw.p_from_z(z1, lat1, geo_strf_dyn_height=0, sea_surface_geopotential=0)
    SA = gsw.SA_from_SP(salt1, p, lon1, lat1)
    CT = gsw.CT_from_pt(SA, temp1)
    N2, pmid = gsw.stability.Nsquared(SA, CT, p, lat=lat1)
    return N2, pmid

def main():
    #=================================================================
    #             Read in lev
    #=================================================================
    dxt, dyt, dxu, dyu, dzt, dzw, xt, xu, yt, yu, zt, zw = read_grid_dta(infile_grid)

    lev = np.array(zt) *1e-2  #  cm --> m

    lon, lat = 0., 0.

    z = lev *(-1.)  #   Vertical coordinate

    salt = np.ones(len(lev)) *S0   #   salinity = constant

    #=================================================================
    #    Compute temp such that the resulting N2 = N2_target
    #=================================================================

    Ttop = Ttop_init

    #----------------------------------------------------------------
    #           Iteration

    temp = np.zeros(len(lev))

    if mixed_layer_depth < 0.:
        k_range = range(len(lev)-1)
        temp[0] = Ttop                             #   Set temperature in the mixed layer depth

    else:
        kml = np.argwhere(lev <= mixed_layer_depth)[:,0]
        temp[kml] = Ttop                           #   Set temperature in the mixed layer depth

        k_range = range(np.max(kml), len(lev)-1, 1)

    for k in k_range:                              #   Depth loop

        #print('')
        #print('k, lev[k] = ',k,lev[k])

        T_set = np.array([Ttop, Ttop-1.])          #   Initial estimate
        salt_tmp, z_tmp = salt[k:k+2], z[k:k+2]

        for i in range(max_iteration):             #   Iteration Loop
            N2_1, _ = calc_N2(np.array([Ttop, T_set[0]]), salt_tmp,z_tmp,lon,lat)
            N2_2, _ = calc_N2(np.array([Ttop, T_set[1]]), salt_tmp,z_tmp,lon,lat)

            diff1, diff2 = N2_1-N2_target, N2_2-N2_target   #   Compute two estimates

            if diff1*diff2 > 0:
                print('k = ',k,', i = ',i,', diff1 * diff2 > 0')
                                                   #   Must be N2_1 < N2_target < N2_2 or
                                                   #           N2_2 < N2_target < N2_1

            if   np.abs(diff1) > np.abs(diff2):    #   Discard the one with larger error
                T_set[0] = np.average(T_set)       #     and update it with the mean
            elif np.abs(diff1) < np.abs(diff2):    #     of the two estimates
                T_set[1] = np.average(T_set)

            #print('i, diff1, diff2 = ',i,diff1,diff2)

            if np.logical_and(np.abs(diff1) < tolerance, np.abs(diff2) < tolerance):
                break                                       #   Exit

        T_final = np.average(T_set)                         #   The result

        #N2_final, tmp = calc_N2([Ttop,T_final], salt_tmp,z_tmp,lon,lat)
        #print('N2_final  = ',N2_final)
        #print('N2_target = ',N2_target)
        #print('T_final = ',T_final)

        temp[k+1] = T_final                                 #   Save
        Ttop = T_final                                     #   Update Ttop

    #print('')
    #print('temp = ',temp)

    N2_res, p_mid = calc_N2(temp,salt,z,lon,lat)

    #=================================================================
    #                  Output to a NetCDF file
    #=================================================================

    N2_target_array = np.ones(len(lev)) *N2_target   #  create output array

    f1 = netCDF4.Dataset(outfile_nc,'w',format='NETCDF4')
    f1.createDimension('lev',   len(lev))
    lev_1 = f1.createVariable('lev', 'float32',('lev'))
    T_1   = f1.createVariable('temp','float32',('lev'))
    S_1   = f1.createVariable('salt','float32',('lev'))
    n2_1  = f1.createVariable('N2',  'float32',('lev'))
    T_1.setncattr( 'units','Cdeg')
    S_1.setncattr( 'units','psu')
    n2_1.setncattr('units','s^-2')
    lev_1[:] = lev
    T_1[:]   = temp
    S_1[:]   = salt
    n2_1[:]  = N2_target_array
    f1.close()

    ##=================================================================
    ##                 Fit a linear function to temp
    ##=================================================================
    #
    #from scipy.optimize import curve_fit
    #
    #def func(x, a, b):
    #    return a*x + b
    #
    #popt, pcov = curve_fit(func, z, temp)
    #
    #T_fit = func(z, *popt)
    #
    #N2_fit, p_mid = calc_N2(T_fit,salt,z,lon,lat)

    #=================================================================
    #                       Plots
    #=================================================================

    fig, axs0 = plt.subplots(nrows=2, ncols=2, figsize=(8,8))
    arr_ij = [(x, y) for x, y in np.ndindex(axs0.shape)]
    axs = [axs0[index] for index in arr_ij]

    axs[0].plot( salt,              lev, 'k-', label='S')
    axs[1].plot( temp,              lev, 'k-', label='T')
    axs[2].plot( N2_target_array,   lev, 'r:', label=r'$N^2$, target')
    axs[2].plot( N2_res,          p_mid, 'k-', label=r'$N^2$')
    axs[2].set_xlim([N2_target*0.99, N2_target*1.01])
    axs[0].set(xlabel='psu')
    axs[1].set(xlabel='Cdeg')
    axs[2].set(xlabel=r'$s^{-2}$')
    fig.delaxes(axs[3])
    for ax in [axs[0],axs[1],axs[2]]:
        ax.legend()
        ax.invert_yaxis()
        ax.set(ylabel='Depth (m)')
    plt.tight_layout()
    plt.savefig(outfile_fig)

    #axs[1].plot( T_fit  , lev  , 'r:', label='T, fitted')
    #axs[2].plot(N2_fit  , p_mid, 'r:', label=r'$N^2$, fitted')

if __name__ == "__main__":
    sys.exit(main())
