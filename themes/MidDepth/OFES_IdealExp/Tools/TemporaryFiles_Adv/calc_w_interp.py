"""========================================================================================
    Interpolate w onto u grids
=========================================================================================="""

import numpy as np
import netCDF4
from scipy.interpolate import RegularGridInterpolator
import sys

#==========================================================================================
#   Read in Input time series
#==========================================================================================

ext_out = sys.argv[1]  # first argument after script name
#print("Input keyword is:", ext_out)

in_u = 'u_eq'+ext_out
in_w = 'w_0p1S-0p1N'+ext_out

f_u = netCDF4.Dataset(in_u, 'r')
f_w = netCDF4.Dataset(in_w, 'r')

xu = np.copy( f_u.variables['lon'] )
yu = np.copy( f_u.variables['lat'] )
zu = np.copy( f_u.variables['lev'] )
tu =          f_u.variables['time']

xw = np.copy( f_w.variables['lon'] )
yw = np.copy( f_w.variables['lat'] )
zw = np.copy( f_w.variables['lev'] )

wvel0 = np.copy( f_w.variables['w'] )

wvel0[np.abs(wvel0) > 1e10] = np.nan

#==========================================================================================
#            Interpolate w onto u grids
#==========================================================================================

zw1   = np.insert(zw, 0, 0.0)                              #  insert the top level
wvel1 = np.concatenate((wvel0[:,0:1,:,:], wvel0), axis=1)

xeast = xw[-1] + (xw[-1] - xw[-2])                         #  xw at len(xw)+1
xw1   = np.insert(xw, len(xw), xeast)                      #  insert the easternmost grid
wvel1 = np.concatenate((wvel1, wvel1[:,:,:,-1:]), axis=3)

t = np.arange(len(tu))

interp1 = RegularGridInterpolator( (t,zw1,yw,xw1), wvel1 )
t1g,z1g,y1g,x1g = np.meshgrid(t,zu,yu,xu, indexing='ij')

wvel = interp1((t1g,z1g,y1g,x1g))

#==========================================================================================
#            Write out to NetCDF files
#==========================================================================================

outfile = 'w_onU_eq'+ext_out

def write_netcdf(file_name, var_name, var_unit_name, var, lon,lat,lev,time):

    print('  calc_w_interp: OUT >>> '+file_name)
    f1 = netCDF4.Dataset(file_name,'w',format='NETCDF4')

    f1.createDimension('lon', len(lon))
    f1.createDimension('lat', len(lat))
    f1.createDimension('lev', len(lev))
    f1.createDimension('time',len(time))

    lon1  = f1.createVariable('lon'   ,'float64',('lon'))
    lat1  = f1.createVariable('lat'   ,'float64',('lat'))
    lev1  = f1.createVariable('lev'   ,'float64',('lev'))
    time1 = f1.createVariable('time'  ,'float64',('time'))
    var1  = f1.createVariable(var_name,'float32',('time','lev','lat','lon'))

    time1.setncattr('units',time.units)
    var1.setncattr( 'units',var_unit_name)

    lon1[:]  = lon
    lat1[:]  = lat
    lev1[:]  = lev
    time1[:] = np.copy(time)
    var1[:]  = var

    f1.close()

write_netcdf(outfile, 'w','cm s^-1', wvel, xu,yu,zu,tu)

