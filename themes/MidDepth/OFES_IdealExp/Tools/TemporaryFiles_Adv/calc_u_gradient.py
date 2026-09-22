"""========================================================================================
    Compute du/dx, du/dy, and du/dz
=========================================================================================="""
import numpy as np
import netCDF4
import sys

#==========================================================================================

#==========================================================================================
#   Constants
#==========================================================================================

deg2cm = 1.1132387e7    #   cm/deg for longitude (from L2643 in OFES_patch/topog.F)
#OneDayInSec = 86400.    #   One day in sec
#dt = 1.                 #   Sampling rate [day]

#==========================================================================================
#   Read in Input time series
#==========================================================================================

ext_out = sys.argv[1]  # first argument after script name
#print("Input keyword is:", ext_out)

in_u = 'u_0p1S-0p1N'+ext_out

f1 = netCDF4.Dataset(in_u, 'r')
uvel = np.copy(f1.variables['u'])
lon  = np.copy(f1.variables['lon'] )
lat  = np.copy(f1.variables['lat'] )
zu   = np.copy(f1.variables['lev'] )
time =         f1.variables['time']

uvel[np.abs(uvel) > 1e10] = np.nan    #   -9.99999e33 --> NaN

uvel_saved = np.copy(uvel)

uvel = uvel[:,:,:,:-1]    #  Cut out Nan at the eastern end grid
                          #  Dimension of longitude of uvel is decreased by 1.

if np.any(np.isnan(uvel)):
    raise ValueError("uvel includes NaN")

#==========================================================================================
#            Compute gradients
#==========================================================================================

jcent = 1

dx_cm = (lon[2]-lon[0]) *deg2cm *np.cos( np.deg2rad(lat[1]) )
dy_cm = (lat[2]-lat[0]) *deg2cm
dz_cm = np.abs( (zu[0:-2]-zu[2:]) *1e2 )  #  *1e2 is for cm

dx_cm_half_west = (lon[ 1]-lon[ 0]) *deg2cm *np.cos( np.deg2rad(lat[1]) )
dx_cm_half_east = (lon[-1]-lon[-2]) *deg2cm *np.cos( np.deg2rad(lat[1]) )

dz_cm_half_top = np.abs( (zu[ 0]-zu[ 1]) *1e2 )  #  *1e2 is for cm
dz_cm_half_bot = np.abs( (zu[-2]-zu[-1]) *1e2 )

#print('dz_cm = ',dz_cm)

#---------------------
#   du/dx

dudx = (uvel[:,:,jcent:jcent+1, 2:] -
        uvel[:,:,jcent:jcent+1, 0:-2]) /dx_cm      #  "jcent:jcent+1" is to keep the lat axis

dudx_west = ( 0.5 *(uvel[:,:,jcent:jcent+1, 1]
                  + uvel[:,:,jcent:jcent+1, 0]) - 0. ) /dx_cm_half_west

dudx_east = ( 0. - 0.5*(uvel[:,:,jcent:jcent+1, -1]
                      + uvel[:,:,jcent:jcent+1, -2]) ) /dx_cm_half_east

dudx = np.concatenate((dudx_west[:,:,:,None], dudx), axis=3)  # Insert to western end
dudx = np.concatenate([dudx, dudx_east[:,:,:,None]], axis=3)  # Insert to eastern end

#---------------------
#   du/dy

dudy = (uvel[:,:,2:3,:] -
        uvel[:,:,0:1,:]) /dy_cm     #  "2:3" and "0:1" are to keep the lat axis

#---------------------
#   du/dz

dudz = (uvel[:, 0:-2, jcent:jcent+1,:] -
        uvel[:, 2:,   jcent:jcent+1,:]) / dz_cm[None, :, None, None]
                                #  "jcent:jcent+1" is to keep the lat axis

dudz_top = ( 0. - 0.5 *(uvel[:, 0,jcent:jcent+1,:]
                      + uvel[:, 1,jcent:jcent+1,:]) ) /dz_cm_half_top

dudz_bot = ( 0.5 *(uvel[:,-2,jcent:jcent+1,:]
                 + uvel[:,-1,jcent:jcent+1,:]) - 0. ) /dz_cm_half_bot

dudz = np.concatenate((dudz_top[:,None,:,:], dudz), axis=1)  #  Insert to top end
dudz = np.concatenate((dudz, dudz_bot[:,None,:,:]), axis=1)  #  Insert to bottom end

#---------------------
#    Add Nan to the eastern end such that the shapes of uvel[:,:,1,:], dudx,
#        dudy, and  are identical.
#    Dimension of longitude of dudx, dudy, dudz is increased by 1.

pads = ((0, 0), (0, 0), (0, 0), (0, 1))
dudx = np.pad( dudx, pads, constant_values=np.nan )
dudy = np.pad( dudy, pads, constant_values=np.nan )
dudz = np.pad( dudz, pads, constant_values=np.nan )

#uvel = uvel_saved
#print("uvel.shape = ", uvel.shape)
#print("dudx.shape = ", dudx.shape)
#print("dudy.shape = ", dudy.shape)
#print("dudz.shape = ", dudz.shape)
#print("np.sum(np.isnan(uvel)) = ", np.sum(np.isnan(uvel)))
#print("np.sum(np.isnan(dudx)) = ", np.sum(np.isnan(dudx)))
#print("np.sum(np.isnan(dudy)) = ", np.sum(np.isnan(dudy)))
#print("np.sum(np.isnan(dudz)) = ", np.sum(np.isnan(dudz)))
#print("uvel[0,0,1,:5],uvel[0,0,1,-5:] = ", uvel[0,0,1,:5],uvel[0,0,1,-5:])
#print("dudx[0,0,0,:5],dudx[0,0,0,-5:] = ", dudx[0,0,0,:5],dudx[0,0,0,-5:])
#print("dudy[0,0,0,:5],dudy[0,0,0,-5:] = ", dudy[0,0,0,:5],dudy[0,0,0,-5:])
#print("dudz[0,0,0,:5],dudz[0,0,0,-5:] = ", dudz[0,0,0,:5],dudz[0,0,0,-5:])

#==========================================================================================
#            Write out to NetCDF files
#==========================================================================================

outfile_ux = 'dudx_eq'+ext_out
outfile_uy = 'dudy_eq'+ext_out
outfile_uz = 'dudz_eq'+ext_out

def write_netcdf(file_name, var_name, var_unit_name, var, lon,lat,lev,time):

    print('  calc_u_gradient: OUT >>> '+file_name)
    f1 = netCDF4.Dataset(file_name,'w',format='NETCDF4')

    f1.createDimension('lon', len(lon))
    f1.createDimension('lat', len(lat))
    f1.createDimension('lev', len(zu))
    f1.createDimension('time',len(time))

    lon1  = f1.createVariable('lon'    ,'float32',('lon'))
    lat1  = f1.createVariable('lat'    ,'float32',('lat'))
    lev1  = f1.createVariable('lev'    ,'float32',('lev'))
    time1 = f1.createVariable('time'   ,'float64',('time'))
    var1  = f1.createVariable(var_name ,'float32',('time','lev','lat','lon'))

    time1.setncattr('units',time.units)
    var1.setncattr( 'units',var_unit_name)

    lon1[:]  = lon
    lat1[:]  = lat
    lev1[:]  = zu
    time1[:] = np.copy(time)
    var1[:,:,:,:] = var

    f1.close()

write_netcdf(outfile_ux, 'dudx','s^-1', dudx, lon,lat[jcent:jcent+1],zu,time)
write_netcdf(outfile_uy, 'dudy','s^-1', dudy, lon,lat[jcent:jcent+1],zu,time)
write_netcdf(outfile_uz, 'dudz','s^-1', dudz, lon,lat[jcent:jcent+1],zu,time)

