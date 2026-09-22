import numpy as np
import netCDF4

def Read_In_taux(infile, dx, dy, xstrt=40., xend=105., ystrt=-20., yend=20.):
    import numpy as np
    import netCDF4

    print("IN << "+infile)

    f1 = netCDF4.Dataset(infile, 'r')
    taux = np.copy(f1.variables['taux'])
    lon  = np.copy(f1.variables['lon'])
    lat  = np.copy(f1.variables['lat'])
    time =         f1.variables['time']

    #   Missing value to NaN

    taux = np.where( np.abs(taux) > 1e10, np.nan, taux )

    #  Round off up to one decimal place, 0.1

    lon = np.round(lon, 1)
    lat = np.round(lat, 1)

    #   Limit the region longitudinally

    i = np.argwhere(np.logical_and( lon >= xstrt, lon <= xend))[:,0]
    taux = taux[:,:,i]
    lon = lon[i]

    #   Limit the region latitudinally

    j = np.argwhere(np.logical_and( lat >= ystrt, lat <= yend ))[:,0]
    taux = taux[:,j,:]
    lat = lat[j]

    #   Box averaging in the zonal direction

    if dx != 0:
        nx1 = (len(lon) //dx) *dx

        tmp1 = np.array(taux[:,:,0:nx1]).reshape(len(time),len(lat),-1,dx)
        taux = np.nanmean(tmp1,axis=3)

        tmp = lon[0:nx1].reshape(-1, dx)

        lon = lon[0:nx1].reshape(-1, dx).mean(axis=1)

    #   Box averaging in the meridional direction

    if dy != 0:
        ny1 = (len(lat) //dy) *dy

        tmp1 = np.array(taux[:,0:ny1,:]).reshape(len(time),-1,dy,len(lon))
        taux = np.nanmean(tmp1,axis=2)

        tmp = lat[0:ny1].reshape(-1, dy)

        lat = lat[0:ny1].reshape(-1, dy).mean(axis=1)

    #print('np.shape(taux) = ',np.shape(taux))
    #print('np.shape(lon)  = ',np.shape(lon))
    #print('np.shape(lat)  = ',np.shape(lat))
    #print('np.shape(time) = ',np.shape(time))

    #   Return
    return taux,lon,lat,time
