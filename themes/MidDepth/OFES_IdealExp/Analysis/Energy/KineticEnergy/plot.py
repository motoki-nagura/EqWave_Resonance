import numpy as np
import netCDF4
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

#-----------------------------------------------------------

f1 = netCDF4.Dataset('sum_mean_merged.nc', 'r')
lev  = np.copy(f1.variables['lev'])
time = np.copy(f1.variables['time'])
var  = np.copy(f1.variables['u'])
f1.close()

var = np.squeeze(var)

time = time - time[0]  #    First day is hour 0
time /= 24.            #    hours into days

var_zavg = np.average(var, axis=1)

#-----------------------------------------------------------

fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(8,8))

iplot=0
axs[iplot].set_title(r'Mean of KE, $(u^2 + v^2)/2$')
axs[iplot].plot(time,var_zavg,'k-')
axs[iplot].set(xlabel='Time (days)', ylabel=r'KE (cm$^2$ s$^{-2}$)')

iplot=1
for iplot in range(1,3):
    axs[iplot].set_title(r'Horizontal Mean of KE, $(u^2 + v^2)/2$')
    cs=axs[iplot].contourf(time,lev,np.transpose(var),np.arange(0,2.1,0.2),cmap='jet',extend='both')
    cbar = plt.colorbar(cs,ax=axs[iplot])
    cbar.set_label('cm$^2$ s$^{-2}$')
    axs[iplot].set(xlabel='Time (days)', ylabel='Depth (m)')
    axs[iplot].invert_yaxis()
    if iplot == 2:
        axs[iplot].set_xlim([0.,1000.])

for ax in axs:
    ax.xaxis.set_minor_locator(AutoMinorLocator(4))

plt.tight_layout()
plt.savefig('fig.pdf')

