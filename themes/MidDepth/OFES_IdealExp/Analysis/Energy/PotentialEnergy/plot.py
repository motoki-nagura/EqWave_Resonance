import numpy as np
import netCDF4
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

#-----------------------------------------------------------

f1 = netCDF4.Dataset('rho_diff_sq_hmean_merged.nc', 'r')
lev  = np.copy(f1.variables['lev'])
time = np.copy(f1.variables['time'])
var  = np.copy(f1.variables['rho'])
f1.close()

var = np.squeeze(var)

time = time - time[0]  #    First day is hour 0
time /= 24.            #    hours into days

var_zavg = np.average(var, axis=1)

#-----------------------------------------------------------

fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(8,8))

ttl1     = r'An index to PE, $\rho^{\prime 2}$'
ttl2     = r'Anomaly of the index to PE, $\rho^{\prime 2}$'
ttl_unit = r'(g$^2$ cm$^{-6}$)'

iplot=0
axs[iplot].set_title(ttl1)
axs[iplot].plot(time,var_zavg,'k-')
axs[iplot].set(xlabel='Time (days)', ylabel=ttl_unit)

iplot=1
for iplot in range(1,3):
    axs[iplot].set_title(ttl2)
    cs=axs[iplot].contourf(time,lev,np.transpose(var),np.arange(0,1.1,0.1)*1e-11, \
                           cmap='jet',extend='both')
    cbar = plt.colorbar(cs,ax=axs[iplot])
    cbar.set_label(ttl_unit)
    axs[iplot].set(xlabel='Time (days)', ylabel='Depth (m)')
    axs[iplot].invert_yaxis()
    if iplot == 2:
        axs[iplot].set_xlim([0.,1000.])

for ax in axs:
    ax.xaxis.set_minor_locator(AutoMinorLocator(4))

plt.tight_layout()
plt.savefig('fig.pdf')

