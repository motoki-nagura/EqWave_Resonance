"""
    Plot the time series of kinetic and available potential energy averaged over the basin
"""

import numpy as np
import netCDF4
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import AutoMinorLocator

#-----------------------------------------------------------

f1 = netCDF4.Dataset('../KineticEnergy/sum_mean_merged.nc', 'r')
lev  = np.copy(f1.variables['lev'])
#time = np.copy(f1.variables['time'])
u2   = np.copy(f1.variables['u'])
f1.close()

f1 = netCDF4.Dataset('../PotentialEnergy/rho_diff_sq_hmean_merged.nc', 'r')
ρ2  = np.copy(f1.variables['rho'])
f1.close()

u2, ρ2 = np.squeeze(u2), np.squeeze(ρ2)

time = np.arange(0, len(u2))   #  time in days

#----------------------------------------------------------
#            Unit Conversion

u2  *= 1e-4                  #  u² + v²,  cm^2 s^-2 ---> m^2 s^-2

ρ0   = 1.035 *1e3            #  Mean sea water density [kg m^-3], From MOM3 manual
grav = 980.6 *1e-2           #  Acceleration due to gravity [m s^-2], From MOM3 manual
N2   = 5.0e-6                #  Buoyancy frequency squared [s^-2],
                             #      From Model/InitialConditions/1D/TS/ConstantN/calc.py

ρ2   *= 1e6                  #  ρ²,  g cm^-3 ---> kg m^-3
b2    = ρ2 *grav**2 /ρ0**2   #  Buoyancy squared [m^2 s^-4]
b2N2  = b2 /N2               #  b²/N² [m^2 s^-2]

u2   *= 0.5                  #  Kinetic energy, (u² + v²)/2 [m^2 s^-2]
b2N2 *= 0.5                  #  Available potential energy, b² /(2N²) [m^2 s^-2]

Ener = u2 + b2N2             #  Total energy. [m^2 s^-2]

Ener *= 1e4                  #  m^2 s^-2 --> cm^2 s^-2
u2   *= 1e4
b2N2 *= 1e4

#----------------------------------------------------------
#             Vertical Average

z0 = 500.

k = np.argwhere(lev > z0)[:][0]
Ener_zavg = np.average(Ener[:,k], axis=1)
u2_zavg   = np.average(u2[  :,k], axis=1)
b2N2_zavg = np.average(b2N2[:,k], axis=1)


#-----------------------------------------------------------

ttl1     = rf'Energy, $(u^2 + v^2)/2 + b^2/(2 N^2)$, below {int(z0):d} m'
ttl2     = rf'Energy, $(u^2 + v^2)/2 + b^2/(2 N^2)$'

#scl_exp  = int(np.log10(scale))
#ttl_unit = fr'$10^{{{scl_exp}}}$'+r' m$^2$ s$^{-2}$'

ttl_unit = r' cm$^2$ s$^{-2}$'


with PdfPages('fig.pdf') as pdf:

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8,3))

    ax.plot(time,Ener_zavg, 'k-', label='E')
    ax.plot(time,  u2_zavg, 'r-', label='KE')
    ax.plot(time,b2N2_zavg, 'b-', label='APE')
    ax.legend()
    ax.set(xlabel='Time (days)', ylabel=ttl_unit)
    ax.xaxis.set_minor_locator(AutoMinorLocator(4))
    ax.yaxis.set_ticks(np.arange(0., 3.1, 1.0))
    ax.yaxis.set_minor_locator(AutoMinorLocator(5))

    plt.tight_layout()
    pdf.savefig()
    plt.close(fig)

    #----------------

    for iloop in range(2):

        fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(8,6))
        iplot=0

        #axes[iplot].set_title(ttl1)
        axes[iplot].plot(time,Ener_zavg, 'k-', label='E')
        axes[iplot].plot(time,  u2_zavg, 'r-', label='KE')
        axes[iplot].plot(time,b2N2_zavg, 'b-', label='APE')
        axes[iplot].legend()
        axes[iplot].set(xlabel='Time (days)', ylabel=ttl_unit)
        axes[iplot].yaxis.set_minor_locator(AutoMinorLocator(2))
        iplot += 1

        axes[iplot].set_title(ttl2)
        cs=axes[iplot].contourf(time,lev,np.transpose(Ener), np.arange(0,5.1,0.2), \
                               cmap='jet',extend='both')
        cbar = plt.colorbar(cs,ax=axes[iplot])
        cbar.set_label(ttl_unit)
        axes[iplot].set(xlabel='Time (days)', ylabel='Depth (m)')
        axes[iplot].invert_yaxis()
        axes[iplot].yaxis.set_minor_locator(AutoMinorLocator(5))

        for ax in axes:
            ax.xaxis.set_minor_locator(AutoMinorLocator(4))
            if iloop == 1:
                ax.set_xlim([0.,1000.])

        plt.tight_layout()
        pdf.savefig()
        plt.close(fig)

