"""   subroutines for rays  """

import numpy as np

from tools.DispersionRelation.VerticallyPropagating.for_mwvn.sub_DispRel_mwvn import DispRel_mwvn
from themes.MidDepth.OFES_IdealExp.Tools.common.Fixed_Parameters import Fixed_Parameters

def ray_K_LongRossby(
        Tper: float, N2 :float, xend :float, zend :float, txt_init :str,
        x0=0., z0=0., nmer=1, nrepeat=4,
        ):
    """
       Compute rays of Kelvin and long Rossby waves considering
          reflection at the lateral and vertical boundaries

        Input

           Tper     : Period [days]
           N2       : Background buoyancy frequency [s^-2]
           xend[2]  : Longitudes of WB and EB [degrees]
           zend[2]  : Depth of the surface and bottom [m]
           txt_init : "K" or "R". The direction of initial propagation
           x0       : Initial location in longitude [degree]
           z0       : Initial location in depth [m]
           nmer     : Meridional mode number of Rossby waves
           nrepeat  : # of reflections

       Output

           xs       :  Zonal end points of rays
           zs       :  Vertical end points of rays
    """
    π = np.pi
    deg2met = Fixed_Parameters.deg2met       #   1 deg in meters
    day2sec = Fixed_Parameters.day2sec       #   1 day in seconds
    #
    ω = 2.*π /(Tper *day2sec)                #   Frequency
    #
    Ktanθ  = - ω /np.sqrt(N2)                #   tanθ for downward Kelvin
    Rtanθ  =   ω /np.sqrt(N2) *(2*nmer+1)    #   tanθ for Rossby
    #
    ver_dir = -1                             #   Vertical propagation index.  -1 = downward
                                             #                                 1 = upward

    if   txt_init == 'K':                    #   Zonal propagation index. 
        zon_dir = 1                          #     1 = Kelvin = eastward
    elif txt_init == 'R':
        zon_dir = -1                         #    -1 = Rossby = westward
    else:
        raise ValueError("Invalid input for txt_init")

    #
    xs, zs = np.zeros((2,nrepeat+1))
    #
    zs[0] = z0                               #   The initial depth
    xs[0] = x0                               #   The initial longitude
    #
    for irepeat in range(nrepeat):
        if    zon_dir == 1:                  #   Eastward propagation
            tanθ = Ktanθ                     #      Kelvin
            xs[irepeat+1] = xend[1]          #      End point is EB
        elif  zon_dir == -1:                 #   Westward propagation
            tanθ = Rtanθ                     #      Rossby
            xs[irepeat+1] = xend[0]          #      End point is WB

        fact = ver_dir *deg2met *tanθ

        zs[irepeat+1] = fact *(xs[irepeat+1] - xs[irepeat]) + zs[irepeat]
                                             #   Compute new z

        if   zs[irepeat+1] >= zend[1]:       #   If hit the bottom
            zs[irepeat+1] = zend[1]          #      New z is at the bottom
            xs[irepeat+1] = (zs[irepeat+1] - zs[irepeat]) /fact + xs[irepeat]
                                             #      New x is re-computed
            ver_dir = ver_dir *(-1.)         #      Flip vertical propagation direction

        elif zs[irepeat+1] <= zend[0]:       #   If hit the surface
            zs[irepeat+1] = zend[0]          #      New z is at the surface
            xs[irepeat+1] = (zs[irepeat+1] - zs[irepeat]) /fact + xs[irepeat]
                                             #      New x is re-computed
            ver_dir = ver_dir *(-1.)         #      Flip vertical propagation direction

        else:                                #   If hit the lateral boundary
            zon_dir = zon_dir *(-1.)         #      Flip zonal propagation direction

    return xs,zs


def ray_Rossby_full(
        Tper :float, N2 :float, kwvn :float,
        nmer=1, x0=0.,x1=10.,z0=0.,
        ):
    """
        Compute rays of Rossby waves using the full (dispersive) dispersion
            relationship of Rossby waves.

         Input

            Tper     : Period [days]
            N2       : Background buoyancy frequency [s^-2]
            kwvn     : Zonal wavenumber [rad meter^-1]
            nmer     : Meridional mode number
            x0       : Initial location in longitude [degree]
            x1       : Final location in longitude [degree]
            z0       : Initial location in depth [m]

        Output

            zend     : Final location in depth [m]
    """
    π = np.pi
    deg2met = Fixed_Parameters.deg2met  #  1 deg in meters
    day2sec = Fixed_Parameters.day2sec  #  1 day in seconds
    β       = Fixed_Parameters.β        #  Meridional gradient of the Coriolis coefficient [m^-1 s^-1]
    #
    ω = 2.*π /(Tper *day2sec)     #  Frequency [rad sec^-1]

    m = DispRel_mwvn(N2,β, ω, kwvn, nmer)  #  Compute Vertical wavenumber [rad meter^-1]

#    m = m[0,0]
    k = kwvn
    n = nmer

    tanθ = -1. *(ω**3 *m**2 + (β + ω *k) *N2 *k) /(N2 *(β + 2.*ω *k) *m)

    fact = -1. *deg2met *tanθ

    zend = fact *(x1 - x0) + z0

    return zend


def draw_arrows_mid_head(axs,xs,zs,clr='k'):
    #
    #        Draw a line and an arrow head at the middle point of the line
    #
    axs.plot(xs, zs, c=clr, ls='-')          #  draw line
    #
    xms = xs[0:-1] + (xs[1:] - xs[0:-1]) / 2.
    zms = zs[0:-1] + (zs[1:] - zs[0:-1]) / 2.
    for xm, ym, x, y in zip(xms,zms,xs,zs):
        axs.annotate(
            "",
            xy=(xm, ym),            # arrow end
            xytext=(x, y),          # arrow start
            arrowprops=dict(
                arrowstyle="->",
                linewidth=1,
                color=clr,
                mutation_scale=20   # ← big arrow head
            )
        )
    return axs
