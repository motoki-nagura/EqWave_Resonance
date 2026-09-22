import numpy as np

from tools.DispersionRelation.VerticallyPropagating.for_mwvn.sub_DispRel_mwvn import (
        DispRel_mwvn,)
from themes.MidDepth.OFES_IdealExp.Tools.common.Fixed_Parameters import (
        equatorial_beta,)

def Overplot_DispersionRelation(
        omega, kwvns, ax, n_in=[-1,1,3],
        N2=5e-6, line_color='k', line_width=1.0
        ):
    """ Superpose the dispersion relation of equatorial wave modes to a figure
         omega  : Frequency [rad second^-1]
         kwvns  : Zonal wavenumber [rad meter^-1]
         ax     : Matplotlib axis
         n_in   : Meridional mode number (default = [-1,1,3])
         N2     : Background buoyancy frequency squared [s^-2] (default = 5e-6)

         no longer used
         EΩ     : Earth's rotation rate [s^-1] (default = 7.292e-5)
         Erad   : Earth's radius [m] (default = 6371e3)
    """
    nmodes = np.array(n_in)
    β = equatorial_beta()  #  Meridional gradient of the Coriolis coefficient [m^-1 s^-1]
    #β      = 2.*EΩ /Erad
    #
    mwvns = np.zeros((len(kwvns),len(nmodes)))*np.nan
    for ik,k in enumerate(kwvns):
        for inm,n in enumerate(nmodes):
            mwvns[ik,inm] = DispRel_mwvn(N2,β, omega,k,n)  #  Vertical wavenumber [rad meter^-1]
    #
    #      -1   0    1    2   3
    lss = ['-','-','--','-.',':']
    for jn in range(len(nmodes)):
        ax.plot(kwvns, mwvns[:,jn], c=line_color, ls=lss[jn], lw=line_width,
                label=f'n={int(nmodes[jn]):1d}')
    #axs_1.legend()
    for jn in range(len(nmodes)):
        ax.plot(kwvns, mwvns[:,jn]*(-1.), c=line_color, ls=lss[jn], lw=line_width)   #   negative m
