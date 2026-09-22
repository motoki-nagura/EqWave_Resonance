import numpy as np
from scipy.optimize import root

def DispRel_mwvn(N2,β, ω, kwvn, nmode, m_initial=0.0): #{{{
    """
-------------------------------------------------------------------------------------
    Solve the dispersion relation for m using scipy.optimize.root
-------------------------------------------------------------------------------------
  Input
     N2        (scalar) : Buoyancy frequency squared [s^-2]
     β         (scalar) : Meridional gradient of the Coriolis coefficient [m^-1 s^-1]
     ω         (scalar) : Frequency [rad s^-1]
     kwvn      (scalar) : Zonal wavenumber [rad m^-1]
     nmode     (scalar) : Meridional mode number
     m_initial (scalar) : Initial guess of vertical wavenumber [rad m^-1]

  Output
     mwvn      (scalar) :  Vertical wavenumber [rad m^-1].
-------------------------------------------------------------------------------------
    """

    def func_NonNegative_n(m, β,N2,k):
        F1 = m**2 - β *np.sqrt(N2) *(2.*nmode+1.) *np.abs(m) /ω**2 \
                  - N2 *k /ω**2 *(β /ω + k)
        return F1
    #---
    if ω < 0.:
        raise RuntimeError("DispRel_mwvn: ω must be positive.   ω = ",ω)

    #---
    if np.logical_or( np.logical_not( np.isscalar(kwvn ) ),
                      np.logical_not( np.isscalar(nmode) ) ):
        raise RuntimeError("DispRel_mwvn: This subroutine was edited, and kwvn and nmode must be scalars. Input kwvn or nmode (or both) is an array. Fix it.")

    #---
    mwvn = np.nan

    if   nmode == 0.:    #   Yanai mode
        tmp_m = np.sqrt(N2) /ω *(kwvn + β/ω)
        if tmp_m < 0.:
            tmp_m = np.nan        #   |m| must be positive
        mwvn = tmp_m

    elif nmode >= 1.:    #   Rossby mode
        if kwvn <= 0.:            #   ω > 0. k must be negative when n >= 1
            res = root(func_NonNegative_n, m_initial, args=(β,N2,kwvn), method="lm", tol=1e-12)
            mwvn = res.x[0]

            if not res.success:
                print("  Success:",res.success,", n=",nmode,", k=",kwvn)
                print(res)
                mwvn = np.nan

    elif nmode == -1:    #   Kelvin mode
        if kwvn >= 0.:            #   ω > 0. k must be positive when n == -1
            mwvn = np.sqrt(N2) *kwvn /ω

    #---
    return mwvn
#}}}

#==========================================================================================
def DispRel_mwvn_ver2(N2,β, ω, kwvn, nmode):  #{{{
    """
-------------------------------------------------------------------------------------
    Solve the dispersion relation for m using numpy.roots
-------------------------------------------------------------------------------------
  Input
     N2     (scalar) : Buoyancy frequency squared [s^-2]
     β      (scalar) : Meridional gradient of the Coriolis coefficient [m^-1 s^-1]
     ω      (scalar) : Frequency [rad s^-1]
     kwvn   (scalar) : Zonal wavenumber [rad m^-1]
     nmode  (scalar) : Meridional mode number

  Output
     mwvn   (scalar) :  Vertical wavenumber [rad m^-1].
-------------------------------------------------------------------------------------
    """

    #---
    if ω < 0.:
        raise RuntimeError("DispRel_mwvn_ver2: ω must be positive.   ω = ",ω)

    #---
    if np.logical_or( np.logical_not( np.isscalar(kwvn ) ),
                      np.logical_not( np.isscalar(nmode) ) ):
        raise RuntimeError("DispRel_mwvn_ver2: This subroutine was edited, and kwvn and nmode must be scalars. Input kwvn or nmode (or both) is an array. Fix it.")

    #---
    mwvn = np.zeros(2) *np.nan

    if nmode == 0.:     #   Yanai mode
        mwvn[0] = np.sqrt(N2) /ω *(kwvn + β/ω)
        mwvn[1] = np.nan

    if nmode >= 1.:     #   Rossby mode
        if kwvn <= 0.:                  #   ω > 0. k must be negative when n >= 1
            p = np.array([ω**2 /N2, \
                          - β *(2.*nmode+1.) /np.sqrt(N2), \
                          - β *kwvn /ω - kwvn**2])
            tmp_m = np.roots(p)
            tmp_m = tmp_m[ np.argsort( np.abs(tmp_m) ) ]
            tmp_m[ tmp_m < 0. ] = np.nan #   m must be positive
            mwvn = tmp_m

    elif nmode == -1:   #   Kelvin mode
        if kwvn >= 0.:                  #   ω > 0. k must be positive when n == -1
            mwvn[0] = np.sqrt(N2) *kwvn /ω
            mwvn[1] = np.nan

    #---
    return mwvn
#}}}

