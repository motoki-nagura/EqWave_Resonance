import numpy as np
from config import Params

def get_params() -> Params:
    #--------
    day_in_sec = 86400.           #  1 day in seconds
    π = np.pi
    #--------
    return Params(
        day_in_sec = day_in_sec,
        N2         = 5e-6,                  #  Buoyancy frequency squared, N^2 [s^-2]
        n1_values = np.array([1]),          #  Meridional mode number for wave 1
        n2_values = np.array([-1,0,1,2,3]), #                             wave 2
        n3_values = np.array([-1,0,1,2,3]), #                             wave 3
        ω1 = 2.*π /(180. *day_in_sec),      #  Frequency for wave 1 [s^-1]. 180 days.
        ω2 = 2.*π /(180. *day_in_sec),      #  Frequency for wave 2 [s^-1]. 180 days.
        ω3 = 2.*π /( 90. *day_in_sec),      #  Frequency for wave 3 [s^-1]. 90 days.  (Must be highest |ω|)
        k1_values  = np.linspace(0.0e-5, 0.5e-5, 61),
                                            #  Zonal wavelength for wave 1 [rad m^-1]
        m1_guess   = -5.0e-3,               #  Initial guess of zonal wavenumber for wave 1 [rad m^-1]
        k2_guesses = 2.*π *np.linspace(-1./400e3, 1./400e3, 150),
                                            #  Initial guesses of zonal wavenumber for wave 2 [rad m^-1]
        m2_guesses = 2.*π *np.linspace(-1./250. , 1./250.,  150),
                                            #  Initial guesses of vertical wavenumber for wave 2 [rad m^-1]
    )
