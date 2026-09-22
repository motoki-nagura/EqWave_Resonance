""" Physical constants and parameter helpers.

    Usage:

      sys.path.insert(0,MY_WORK_ROOT+'/code/themes/MidDepth/OFES_IdealExp/Tools/subs')
      from fixed_parameters import Fixed_Parameters

      def func():
          β = fixed_parameters.β
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


EARTH_ROTATION_RATE = 7.292e-5
"""Earth rotation rate Omega in s^-1. (from MOM3 manual)"""

EARTH_RADIUS = 6.371e6
"""Mean Earth radius R in m. (from MOM3 manual)"""

ONE_DAY_IN_SECONDS = 86400.
"""1 day in seconds"""


def equatorial_beta(
    EΩ:  float = EARTH_ROTATION_RATE,
    Erad: float = EARTH_RADIUS,
) -> float:
    """Compute beta = 2 EΩ / R in s^-1 m^-1."""
    if Erad <= 0.0:
        raise ValueError("Earth radius must be positive.")
    return 2.0 * EΩ / Erad


def degrees_to_meters(
    Erad: float = EARTH_RADIUS,
) -> float:
    """Compute 1 degree in meters = 2π Erad /360."""
    π = np.pi
    if Erad <= 0.0:
        raise ValueError("Earth radius must be positive.")
    return 2. *π *Erad /360.


@dataclass(frozen=True)
class Fixed_Parameters:
    """Shared parameters."""

    β:       float = equatorial_beta()
    deg2met: float = degrees_to_meters()
    day2sec: float = ONE_DAY_IN_SECONDS

    def __post_init__(self) -> None:
        if self.β       <= 0.0:
            raise ValueError("β must be positive.")
        if self.deg2met <= 0.0:
            raise ValueError("deg2met must be positive.")
        if self.day2sec <= 0.0:
            raise ValueError("day2sec must be positive.")
