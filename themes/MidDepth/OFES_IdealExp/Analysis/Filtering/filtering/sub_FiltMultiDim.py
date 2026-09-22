def ApplyFilter(var,W):
    """
       Apply filtering in ω-kx-kz space

           Input
                var    :  Array[{time},{lev},{lon}]    Input time series
                W      :  Array[{ω},{kz},{kx}]         Filter function in frequency space

           Output
                var_filt:  Array[{time},{lev},{lon}]   Filtered time series

    """
    import numpy as np
    import pyfftw

    # ==========================================================
    # 1. FFT
    # ==========================================================
    spec = pyfftw.interfaces.numpy_fft.fftn(var, axes=(0, 1, 2), threads=8)
    #spec = np.fft.fftn(var, axes=(0, 1, 2))

    # ==========================================================
    # 2. Spectral filter
    # ==========================================================
    spec_filt = spec * W

    # ==========================================================
    # 3. Inverse FFT
    # ==========================================================
    var_filt = pyfftw.interfaces.numpy_fft.ifftn(spec_filt, axes=(0, 1, 2), threads=8).real
    #var_filt = np.fft.ifftn(spec_filt, axes=(0, 1, 2)).real

    return var_filt

#-------------------------------------------------------------------------

def MakeFilter (x,z,t, ω_cent,ω_hwid,ω_tap, kx_cent,kx_hwid,kx_tap, kz_cent,kz_hwid,kz_tap):
    """
   Generate a filter function

     Input
       x       : Array[{lon}]  :  lon coordinate
       z       : Array[{lev}]  :  lev coordinate
       t       : Array[{time}] :  time coordinate

       ω_cent  : Scalar        :  Center of the filter in ω [rad time^-1]
       ω_hwid  : Scalar        :  Half width of the filter in ω [rad time^-1]
       ω_tap   : Scalar        :  Width of cosine taper on each side [rad time^-1]
       kx_cent, kx_hwid, kx_tap : Scalars :  Same as ω_cent, ω_hwid, and ω_tap, but for kx [rad lon^-1]
       kz_cent, kz_hwid, kz_tap : Scalars :  Same as ω_cent, ω_hwid, and ω_tap, but for kz [rad lev^-1]

    """
    import numpy as np

    def cosine_band(x, xcent, width, taper):
        """
        Smooth band-pass window centered at xcent.

        Parameters
        ----------
        x      : array
        xcent  : center
        width  : half-width of flat passband
        taper  : width of cosine taper on each side
        """
        w = np.zeros_like(x, dtype=float)
        d = np.abs(x - xcent)

        # flat passband
        core = d <= width
        w[core] = 1.0

        # cosine taper
        edge = (d > width) & (d <= width + taper)
        w[edge] = 0.5 * (
            1 + np.cos(np.pi * (d[edge] - width) / taper)
        )

        return w

    # ==========================================================
    # 1. Coordinates
    # ==========================================================
    nx,nz,nt = len(x),    len(z),    len(t)
    dx,dz,dt = x[1]-x[0], z[1]-z[0], t[1]-t[0]

    omega = np.fft.fftfreq(nt, d=dt) * 2 * np.pi
    kz    = np.fft.fftfreq(nz, d=dz) * 2 * np.pi
    kx    = np.fft.fftfreq(nx, d=dx) * 2 * np.pi

    omega *= -1.   #  Sign of omega flipped, see fft_spectrum_3d in subs.py 

    Ω, KZ, KX = np.meshgrid(omega, kz, kx, indexing="ij")

    # ==========================================================
    # 2. Filter
    # ==========================================================

    # Frequency window
    W_omega = cosine_band(Ω,   ω_cent,  ω_hwid,  ω_tap)

    # Zonal wavenumber window
    W_kx    = cosine_band(KX, kx_cent, kx_hwid, kx_tap)

    # Vertical wavenumber window
    W_kz    = cosine_band(KZ, kz_cent, kz_hwid, kz_tap)

    # Generate a 3D tapered filter
    W = W_omega * W_kx * W_kz

    return W, omega,kz,kx

