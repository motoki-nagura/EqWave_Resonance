import numpy as np
import pyfftw

def fft_spectrum_3d ( var, xcrd, ycrd, tcrd, verbose_check=False ):
    """
----------------------------------------------------------------------------
     3D Power spectrum by FFT
----------------------------------------------------------------------------

----------------------------------------------------------------------------
 Input

     var[{time},{lat},{lon}] :  Input data
    xcrd[{lon}]              :  X coordinate
    ycrd[{lat}]              :  Y coordinate
    tcrd[{time}]             :  T coordinate

 Output

    psd[{freq_t},{freq_y},{freq_x}] :  Power spectrum density;  Unit is (var**2) *xcrd *ycrd *tcrd
    freq_x[{freq_x}]                :  Frequency in X space;  Unit is 1/xcrd
    freq_y[{freq_y}]                :  Frequency in Y space;  Unit is 1/ycrd
    freq_t[{freq_t}]                :  Frequency in T space;  Unit is 1/tcrd

----------------------------------------------------------------------------"""

    dx,   dy,   dt   = xcrd[ 1]-xcrd[0], ycrd[ 1]-ycrd[0], tcrd[ 1]-tcrd[0]   #  grid intervals
    xlen, ylen, tlen = xcrd[-1]-xcrd[0], ycrd[-1]-ycrd[0], tcrd[-1]-tcrd[0]   #  length of coordinates

    Yk = pyfftw.interfaces.numpy_fft.fftn(var, axes=(0,1,2), threads=8) *dx *dy *dt  #  multi-threaded (sped up by 20%)

    #Yk = np.fft.fftn(var, axes=(0,1,2)) *dx *dy *dt

    #     dx and dt are multiplied to Yk above. This is due to the discrete sampling
    #       (see Eq. (5.24b) and the discussion below the equation in Thomson and Emery 2014)
    #       and the formulation of numpy fft2 (see
    #       https://numpy.org/doc/stable/reference/routines.fft.html#higher-dimensions)

    freq_x = np.fft.fftfreq(len(xcrd), d=dx)
    freq_y = np.fft.fftfreq(len(ycrd), d=dy)
    freq_t = np.fft.fftfreq(len(tcrd), d=dt)

    psd = np.abs(Yk)**2 /xlen /ylen /tlen                 #  Eq.(5.33a) in Thomson and Emery (2014)

    ix_1neg = np.min(np.argwhere( freq_x < 0. )[:,0])     #  first negative frequency in freq_x
    freq_x = np.roll(freq_x, -ix_1neg)                    #  shift the array
    psd    = np.roll(psd,    -ix_1neg, axis=2)
    del ix_1neg

    iy_1neg = np.min(np.argwhere( freq_y < 0. )[:,0])     #  first negative frequency in freq_y
    freq_y = np.roll(freq_y, -iy_1neg)                    #  shift the array
    psd    = np.roll(psd,    -iy_1neg, axis=1)
    del iy_1neg

    lt_1neg = np.min(np.argwhere( freq_t < 0. )[:,0])     #  first negative frequency in freq_t
    freq_t = np.roll(freq_t, -lt_1neg)                    #  shift the array
    psd    = np.roll(psd,    -lt_1neg, axis=0)
    del lt_1neg

    freq_t = freq_t*(-1.)        #  The sign of frequency in time is inverted. Eastward propagating
                                 #    signals are in the quadrant with positive k and positive omega.

    #    Check
    if verbose_check:
      lhs = dx *dy *dt *np.sum( np.abs(var)**2 ) #  Perseval's theorem  (below Eq. 5.29 in
      rhs = np.sum(psd)                          #   Thomson and Emery 2014)

      psd_df = np.sum(psd)/(xlen*ylen*tlen)      #  must be equal to signal variance.
                                                 #   Eq.(5.35) in Thomson and Emery (2014)

      print(f"Check Perseval's theorem: lhs, rhs    = {lhs:10.4e}, {rhs:10.4e}")
      print(f'Check variance: psd*df, np.var(input) = {psd_df:10.4e}, {np. var(var):10.4e}')

    return psd, freq_x, freq_y, freq_t
