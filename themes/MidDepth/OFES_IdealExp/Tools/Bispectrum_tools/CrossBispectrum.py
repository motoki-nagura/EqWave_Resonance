import numpy as np
from numba import njit, prange
from netCDF4 import Dataset
from scipy.signal import detrend

# =========================
# FFT utility
# =========================
def compute_fft(A, dt, dz, dx):
    """
 Compute FFT

     Input:
        A           : Array[{t},{z},{x}] : Input time series
        dt, dz, dx  : Scalars            : Grid intervals

     Output:
        ω           : Array[{ω}]         : Frequency [rad {t}^-1]
        m           : Array[{m}]         : Vertical wavenumber [rad {z}^-1]
        k           : Array[{k}]         : Zonal wavenumber [rad {x}^-1]
        Ahat        : Array[{ω},{m},{k}] : Fourier transform

    """

    nω, nm, nk = A.shape
    Ahat = np.fft.fftn(A, axes=(0,1,2))

    Ahat *= dt *dz *dx
    #   (unit of Ahat) = (unit of A) /(unit of ω) /(unit of m) /(unit of k)

    ω = 2*np.pi*np.fft.fftfreq(nω, d=dt)
    m = 2*np.pi*np.fft.fftfreq(nm, d=dz)
    k = 2*np.pi*np.fft.fftfreq(nk, d=dx)

    return ω, m, k, Ahat


# =========================
# Build lookup tables
# =========================
def build_lookup_tables(ω, m, k):
    nω, nm, nk = len(ω), len(m), len(k)

    iω3_lookup = np.zeros((nω, nω), dtype=np.int64)
    im3_lookup = np.zeros((nm, nm), dtype=np.int64)
    ik3_lookup = np.zeros((nk, nk), dtype=np.int64)

    Δω, Δm, Δk = np.abs(ω[1]-ω[0]), np.abs(m[1]-m[0]), np.abs(k[1]-k[0])

    for iω1 in range(nω):
        for iω2 in range(nω):
            ω3 = ω[iω1] + ω[iω2]
            iω3_lookup[iω1, iω2] = np.argmin(np.abs(ω - ω3))

            if np.min(np.abs(ω - ω3)) > Δω*1e-2: #    no match
                iω3_lookup[iω1, iω2] = -999      #    arbitrary negative number

    for jm1 in range(nm):
        for jm2 in range(nm):
            m3 = m[jm1] + m[jm2]
            im3_lookup[jm1, jm2] = np.argmin(np.abs(m - m3))

            if np.min(np.abs(m - m3)) > Δm*1e-2: #    no match
                im3_lookup[jm1, jm2] = -999      #    arbitrary negative number

    for lk1 in range(nk):
        for lk2 in range(nk):
            k3 = k[lk1] + k[lk2]
            ik3_lookup[lk1, lk2] = np.argmin(np.abs(k - k3))

            if np.min(np.abs(k - k3)) > Δk*1e-2: #    no match
                ik3_lookup[lk1, lk2] = -999      #    arbitrary negative number

    return iω3_lookup, im3_lookup, ik3_lookup


# ==================================================
# Parallel triad kernel
# ==================================================
@njit(parallel=True)
def triad_kernel_physical_range(
        Ahat, Bhat, Chat, iω3_lookup, im3_lookup, ik3_lookup,
        bicoherence_option=None):
    """
 Compute bispectrum for the specific range

  Input:
     Ahat, Bhat, Chat   : Array[{ω},{m},{k}] : Fourier transforms of A, B, and C
     iω3_lookup         : Array[{ω},{ω}]     : Lookup tables such that ω1+ω2=ω3
     im3_lookup         : Array[{m},{m}]     : Lookup tables such that m1+m2=m3
     ik3_lookup         : Array[{k},{k}]     : Lookup tables such that k1+k2=k3
     bicoherence_option : String             : "Hinich_and_Wolinsky_2005" = Standard normalization
                                                    described in Hinich and Wolinsky (2005)
                                               "Kim_and_Powers_1979" = Kim and Powers (1979)

  Output:
     Bispec : Array[{ω1},{m1},{k1},{ω2},{m2},{k2}]
                  : Ahat[ω1,m1,k1] Bhat[ω2,m2,k2] Chat[ω1+ω2,m1+m2,k1+k2]^*,
                      where ^* is complex conjugate

     P1, P2, P3 : Array[{ω1},{m1},{k1},{ω2},{m2},{k2}]

       If bicoherence_option == "Hinich_and_Wolinsky_2005",
          P1 : Abs( Ahat[ω1,m1,k1] )**2
          P2 : Abs( Bhat[ω2,m2,k2] )**2
          P3 : Abs( Chat[ω1+ω2,m1+m2,k1+k2] )**2

       If bicoherence_option == "Kim_and_Powers_1979",
          P1 : Abs( Ahat[ω1,m1,k1] Bhat[ω2,m2,k2] )**2
          P2 : Abs( Chat[ω1+ω2,m1+m2,k1+k2] )**2
          P3 : Unity (dummy)
    """

    nω, nm, nk = Ahat.shape

    Bispec = np.zeros((nω,nm,nk, nω,nm,nk), dtype=np.complex128)
    P1     = np.zeros((nω,nm,nk, nω,nm,nk), dtype=np.float64)
    P2     = np.zeros((nω,nm,nk, nω,nm,nk), dtype=np.float64)
    P3     = np.zeros((nω,nm,nk, nω,nm,nk), dtype=np.float64)

    for iω1 in prange(nω):
        for jm1 in range(nm):
            for lk1 in range(nk):
                A = Ahat[iω1, jm1, lk1]

                for iω2 in range(nω):
                    iω3 = iω3_lookup[iω1, iω2]

                    for jm2 in range(nm):
                        jm3 = im3_lookup[jm1, jm2]

                        for lk2 in range(nk):
                            lk3 = ik3_lookup[lk1, lk2]

                            B = Bhat[iω2, jm2, lk2]
                            C = Chat[iω3, jm3, lk3]

                            if (iω3 >= 0) and (jm3 >= 0) and (lk3 >= 0):
                                Bispec[iω1,jm1,lk1, iω2,jm2,lk2] = A *B *np.conj(C)

                                if   bicoherence_option == "Hinich_and_Wolinsky_2005":
                                    P1[iω1,jm1,lk1, iω2,jm2,lk2] = np.abs( A )**2
                                    P2[iω1,jm1,lk1, iω2,jm2,lk2] = np.abs( B )**2
                                    P3[iω1,jm1,lk1, iω2,jm2,lk2] = np.abs( C )**2
                                elif bicoherence_option == "Kim_and_Powers_1979":
                                    P1[iω1,jm1,lk1, iω2,jm2,lk2] = np.abs( A*B )**2
                                    P2[iω1,jm1,lk1, iω2,jm2,lk2] = np.abs( C )**2
                                    P3[iω1,jm1,lk1, iω2,jm2,lk2] = 1. # dummy

    return Bispec, P1, P2, P3


# ==================================================
# Block-averaged bispectrum (50% overlap)
# ==================================================
def block_averaged_bispectrum(A, B, C, dt, dz, dx, block_size, ω_range, m_range, k_range,
                              window=None, save_lookup_table=False, bicoherence_option=None,
                              test_synthetic=False
                              ):
    """
 Outermost wrapper.

  Input:
     A, B, C    : Array[{t},{z},{x}] : Input time series
     dt, dz, dx : Scalars            : Grid intervals
     block_size : Scalar             : Block size
     ω_range    : Array              : Range of ω [rad {t}^-1]
                                         Length must be 2 or 4.
                                         If len=2, ω_range[0] < ω < ω_range[1].
                                         If len=4, ω_range[0] < ω < ω_range[1] or 
                                                   ω_range[2] < ω < ω_range[3].
     m_range    : Array[2]           : Range of m [rad {z}^-1]. Min and max.
     k_range    : Array[2]           : Range of k [rad {x}^-1]. Min and max.
     window     : String             : Window for time domain. Only "hanning" is available.

     save_lookup_table  : Logical    : if =True, save lookup tables in a netCDF file.

     bicoherence_option : String     :
         * "Hinich_and_Wolinsky_2005" = Standard normalization described in Hinich and
            Wolinsky (2005, Journal of Statistical Planning and Inference, 130, 405–411)
         * "Kim_and_Powers_1979" = Kim and Powers (1979, IEEE Trans. Plasma Sci. PS-7 (2),
            120–131)

     test_synthetic : Logical or Int:
         * If =Int, a synthetic data test is carried out by randomizing the phase of
             Chat following Claret et al. (2026, JPO, Vol.56, pp. 1299-1318).
             The input integer will be the number of samples.
         * If =False, calculation is performed without testing.

  Output:
     ω           : Array[{ω}]               : Frequency [rad {t}^-1]
     m           : Array[{m}]               : Vertical wavenumber [rad {z}^-1]
     k           : Array[{k}]               : Zonal wavenumber [rad {x}^-1]
     Bispec_mean : Array[{ω1},{m1},{k1},{ω2},{m2},{k2}] : Block averages of bispectrum
     Bicoherence : Array[{ω1},{m1},{k1},{ω2},{m2},{k2}] : Bicoherence
    """

    if A.shape == B.shape == C.shape:
        pass
    else:
        raise ValueError("block_averaged_bispectrum: Input A, B, and C are not in the same size.")

    nt, nz, nx = A.shape

    step = block_size // 2
    nblocks = (nt - block_size) // step + 1

    Bispec_acc = None
    P1_acc     = None
    P2_acc     = None

    if window == None:
        weight_window = np.ones(block_size)
    elif window == "hanning":
        weight_window = np.hanning(block_size)
        raise ValueError('Application of the Hanning window significantly distorts '+
                         'bicoherence. Avoid using it.')
    else:
        raise ValueError('"window" must be None or "hanning".')

    if ( (not bicoherence_option == "Hinich_and_Wolinsky_2005") and
         (not bicoherence_option == "Kim_and_Powers_1979") ):
        raise ValueError('bicoherence_option must be "Hinich_and_Wolinsky_2005" or '+
                         '"Kim_and_Powers_1979"'
                         )

    if test_synthetic is False:
        Do_Test_Synthetic = False
        n_Test_Synthetic = 0
    elif isinstance(test_synthetic, int) and not isinstance(test_synthetic, bool):
        Do_Test_Synthetic = True
        n_Test_Synthetic = test_synthetic
    else:
        raise TypeError("test_synthetic must be False or an integer")

    #    <--------- Block loop -------->
    for ib in range(nblocks):

        #    Index range for a block

        start = ib * step
        end = start + block_size
        sl = slice(start, end)

        #    Apply window

        A_input = detrend(A[sl,:,:], type='constant', axis=0)  #  subtract the mean
        B_input = detrend(B[sl,:,:], type='constant', axis=0)
        C_input = detrend(C[sl,:,:], type='constant', axis=0)

        A_input = A_input *weight_window[:,np.newaxis,np.newaxis]
        B_input = B_input *weight_window[:,np.newaxis,np.newaxis]
        C_input = C_input *weight_window[:,np.newaxis,np.newaxis]

        #    Compute Fourier transform

        ω, m, k, Ahat = compute_fft(A_input, dt, dz, dx)
        _, _, _, Bhat = compute_fft(B_input, dt, dz, dx)
        _, _, _, Chat = compute_fft(C_input, dt, dz, dx)

        #  Amplitude adjustment for the Hanning window

        if window == "hanning":
            Ahat = Ahat *np.sqrt(8./3.)
            Bhat = Bhat *np.sqrt(8./3.)
            Chat = Chat *np.sqrt(8./3.)

        #    In the first computation (ib=0), prepare mask and the lookup tables

        if ib == 0:
            if len(ω_range) == 2:
                ω_mask = (ω >= ω_range[0]) & (ω <= ω_range[1])
            elif len(ω_range) == 4:
                ω_mask = ( ( (ω_range[0] < ω) & (ω < ω_range[1]) ) | 
                           ( (ω_range[2] < ω) & (ω < ω_range[3]) ) )
                #print('ω_range   = ',ω_range)
                #print('ω[ω_mask] = ',ω[ω_mask])
            else:
                raise ValueError(f"len(ω_range) is must be 2 or 4. Now it is {len(ω_range):d}. Stop.")

            m_min = m[m <= m_range[0]].max()
            m_max = m[m >= m_range[1]].min()
            m_mask = (m >= m_min) & (m <= m_max)   #   End points included

            k_min = k[k <= k_range[0]].max()
            k_max = k[k >= k_range[1]].min()
            k_mask = (k >= k_min) & (k <= k_max)   #   End points included

            #m_mask = (m >= m_range[0]) & (m <= m_range[1])
            #k_mask = (k >= k_range[0]) & (k <= k_range[1])

            ω_idx = np.where(ω_mask)[0]
            m_idx = np.where(m_mask)[0]
            k_idx = np.where(k_mask)[0]

            iω3_lookup, im3_lookup, ik3_lookup = build_lookup_tables(ω, m, k)

            iω3_lookup_r = iω3_lookup[ω_idx,:][:,ω_idx]
            im3_lookup_r = im3_lookup[m_idx,:][:,m_idx]
            ik3_lookup_r = ik3_lookup[k_idx,:][:,k_idx]

            if save_lookup_table:
                save_lookup_tables(ω,m,k,
                                   ω[ω_idx],m[m_idx],k[k_idx],
                                   iω3_lookup_r,im3_lookup_r,ik3_lookup_r)

            #print(f'{np.min(ω):.3e} < ω < {np.max(ω):.3e}')
            #print(f'{np.min(m):.3e} < m < {np.max(m):.3e}')
            #print(f'{np.min(k):.3e} < k < {np.max(k):.3e}')
            #print(f'  ω_range = {ω_range[0]:.3e} to {ω_range[1]:.3e},  {np.sum(ω_mask):d} grids')
            #print(f'  m_range = {m_range[0]:.3e} to {m_range[1]:.3e},  {np.sum(m_mask):d} grids')
            #print(f'  k_range = {k_range[0]:.3e} to {k_range[1]:.3e},  {np.sum(k_mask):d} grids')
            #for ll in range(len(ω)):
            #    print(f'll={ll:d}, ω={ω[ll]:.3e}, ω_mask={ω_mask[ll]:d}')
            #raise ValueError('This is a test and the program stops here')

        #    Limit to the prescribed range

        Ahat_r = Ahat[ω_idx,:,:][:,m_idx,:][:,:,k_idx]  #  range is limited
        Bhat_r = Bhat[ω_idx,:,:][:,m_idx,:][:,:,k_idx]  #  range is limited
        Chat_r = Chat                                   #  range is not limited

        ω = ω[ω_idx]
        m = m[m_idx]
        k = k[k_idx]

        #    Compute the bispectrum

        Bispec, P1, P2, P3 = triad_kernel_physical_range(
            Ahat_r, Bhat_r, Chat_r,
            iω3_lookup_r, im3_lookup_r, ik3_lookup_r,
            bicoherence_option,
        )

        #test_check_test(Bispec,ω,m,k)  #  CHECK, if necessary

        #    Optional synthetic data test (Bispec overwritten)

        if Do_Test_Synthetic:
            Bispec, Bispec_Re_sq = \
            triad_kernel_physical_range_For_Test_Synthetic(
                    Ahat_r, Bhat_r, Chat_r,
                    iω3_lookup_r, im3_lookup_r, ik3_lookup_r,
                    bicoherence_option,
                    n_Test_Synthetic,
                    )

        #    Accumulate results

        if Bispec_acc is None:
            Bispec_acc = Bispec
            P1_acc     = P1
            P2_acc     = P2
            P3_acc     = P3
            if Do_Test_Synthetic:
                Bispec_Re_sq_acc = Bispec_Re_sq
        else:
            Bispec_acc += Bispec
            P1_acc     += P1
            P2_acc     += P2
            P3_acc     += P3
            if Do_Test_Synthetic:
                Bispec_Re_sq_acc += Bispec_Re_sq

    #    Compute block average

    Bispec_mean = Bispec_acc / nblocks
    P1_mean     = P1_acc     / nblocks
    P2_mean     = P2_acc     / nblocks
    P3_mean     = P3_acc     / nblocks

    if Do_Test_Synthetic:
        Bispec_Re_sq_mean = Bispec_Re_sq_acc / nblocks

    #test_check_test(Bispec_mean,ω,m,k)  #  CHECK

    #    Compute bicoherence (Eq.(2.3) in Hinich and Wolinsky 2005;
    #                         Eq.(2) in Elgar and Guza 1988)

    Bicoherence = \
        np.abs(Bispec_mean) / (np.sqrt(P1_mean * P2_mean * P3_mean) + 1e-12)

    #  Sort in terms of ω, k, and m

    iω_sort = np.argsort(ω)
    im_sort = np.argsort(m)
    ik_sort = np.argsort(k)

    ω = ω[iω_sort]
    m = m[im_sort]
    k = k[ik_sort]

    Bispec_mean = Bispec_mean[iω_sort,:,:, :,:,:][:,im_sort,:, :,:,:][:,:,ik_sort, :,:,:]# For ω1,k1,m1
    Bicoherence = Bicoherence[iω_sort,:,:, :,:,:][:,im_sort,:, :,:,:][:,:,ik_sort, :,:,:]

    Bispec_mean = Bispec_mean[:,:,:, iω_sort,:,:][:,:,:, :,im_sort,:][:,:,:, :,:,ik_sort]# For ω2,k2,m2
    Bicoherence = Bicoherence[:,:,:, iω_sort,:,:][:,:,:, :,im_sort,:][:,:,:, :,:,ik_sort]

    if Do_Test_Synthetic:
        Bispec_Re_sq_mean = \
        Bispec_Re_sq_mean[iω_sort,:,:, :,:,:][:,im_sort,:, :,:,:][:,:,ik_sort, :,:,:]# For ω1,k1,m1

        Bispec_Re_sq_mean = \
        Bispec_Re_sq_mean[:,:,:, iω_sort,:,:][:,:,:, :,im_sort,:][:,:,:, :,:,ik_sort]# For ω2,k2,m2

    #  Return

    if Do_Test_Synthetic:
        return ω, m, k, Bispec_mean, Bicoherence, Bispec_Re_sq_mean
    else:
        return ω, m, k, Bispec_mean, Bicoherence


# =========================
# netCDF output
# =========================
def save_to_netcdf(filename, ω, m, k, Bispec, Bicoherence):
    nc = Dataset(filename, "w")

    nc.createDimension("ω1", len(ω))
    nc.createDimension("m1", len(m))
    nc.createDimension("k1", len(k))

    nc.createDimension("ω2", len(ω))
    nc.createDimension("m2", len(m))
    nc.createDimension("k2", len(k))

    nc.createVariable("ω", "f8", ("ω1",))[:] = ω
    nc.createVariable("m", "f8", ("m1",))[:] = m
    nc.createVariable("k", "f8", ("k1",))[:] = k

    nc.createVariable("bispectrum_real", "f8", ("ω1","m1","k1","ω2","m2","k2"))[:] = Bispec.real
    nc.createVariable("bispectrum_imag", "f8", ("ω1","m1","k1","ω2","m2","k2"))[:] = Bispec.imag
    nc.createVariable("bicoherence",     "f8", ("ω1","m1","k1","ω2","m2","k2"))[:] = Bicoherence

    nc.close()


def save_lookup_tables(ω,m,k, ω_r,m_r,k_r, iω3_lookup_r,im3_lookup_r,ik3_lookup_r,
                       filename='lookup_tables.nc'):
    nc = Dataset(filename, "w")

    nc.createDimension("ω", len(ω))
    nc.createDimension("m", len(m))
    nc.createDimension("k", len(k))

    nc.createDimension("ω_r", len(ω_r))
    nc.createDimension("m_r", len(m_r))
    nc.createDimension("k_r", len(k_r))

    nc.createVariable("ω", "f8", ("ω",))[:] = ω            #   All grids of ω
    nc.createVariable("m", "f8", ("m",))[:] = m
    nc.createVariable("k", "f8", ("k",))[:] = k

    nc.createVariable("ω_r", "f8", ("ω_r",))[:] = ω_r      #   ω in the specified range
    nc.createVariable("m_r", "f8", ("m_r",))[:] = m_r
    nc.createVariable("k_r", "f8", ("k_r",))[:] = k_r

    nc.createVariable("iω3_lookup", "int", ("ω_r","ω_r",))[:] = iω3_lookup_r   #    lookup in the range
    nc.createVariable("im3_lookup", "int", ("m_r","m_r",))[:] = im3_lookup_r
    nc.createVariable("ik3_lookup", "int", ("k_r","k_r",))[:] = ik3_lookup_r

    nc.close()


# =========================
# text output
# =========================
def test_check_test(Bispec,ω,m,k):
    #ttl, ω1,ω2,k1,m1,k2,m2 = (
    #        'R',  -2.*np.pi/(180.*86400.), -2.*np.pi/( 90.*86400.), \
    #        -2.5e-6, 4.0e-3, -6.5e-6, -3.5e-3 ) # resonant
    ttl, ω1,ω2,k1,m1,k2,m2 = (
           'OR', -2.*np.pi/(180.*86400.), -2.*np.pi/( 90.*86400.), \
            -2.5e-6, 4.0e-3, -5.0e-6, 8.0e-3 ) # off resonant
    iω1 = np.argmin(np.abs( ω - ω1 ))
    iω2 = np.argmin(np.abs( ω - ω2 ))
    ik1 = np.argmin(np.abs( k - k1 ))
    im1 = np.argmin(np.abs( m - m1 ))
    ik2 = np.argmin(np.abs( k - k2 ))
    im2 = np.argmin(np.abs( m - m2 ))
    print(ttl+"  "+
          f"Bispec = {Bispec[iω1,im1,ik1,iω2,im2,ik2]:.4e} at "+
          f"2π/ω1={2.*np.pi/ω[iω1]/86400:.1f} dy, k1={k[ik1]:.3e}, m1={m[im1]:.3e}, "+
          f"2π/ω2={2.*np.pi/ω[iω2]/86400:.1f} dy, k2={k[ik2]:.3e}, m2={m[im2]:.3e}"
         )
    return


# ==================================================
#  For synthetic data test
# ==================================================
def Randomize_Phase(var):
    phase = np.random.uniform(0, 2*np.pi, size=var.shape)
    var_random_phase = np.abs(var) * np.exp(1j * phase)
    return var_random_phase


def triad_kernel_physical_range_For_Test_Synthetic(
        Ahat_r, Bhat_r, Chat_r,
        iω3_lookup_r, im3_lookup_r, ik3_lookup_r,
        bicoherence_option,
        n_Test_Synthetic,
        ):

    Bispec_test = None

    for iloop_test in range(0,n_Test_Synthetic):
        Chat_r_randomized = Randomize_Phase(Chat_r)   #  Synthetic

        Bispec, _, _, _ = triad_kernel_physical_range(
            Ahat_r, Bhat_r, Chat_r_randomized,
            iω3_lookup_r, im3_lookup_r, ik3_lookup_r,
            bicoherence_option,
        )

        if Bispec_test is None:
            Bispec_test       = Bispec
            Bispec_Re_sq_test = Bispec.real **2
        else:
            Bispec_test       += Bispec
            Bispec_Re_sq_test += Bispec.real **2

    Bispec       = Bispec_test       / n_Test_Synthetic #  Results is from randomized Chat
    Bispec_Re_sq = Bispec_Re_sq_test / n_Test_Synthetic

    return Bispec, Bispec_Re_sq

# =========================
# Example usage
# =========================
if __name__ == "__main__":
#{{{

    # Grids
    nt, nz, nx = 64, 16, 16
    dt, dz, dx = 0.1, 0.25, 0.3
    block_size = 32

    # Frequency parameters
    π = np.pi
    dω = 2.*π/(block_size*dt)   #  Frequency grid interval
    dm = 2.*π/(nz*dz)
    dk = 2.*π/(nx*dx)

    ω0, m0, k0 = dω*3, -dm*3, dk*5
    print(f'ω0 = {ω0:.2e}; m0 = {m0:.2e}; k0 = {k0:.2e}')

    # Spectral ranges to compute
    ω_range = (-dω*5, dω*5)
    m_range = (-dm*5, dm*5)
    k_range = (-dk*8, dk*8)

    # Generate grids
    t = np.arange(nt)*dt
    z = np.arange(nz)*dz
    x = np.arange(nx)*dx
    T, Z, X = np.meshgrid(t, z, x, indexing="ij")

    # Generate input time series
    Ain = np.sin(T*ω0) * np.cos(Z*m0)
    Bin = np.cos(T*ω0) * np.sin(X*k0)
    Cin = 2. *Ain *Bin

    # Call the subroutine
    ω, m, k, Bispec, Bicoherence = block_averaged_bispectrum(
        Ain, Bin, Cin,
        dt, dz, dx,
        block_size,
        ω_range, m_range, k_range,
        bicoherence_option="Hinich_and_Wolinsky_2005",
        test_synthetic=False,
    )
#        window='hanning'

    # Write out
    save_to_netcdf("out.nc", ω, m, k, Bispec, Bicoherence)

    print('np.max(Bicoherence) = ',np.max(Bicoherence))
    print('ω = ',ω)
    print('m = ',m)
    print('k = ',k)

    #------------------------------#
    #    Computation for figure
    #------------------------------#

    #   Indices

    iω1 = np.argmin(np.abs(ω - ω0))
    im1 = np.argmin(np.abs(m - m0))
    ik1 = np.argmin(np.abs(k - 0.))

    iω2 = iω1
    im2 = np.argmin(np.abs(m - 0.))
    ik2 = np.argmin(np.abs(k - k0))

    iω3 = np.argmin(np.abs(ω - 2.*ω0))

    #   Power spectrum

    nt, nz, nx = Ain.shape
    step = block_size // 2
    nblocks = (nt - block_size) // step + 1
    #
    Ahat = None
    #
    for ib in range(nblocks):
        start = ib * step
        end = start + block_size
        sl = slice(start, end)

        ω_ps, m_ps, k_ps, Ahat_1 = compute_fft(Ain[sl], dt, dz, dx)
        _,    _,    _,    Bhat_1 = compute_fft(Bin[sl], dt, dz, dx)
        _,    _,    _,    Chat_1 = compute_fft(Cin[sl], dt, dz, dx)

        if Ahat is None:
            Ahat  = Ahat_1
            Bhat  = Bhat_1
            Chat  = Chat_1
        else:
            Ahat += Ahat_1
            Bhat += Bhat_1
            Chat += Chat_1

    Ahat = Ahat / nblocks
    Bhat = Bhat / nblocks
    Chat = Chat / nblocks

    iω_sort = np.argsort(ω_ps)
    im_sort = np.argsort(m_ps)
    ik_sort = np.argsort(k_ps)

    ω_ps = ω_ps[iω_sort]
    m_ps = m_ps[im_sort]
    k_ps = k_ps[ik_sort]

    Ahat = Ahat[iω_sort,:,:][:,im_sort,:][:,:,ik_sort]
    Bhat = Bhat[iω_sort,:,:][:,im_sort,:][:,:,ik_sort]
    Chat = Chat[iω_sort,:,:][:,im_sort,:][:,:,ik_sort]

    #   Indices

    iω1_ps = np.argmin(np.abs(ω_ps - ω0))
    iω2_ps = iω1_ps
    iω3_ps = np.argmin(np.abs(ω_ps - 2.*ω0))

    #---------------#
    #    Figures
    #---------------#
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    with PdfPages('fig.pdf') as pdf:

        #---------------------
        #    Time series

        fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(8,8))
        axs = axs.flatten()     #   row major
        iplot = 0
        #
        i0,k0 = int(nx*0.3), int(nz*0.5)
        #
        tmp = f'[:,k0={k0:d},i0={i0:d}]'
        for var, ttl in zip([Ain[:,k0,i0],Bin[:,k0,i0],Cin[:,k0,i0]],
                            [r"$A$"+tmp, r"$B$"+tmp, r"$C$"+tmp]):
            axs[iplot].set_title(ttl)
            axs[iplot].plot(t,var,'k-')
            axs[iplot].set_xlabel("Time")
            iplot += 1
        #
        plt.tight_layout()
        pdf.savefig()
        plt.close(fig)

        #---------------------
        #    Spectral Power

        fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(8,8))
        axs = axs.flatten()     #   row major
        iplot = 0
        #
        #
        for var, ttl in zip([Ahat[iω1_ps,:,:], Bhat[iω2_ps,:,:], Chat[iω3_ps,:,:]],
                            [r"$A$ at $\omega=\omega_0$", r"$B$ at $\omega=\omega_0$",
                             r"$C$ at $\omega=2\omega_0$"]):
            axs[iplot].set_title("PSD of "+ttl)
            psd = np.abs(var)**2 /(nt*dt) /(nz*dz) /(nx*dx)
            pcm = axs[iplot].pcolormesh(k_ps,m_ps, np.abs(var)**2, shading="auto")
            fig.colorbar(pcm, ax=axs[iplot], label="Spectral Power")
            axs[iplot].set_xlabel(r"$k$")
            axs[iplot].set_ylabel(r"$m$")
            axs[iplot].set_xlim(np.array(k_range))
            axs[iplot].set_ylim(np.array(m_range))
            iplot += 1
        #
        fig.delaxes(axs[iplot])
        #
        plt.tight_layout()
        pdf.savefig()
        plt.close(fig)

        #------------------
        #    Bicoherence

        fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(8,8))
        axs = axs.flatten()     #   row major
        iplot = 0
        #
        axs[iplot].set_title(r"$\omega_1=\omega_0; \omega_2=\omega_0; m_2=0; k_2=k_0$")
        pcm = axs[iplot].pcolormesh(k,m, Bicoherence[iω1,:,:, iω2,im2,ik2], shading="auto")
        fig.colorbar(pcm, ax=axs[iplot], label="Bicoherence")
        axs[iplot].set_xlabel(r"$k_1$")
        axs[iplot].set_ylabel(r"$m_1$")
        axs[iplot].axvline(x=k[ik1], c='w', ls=':')
        axs[iplot].axhline(y=m[im1], c='w', ls=':')
        axs[iplot].set_xlim(np.array(k_range))
        axs[iplot].set_ylim(np.array(m_range))
        iplot += 1
        #
        axs[iplot].set_title(r"$\omega_1=\omega_0; m_1=m_0; k_1=0; \omega_2=\omega_0$")
        pcm = axs[iplot].pcolormesh(k,m, Bicoherence[iω1,im1,ik1, iω2,:,:], shading="auto")
        fig.colorbar(pcm, ax=axs[iplot], label="Bicoherence")
        axs[iplot].set_xlabel(r"$k_2$")
        axs[iplot].set_ylabel(r"$m_2$")
        axs[iplot].axvline(x=k[ik2], c='w', ls=':')
        axs[iplot].axhline(y=m[im2], c='w', ls=':')
        axs[iplot].set_xlim(np.array(k_range))
        axs[iplot].set_ylim(np.array(m_range))
        iplot += 1
        #
        axs[iplot].set_title(r"$m_1=m_0; k_1=0; m_2=0; k_2=k_0$")
        pcm = axs[iplot].pcolormesh(ω,ω, Bicoherence[:,im1,ik1, :,im2,ik2], shading="auto")
        fig.colorbar(pcm, ax=axs[iplot], label="Bicoherence")
        axs[iplot].set_xlabel(r"$\omega_2$")
        axs[iplot].set_ylabel(r"$\omega_1$")
        axs[iplot].axvline(x=ω[iω2], c='w', ls=':')
        axs[iplot].axhline(y=ω[iω1], c='w', ls=':')
        axs[iplot].set_xlim(np.array(ω_range))
        axs[iplot].set_ylim(np.array(ω_range))
        iplot += 1
        #
        fig.delaxes(axs[iplot])
        #
        plt.tight_layout()
        pdf.savefig()
        plt.close(fig)

    #--------------------------------------#
    #               Text Check
    #--------------------------------------#

    ## Check if the resonance condition is satisfied
    #
    #iω3_lookup, im3_lookup, ik3_lookup = build_lookup_tables(ω, m, k)
    #nω = len(ω)
    #ω_mask = (ω >= ω_range[0]) & (ω <= ω_range[1])
    #for iω1 in prange(nω):
    #    for iω2 in range(nω):
    #        iω3 = iω3_lookup[iω1, iω2]
    #        if not ω_mask[iω3]:
    #            continue
    #        tmp_max = np.max([np.abs(ω[iω1]), np.abs(ω[iω2]), np.abs(ω[iω3])])
    #        print(f'iω1={iω1:2d}, iω2={iω2:2d}, iω3={iω3:2d}; '+\
    #              f'ω[iω1]={ω[iω1]:.3e}, ω[iω2]={ω[iω2]:.3e}, ω[iω3]={ω[iω3]:.3e}; '+\
    #              f'(ω[iω1]+ω[iω2]-ω[iω3])/max(abs)={(ω[iω1]+ω[iω2]-ω[iω3])/tmp_max:.3e}'
    #              )
#}}}
