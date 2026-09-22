from dataclasses import dataclass
import numpy
import pycwt as wavelet
from pycwt.helpers import find


@dataclass
class WaveletResult:
    # Preserve full output as attributes
    wave: numpy.ndarray
    scales: numpy.ndarray
    freqs: numpy.ndarray
    coi: numpy.ndarray
    fft: numpy.ndarray
    fftfreqs: numpy.ndarray
    iwave: numpy.ndarray
    power: numpy.ndarray
    fft_power: numpy.ndarray
    fft_theor: numpy.ndarray
    sig95: numpy.ndarray
    glbl_power: numpy.ndarray
    glbl_signif: numpy.ndarray
    period: numpy.ndarray
    scale_avg: numpy.ndarray
    scale_avg_signif: numpy.ndarray
    mother: object


def run_wavelet_analysis(dat, dt, dj=1/12, s0_factor=2, J_factor=7, Tmin=None, Tmax=None):
    """
    Performs the full PyCWT wavelet analysis exactly as in the original script.
    All comment lines are preserved.

       Input :   dat[:]  :  Time series
                 dt      :  Sampling interval
                 Tmin[:] :  Lower limit of the period for averaging
                 Tmax[:] :  Higher limit of the period for averaging

       Output : many
    """

    #         wavelet function

    mother = wavelet.Morlet(6)
    s0 = s0_factor * dt             # Starting scale
    J = J_factor / dj               # Seven powers of two with dj sub-octaves
    alpha, _, _ = ar1(dat)          # Lag-1 autocorrelation for red noise. A straightforward estimate
    #print('alpha = ',alpha)         #    is used when the sophisticated method does not work.
    #alpha, _, _ = wavelet.ar1(dat)  # Lag-1 autocorrelation for red noise. Original.

    #        wavelet transform

    # The following routines perform the wavelet transform and inverse wavelet
    # transform using the parameters defined above. Since we have normalized our
    # input time-series, we multiply the inverse transform by the standard
    # deviation.
    wave_, scales, freqs, coi, fft, fftfreqs = wavelet.cwt(dat, dt, dj, s0, J, mother)
    iwave = wavelet.icwt(wave_, scales, dt, dj, mother)   #  inverse

    # We calculate the normalized wavelet and Fourier power spectra, as well as
    # the Fourier equivalent periods for each wavelet scale.
    power = (numpy.abs(wave_)) ** 2
    fft_power = numpy.abs(fft) ** 2

    # We could stop at this point and plot our results. However we are also
    # interested in the power spectra significance test. The power is significant
    # where the ratio ``power / sig95 > 1``.
    N = len(dat)
    signif, fft_theor = wavelet.significance(
        1.0, dt, scales, 0, alpha,
        significance_level=0.95,
        wavelet=mother
    )
    sig95 = numpy.ones([1, N]) * signif[:, None]
    sig95 = power / sig95

    # Then, we calculate the global wavelet spectrum and determine its
    # significance level.
    glbl_power = power.mean(axis=1)
    dof = N - scales  # Correction for padding at edges
    var = numpy.var(dat)
    glbl_signif, _ = wavelet.significance(
        var, dt, scales, 1, alpha,
        significance_level=0.95, dof=dof,
        wavelet=mother
    )

    # We also calculate the scale average over a period band, and its
    # significance level.
    period = 1 / freqs
    Cdelta = mother.cdelta
    scale_avg        = numpy.zeros((len(Tmin), N))
    scale_avg_signif = numpy.zeros((len(Tmin)))

    for p in range(len(Tmin)):
        sel = find((period >= Tmin[p]) & (period < Tmax[p]))
        scale_avg_1 = (scales * numpy.ones((N, 1))).transpose()
        scale_avg_1 = power / scale_avg_1  # As in Torrence and Compo (1998) equation 24
        scale_avg_1 = var * dj * dt / Cdelta * scale_avg_1[sel, :].sum(axis=0)
        scale_avg_signif_1, _ = wavelet.significance(
            var, dt, scales, sigma_test=2, alpha=alpha,
            significance_level=0.95,
            dof=[scales[sel[0]], scales[sel[-1]]],
            wavelet=mother
        )
        scale_avg[       p,:] = scale_avg_1
        scale_avg_signif[p]   = scale_avg_signif_1

    # Return results in a convenient dataclass
    return WaveletResult(
        wave=wave_, scales=scales, freqs=freqs, coi=coi,
        fft=fft, fftfreqs=fftfreqs, iwave=iwave,
        power=power, fft_power=fft_power, fft_theor=fft_theor, sig95=sig95,
        glbl_power=glbl_power, glbl_signif=glbl_signif, period=period,
        scale_avg=scale_avg, scale_avg_signif=scale_avg_signif,
        mother=mother
    )


def ar1(x):
    r"""Lag-1 autocorrelation coefficient.

    In an Allen and Smith AR(1) model,

    $$  x(t) - <x> = \gamma(x(t-1) - <x>) + \alpha z(t) ,$$

    where $<x>$ is the process mean, $\gamma$ and $\alpha$ are process
    parameters and $z(t)$ is a Gaussian unit-variance white noise.

    Parameters
    ----------
    x : numpy.ndarray, list
        Univariate time series

    Returns
    -------
    g : float
        Estimate of the lag-one autocorrelation.
    a : float
        Estimate of the noise variance [var(x) ~= a**2/(1-g**2)]
    mu2 : float
        Estimated square on the mean of a finite segment of AR(1)
        noise, mormalized by the process variance.

    References
    ----------
    [1] Allen, M. R. and Smith, L. A. Monte Carlo SSA: detecting
        irregular oscillations in the presence of colored noise.
        *Journal of Climate*, **1996**, 9(12), 3373-3404.
        <http://dx.doi.org/10.1175/1520-0442(1996)009<3373:MCSDIO>2.0.CO;2>
    [2] http://www.madsci.org/posts/archives/may97/864012045.Eg.r.html


    ** Copied from https://pycwt.readthedocs.io/en/development/reference/#pycwt.ar1 **

    """
    x = numpy.asarray(x)
    N = x.size
    xm = x.mean()
    x = x - xm

    # Estimates the lag zero and one covariance
    c0 = x.transpose().dot(x) / N
    c1 = x[0 : N - 1].transpose().dot(x[1:N]) / (N - 1)

    # According to A. Grinsteds' substitutions
    B = -c1 * N - c0 * N**2 - 2 * c0 + 2 * c1 - c1 * N**2 + c0 * N
    A = c0 * N**2
    C = N * (c0 + c1 * N - c1)
    D = B**2 - 4 * A * C

    if D > 0:
        g = (-B - D**0.5) / (2 * A)
    else:
        #print(                                #  --> Added by nagura (2025/11/25)
        #    "Cannot place an upperbound on the unbiased AR(1). "
        #    "Series is too short or trend is too large.")
        import numpy as np
        g = np.corrcoef(x[:-1], x[1:])[0,1]   #  Substitute a simple solution
        #                                     #  <-- Added by nagura (2025/11/25)
        #raise Warning(                       #  --> Commented out by nagura (2025/11/25)
        #    "Cannot place an upperbound on the unbiased AR(1). "
        #    "Series is too short or trend is too large."
        #)                                    #  <-- Commented out by nagura (2025/11/25)

    # According to Allen & Smith (1996), footnote 4
    mu2 = -1 / N + (2 / N**2) * (
        (N - g**N) / (1 - g) - g * (1 - g ** (N - 1)) / (1 - g) ** 2
    )
    c0t = c0 / (1 - mu2)
    a = ((1 - g**2) * c0t) ** 0.5

    return g, a, mu2
