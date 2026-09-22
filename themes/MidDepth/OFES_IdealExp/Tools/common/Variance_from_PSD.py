import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import simpson

def Variance_from_PSD(freq, Pxx, f_low, f_high):
    """
      Build the new frequency grid:
      - include f_low and f_high exactly
      - include all original freq points lying between them

      freq   : Array[:] : Frequency grid points
      Pxx    : Array[:] : Power spectral density
      f_low  : Scalar   : lower limit of frequency
      f_high : Scalar   : higher limit of frequency
    """
    mask = (freq > f_low) & (freq < f_high)
    freq_new = np.concatenate(([f_low], freq[mask], [f_high]))

    # ------------------------
    # Interpolate var onto the new grid
    # ------------------------
    interp_fun = interp1d(freq, Pxx, kind='linear', fill_value="extrapolate")
    tmp_Pxx_new = interp_fun(freq_new)

    # ------------------------
    # Compute the definite integral
    # ------------------------
    Variance = simpson(tmp_Pxx_new, x=freq_new)

    ## ------------------------
    ## Optionally, interpolate onto the target frequency
    ## ------------------------
    #Pxx[p,k,i] = interp_fun(1./T_target[p])

    ## ------------------------
    ## Check frequency range
    ## ------------------------
    #if show_message:
    #    tmp_list = 1./freq_new
    #    print(f'T_target = {int(T_target[p]):d}, '+\
    #           '1/freq_new = '+', '.join(f'{x:.2f}' for x in tmp_list))

    return Variance
