import numpy as np
from scipy import ndimage

def sign_flip__omega_extract__smoothing(var,ω,ω1,ω2, sign_flip=False, smoothing=None):
    """
      var[{ω1},{m1},{k1},{ω2},{m2},{k2}] :  Input
      ω                                  :  ω grid
      ω1                                 :  Target frequency of wave 1
      ω2                                 :  Target frequency of wave 2
      smoothing                          : = None     :  No Smoothing
                                         : = (3,3,3,3):  Boxcar filter
                                         : = "121"    :  1-2-1 filter for all the dimensions
                                         : = "12321"  :  1-2-3-2-1 filter
    """
    #-----------------------------------------------------------
    #   Sign is flipped. Positive = Energy gain. See note.
    #-----------------------------------------------------------
    if sign_flip:
        var_out = np.negative(var)
    else:
        var_out = np.copy(var)

    #--------------------------------------------------------------------------------------------
    #  Pick ω1 and ω2
    #
    #    Frequency is *Negative* here such that k < 0 means westward propagaiton
    #    (Usual Fourier transform converts du/dt to iω u)
    #
    #--------------------------------------------------------------------------------------------

    iω1 = np.argmin( np.abs(ω - (- ω1)) )
    iω2 = np.argmin( np.abs(ω - (- ω2)) )

    var_out = var_out[iω1,:,:, iω2,:,:]             #  ω' = ω1,  ω'' = ω2

    #-----------------------------------------------------------
    #   Smoothing over k1, m1, k2, and m2
    #-----------------------------------------------------------
    if smoothing != None:

        if smoothing is True:
            raise ValueError("Currently smoothing = True. Specify filter type.")

        elif ( isinstance(smoothing, tuple)
               and len(smoothing) == 4
               and all(isinstance(v, int) for v in smoothing)
             ):
            var_out = ndimage.uniform_filter(var_out, size=smoothing)    #   boxcar filter

        elif smoothing == "121":
            w = np.array([1, 2, 1], dtype=np.float64)   #  1-2-1 filter
            w /= w.sum()   # [0.25, 0.5, 0.25]

            for axis in range(var_out.ndim):
                var_out = ndimage.convolve1d(var_out, w, axis=axis, mode='reflect')

        elif smoothing == "12321":
            w = np.array([1,2,3,2,1], dtype=np.float64)   #  1-2-3-2-1 filter
            w /= w.sum()

            for axis in range(var_out.ndim):
                var_out = ndimage.convolve1d(var_out, w, axis=axis, mode='reflect')


    return var_out
