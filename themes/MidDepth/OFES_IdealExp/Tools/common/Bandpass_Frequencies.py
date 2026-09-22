def Bandpass_Frequencies(Tper, df_log10 = 0.125):
    """ Compute low-cut and high-cut frequencies
             Tper     :  Target period [day]
             df_log10 :  Separation in log_10 frequency
    """
    df = (10**df_log10 -1.) /(10**df_log10 +1.) *1./Tper
    low_cut, high_cut = 1./Tper-df, 1./Tper+df
    return low_cut, high_cut
