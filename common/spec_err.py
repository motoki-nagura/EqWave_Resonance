from scipy.stats import chi2

def spec_err(edof,alpha):
    """
    Compute error for power spectrum

 

    Input edof :  Effective degree of freedom
          alpha : Percentage for significance
 
    Output : value_low, value_high
 
       edof should be determined based on
         Table 5.6.4 in Emery and Thomson (when the use of a digital filter)
           or pp. 450-451 in Emery and Thomson (when the use of band averaging)
         If a box averaging of 'nw' grid points is applied to the spectrum (in frequency
           space) 'nm' times, edof = 2*(nw**nm) (pp.450-451 in Emery and Thomson)
    
       alpha should be 0 < alpha < 1. Usually, alpha = 0.95.
 
       Actual error bar is computed as power *value_low or power *value_high
         (It depends on spectral power; See Eq. (5.6.92) in Emery and Thomson)
 
 
    --- Sample program to plot error bars ---
 
        x1, y1, alpha = freq[i], vps[i], 0.95
        relerr_low, relerr_high = spec_err(edof, alpha)
        err_low, err_hig = y1 *relerr_low, y1 *relerr_hig
        axs[iplot].plot(x1,y1,'ko')
        axs[iplot].plot([x1,x1],[err_low,err_hig],c='k')
        axs[iplot].text(x1*1.1, y1, f'{alpha_spec*1e2:3.1f}%')
 
    --- Consistency check for Table D.2 in Emergy and Thomson ---
        print(chi2.ppf(q=0.025, df=10))  #   3.2469
        print(chi2.ppf(q=0.500, df=20))  #  19.337
        print(chi2.ppf(q=0.950, df=25))  #  37.652
        print(chi2.ppf(q=0.010, df=25))  #  11.523
    """
    alpha1 = 1. - alpha   #  see one line below Eq.(5.6.92) in Emery and Thomson

    value_low = edof /chi2.ppf(q=1.-alpha1/2., df=edof)   #  Eq.(5.6.92) in Emery and Thomson
    value_hig = edof /chi2.ppf(q=alpha1/2.,    df=edof)

    return value_low, value_hig
#}}}


