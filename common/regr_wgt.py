"""
      Regression coefficients with weights
"""
import numpy as np
from scipy.stats import t as stut

def regr_wgt ( x, y, sig, mwt, edof, alpha ):
#-------------------------------------------------------------------------------------
#   Citation : Numerical Recipes in Fortan 77, by W.H. Press, S.A. Teukolsky,
#                   W.T. Vetterling, B.P. Flannery, Cambridge University Press, 1992
#
#
#   Given a set of data points x(1:ndata), y(1:ndata) with individual standard
#    deviations sig(1:ndata), fit them to a straight line y = a + bx by minimizing
#    chi squared. Returned are a, b and their respective probable uncertainties
#    siga and sigb, the chi-square chi2, and the goodness-of-fit probability q
#    (that the fit would have chi squared this large or larger). If mwt = 0 on
#    input, then the standard deviations are assumed to be unavailable: q is returned
#    as 1.0 and the normalization of chi2 is to unit standard deviation on all points.
#
#
#    INPUT
#
#      x : Real(1:ndata)        Given data points
#      y : Real(1:ndata)        Given data
#    sig : Real(1:ndata)        Standard deviations of the given data
#    mwt : Logical              If False, the "sig" is unavailable.
#                                 If True, the "sig" is available.
#   edof : Real(1)              Effective degree of freedom
#                                  (Should be estimated by eff_DOF_cor or eff_DOF_trend)
#  alpha : Real(1)              Significance level ( 0. < alpha < 1. )
#
#
#    OUTPUT
#
#      a : Real(1)             Intersect of the regression line
#      b : Real(1)             Slope of the regression line
#   errb : Real(1)             Standard error of the slope "b", computed with Eq.(3.12.9)
#                                in "Data analysis methods in Physical Oceanography" authored
#                                by Emery and Thomson, 2004, pages 238, 261-262.
#
#-------------------------------------------------------------------------------------

  x = np.array(x).astype('double')
  y = np.array(y).astype('double')

  ndata = np.size(x)

  sx=0.                                  #   Initialize sums to zero
  sy=0.
  st2=0.
  b=0.

  if (mwt):                              # Accumulate sums ...
    ss = 0.
    for i in range(ndata):                     # ... with weights
      wt = 1./(sig[i]**2)
      ss = ss + wt
      sx = sx + x[i]*wt
      sy = sy + y[i]*wt
  else:
    sx = np.sum(x)                     # .. or without weights.
    sy = np.sum(y)
    ss = float(ndata)

  sxoss = sx/ss

  if (mwt):
    for i in range(ndata):
      t = (x[i]-sxoss)/sig[i]
      st2 = st2 + t*t
      b = b + t*y[i]/sig[i]
  else:

    st2 = np.sum((x-sxoss)**2)
    b = np.sum((x-sxoss) *y)

  b = b/st2                                        #  Solve for a, b, sigma_a and sigma_b.
  a = (sy-sx*b)/ss

# -----  Error estimate

  avex = np.sum(x) / float(ndata)
  sx2 = np.sum( (x-avex)**2 ) /float(ndata-1)

  sse = np.sum( ( y - (b*x + a) )**2 )
  s_epsilon = np.sqrt( sse/(float(ndata)-2.) )

  alpha1 = 1. - (1.-alpha)*0.5   #  two sided
  errb = s_epsilon * stut.ppf ( alpha1, edof ) / np.sqrt( (edof-1.) * sx2 )

  if errb < 0.:
    print('s_epsilon = ',s_epsilon)
    print('Students_T ( alpha, edof ) = ',stut.ppf ( alpha, edof ))
    print('sqrt((edof-1.) * sx2 ) = ',np.sqrt((edof-1.) * sx2 ))

  return a,b,errb
