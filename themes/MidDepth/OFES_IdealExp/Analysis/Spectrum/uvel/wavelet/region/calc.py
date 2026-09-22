"""
  Compute wavelet power averaged over specific period ranges in a limited spatial domain
"""
from rich import traceback
traceback.install()

import numpy as np
import xarray as xr
import netCDF4
from tqdm import tqdm
from joblib import Parallel, delayed

import os, pathlib, sys

from common.wavelet_subroutine import run_wavelet_analysis
from themes.MidDepth.OFES_IdealExp.Tools.common.Bandpass_Frequencies import (
        Bandpass_Frequencies, )


#=================================
#          Parameters
#=================================

EXP_NAME = "180dy+90dy"
#EXP_NAME = "180dy"
#EXP_NAME = "90dy"
#EXP_NAME = "BottomDamp_180dy"

MY_WORK_ROOT = str( pathlib.Path(os.environ["MY_WORK_ROOT"]) )

INPUT_FILE = MY_WORK_ROOT+"/data/ofes_exp/results/box_bounded/compiled/"+\
         "u_eq_box_bounded_"+EXP_NAME+"_2000-2029.nc"

#------
T_target = [180., 90., 60.]    #   Target periods

x0, x1 = 20., 40.              #   Longitudinal range
j0     = 0                     #   Latitudinal index
z0, z1 = 2000., 3200.          #   Depth range

S0_FACTOR = 5.     #  Parameter for wavelet analysis

df_log10 = {"default" : 0.125,   #    Separation in log_10 frequency
            "wide" :    0.175,
            "narrow" :  0.100,} ["default"]

#------
current_dir = os.getcwd()                      #  Get current directory path
OUT_DIR = current_dir.replace("code", "data")  #  Replace "code" with "data"

ext = ""
if df_log10 != 0.125:
    ext = f"_df{df_log10:.3f}".replace(".","p")

OUTPUT_FILE = OUT_DIR+"/out_"+EXP_NAME+ext+".nc"


#=================================
#          Functions
#=================================
def main():
    #=================================
    #       Input Time Series
    #=================================
    print('IN << '+INPUT_FILE)
    var  = xr.open_dataset(INPUT_FILE)['u'].\
            isel(lat=j0).sel(lon=slice(x0,x1),lev=slice(z0,z1))

    lon, lev = np.array(var.lon), np.array(var.lev)
    var      = np.array(var)

    print(f"var size: {var.nbytes/1024**3:.2f} GB")
    print(f"var == nan: {int(np.sum(np.isnan(var))):d}")
    print(f"var == inf: {int(np.sum(np.isinf(var))):d}")

    #   Pad 0 at the beginning in time
    Npad = 2000
    var = np.pad(
        var,
        pad_width=((Npad, 0), (0, 0), (0, 0)),  # Pad 0 N times in time
        mode='constant',
        constant_values=0
        )

    time = np.arange(0., len(var[:,0,0])) - Npad  #  time in days
    dt = 1.0                                      #  sampling interval [days]

    #=================================
    #       Find Period Bands
    #=================================
    Tmin, Tmax = np.zeros((2,len(T_target)))

    for i,T in enumerate(T_target):
        lc, hc = Bandpass_Frequencies(T, df_log10=df_log10)
        Tmin[i] = 1./hc
        Tmax[i] = 1./lc

    #=================================
    #          RUN SUBROUTINE
    #=================================

    #--------------------------------------
    #          Parallel version
    #--------------------------------------

    # ---- Pre-compute sizes ----
    nt = var.shape[0]

    # Test call to get output sizes
    test = run_wavelet_analysis(
            var[:, 0, 0], dt, s0_factor=S0_FACTOR, Tmin=Tmin, Tmax=Tmax)

    nper = len(test.period)

    # Allocate output arrays
    power_avg_out = np.zeros((nper, nt))
    scale_avg_out = np.full((len(Tmin), nt, len(lev), len(lon)), np.nan)

    period = None
    coi    = None
    sig95  = None

    # ---- Worker function ----
    def wavelet_worker(k, i):
        res = run_wavelet_analysis(
            var[:, k, i],
            dt,
            s0_factor=S0_FACTOR,
            Tmin=Tmin,
            Tmax=Tmax
        )
        return k, i, res.power, res.scale_avg, res.period, res.coi, res.sig95

    # ---- Parallel execution with progress bar ----
    tasks = [(k, i) for i in range(len(lon)) for k in range(len(lev))]

    results = Parallel(n_jobs=-1)(
        delayed(wavelet_worker)(k, i)
        for k, i in tqdm(tasks, desc="Wavelet analysis")
    )

    # --- Retrieve results ---
    Nsum = 0

    for k, i, power, scale_avg, p, c, s in results:
        power_avg_out += power
        Nsum += 1

        scale_avg_out[:, :, k, i] = scale_avg

        # Save representative values once
        if period is None:
            period = p
        if coi is None:
            coi = c
        if sig95 is None:
            sig95 = s

    # ---- Normalize ----
    power_avg_out /= Nsum

    ##--------------------------------------
    ##          Unparallel version {{{
    ##--------------------------------------
    #
    ## Pre-compute sizes
    #nt = var.shape[0]
    #
    ## Run one test call to know output sizes
    #test = run_wavelet_analysis(var[:,0,0], dt, s0_factor=5., Tmin=Tmin, Tmax=Tmax)
    #
    #nper = len(test.period)
    #
    ## Allocate output arrays
    #power_avg_out = np.full((nper, nt), 0.)
    #scale_avg_out = np.full((len(Tmin), nt, len(lev), len(lon)), np.nan)
    #
    #period = None
    #coi    = None
    #sig95  = None
    #
    ## ---- execution ----
    #for i in range(len(lon)):
    #    for k in range(len(lev)):
    #        print("i,k = ",i,k)
    #
    #        res = run_wavelet_analysis(var[:,k,i],dt,s0_factor=5., Tmin=Tmin, Tmax=Tmax)
    #
    #        # ---- Store results ----
    #        power_avg_out[:,:]    += res.power   #    sum up
    #        scale_avg_out[:,:,k,i] = res.scale_avg
    #
    #        # Save representative values
    #        if period is None:
    #            period = res.period
    #        if coi is None:
    #            coi = res.coi
    #        if sig95 is None:
    #            sig95 = res.sig95
    #
    #
    #Ngrids = len(lon) *len(lev)
    #power_avg_out /= Ngrids      #   Divide by the number of grid points
    #}}}

    #======================================
    #      Write Out to a NetCDF file
    #======================================

    print("OUT >> "+OUTPUT_FILE)
    f1 = netCDF4.Dataset(OUTPUT_FILE,'w',format='NETCDF4')

    f1.createDimension('lon',   len(lon))
    f1.createDimension('lev',   len(lev))
    f1.createDimension('time',  len(time))
    f1.createDimension('period',len(period))
    f1.createDimension('Tband', len(T_target))
    f1.createDimension('param', 1)

    lon_1       = f1.createVariable('lon',       'float32',('lon'))
    lev_1       = f1.createVariable('lev',       'float32',('lev'))
    time_1      = f1.createVariable('time',      'float64',('time'))
    power_avg_1 = f1.createVariable('power_avg', 'float32',('period','time'))
    scale_avg_1 = f1.createVariable('scale_avg', 'float32',('Tband', 'time','lev','lon'))
    period_1    = f1.createVariable('period',    'float32',('period'))
    coi_1       = f1.createVariable('coi',       'float32',('time'))
    sig95_1     = f1.createVariable('sig95',     'float32',('period','time'))
    T_target_1  = f1.createVariable('T_target',  'float32',('Tband'))
    x0_1        = f1.createVariable('x0',        'float32',('param'))
    x1_1        = f1.createVariable('x1',        'float32',('param'))
    z0_1        = f1.createVariable('z0',        'float32',('param'))
    z1_1        = f1.createVariable('z1',        'float32',('param'))

    time_1.setncattr('units','Days since 2000-01-01')

    lon_1[          :] = lon
    lev_1[          :] = lev
    time_1[         :] = time
    power_avg_1[  :,:] = power_avg_out
    scale_avg_1[:,:,:] = scale_avg_out
    period_1[       :] = period
    coi_1[          :] = coi
    sig95_1[      :,:] = sig95
    T_target_1[     :] = T_target
    x0_1[           0] = x0
    x1_1[           0] = x1
    z0_1[           0] = z0
    z1_1[           0] = z1

    f1.close()

# Script entry point
if __name__ == "__main__":
    sys.exit(main())
