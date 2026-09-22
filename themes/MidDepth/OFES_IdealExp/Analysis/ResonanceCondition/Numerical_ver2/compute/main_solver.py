"""
 Solve the resonance conditions.

   ----------------------------------
            Descriptions
   ----------------------------------

     The resonance conditions are

          ω1 + ω2 = ω3,
          k1 + k2 = k3,
          m1 + m2 = m3,
          n1 + n2 + n3 = odd,
          Ω1(n1,ω1,k1,m1) = 0,
          Ω2(n2,ω2,k2,m2) = 0,
          Ω3(n3,ω3,k3,m3) = 0,

     where

       n : Meridional mode
       ω : Frequency [rad s^-1]
       k : Zonal wavenumber [m s^-1]
       m : Vertical wavenumber [m s^-1]
       Ω(n,ω,k,m) = 0 : Dispersion relation

     In this program, I prescribe n1, n2, n3, ω1, ω2, ω3, and k1. First I solve

          Ω1(n1,ω1,k1,m1) = 0,

     and obtain m1. Then the remaining conditions are written as

          Ω2(n2,ω2,k2,m2) = 0,
          Ω3(n3,ω3,k3,m3) = 0.

    where k3 = k1 + k2 and m3 = m1 + m2. I solve these equations by a root finding
    method in terms of k2 and m2.

   ----------------------------------
            Input parameters
   ----------------------------------

     Input parameters are defined in ./params/params_XXX.py. A list is

                 N2 : Scalar
                 ω1 : Scalar
                 ω2 : Scalar
                 ω3 : Scalar
         n1_grids   : 1D array
         n2_grids   : 1D array
         n3_grids   : 1D array
         k1_grids   : 1D array
         k2_guesses : 1D array
         m1_guess   : Scalar
         m2_guesses : 1D array

   ----------------------------------
              Output
   ----------------------------------

     The results are saved in 'OUTFILE_NETCDF' and write out as texts in 'OUTFILE_TEXT'.

   ----------------------------------
                Note
   ----------------------------------

     n_jobs = 4  :   MacBook Air
     n_jobs = 62 :   EA and ES-Moon

"""
import numpy as np
from root_solver_utils import (
        solve_system, remove_duplicate_solutions,
        check_k_sign, scaling, unscaling, compute_m1, is_sum_zero,
        DipersionRelationTerms, )
import time
import sys
from joblib import Parallel, delayed
from tqdm import tqdm
from tqdm_joblib import tqdm_joblib

from importlib.util import spec_from_file_location, module_from_spec
from pathlib import Path
import os

def main():
    #=========================================================
    #                  User defined parameters
    #=========================================================
    #   Options:  tag                directory                     out_XXX.nc
    settings = { "std":            ("N2_5e-6/180-90-60dy/n1_1"  , ""),
                 "std,wk":         ("N2_5e-6/180-90-60dy/n1_1"  , "_wide_k1"),
                 "std,m1<0":       ("N2_5e-6/180-90-60dy/n1_1"  , "_negative_m1"),
                 "std,k1~=0":      ("N2_5e-6/180-90-60dy/n1_1"  , "_near0_k1"),
                 "n1=-1":          ("N2_5e-6/180-90-60dy/n1_m1" , ""),
                 "n111":           ("N2_5e-6/180-90-60dy/n111"  , ""),
                 "n111,sp":        ("N2_5e-6/180-90-60dy/n111"  , "_specific"),
                 "std,112":        ("N2_5e-6/180-180-90dy/n1_1" , ""),
                 "std,112,m1<0":   ("N2_5e-6/180-180-90dy/n1_1" , "_negative_m1"),
                 "std,112,n1=-1":  ("N2_5e-6/180-180-90dy/n1_m1", ""),
                 "std,112,n1=3":   ("N2_5e-6/180-180-90dy/n1_3" , ""),
                 "140":            ("N2_5e-6/140-70-47dy/n1_1"  , ""),
                 "220":            ("N2_5e-6/220-110-55dy/n1_1" , ""),
                }["std,112,n1=3"]
    #                        n1_1 : n1=1;   n1_m1 : n1=-1

    #   Paths
    PATH_INPUT_PARAM = Path("params/"+settings[0]+"/params"+settings[1]+".py")

    current_dir = os.getcwd()                      #  Get current directory path
    OUT_DIR = current_dir.replace("code", "data")  #  Replace "code" with "data"

    OUTFILE_NETCDF = OUT_DIR+"/results/"+settings[0]+"/out"+settings[1]+".nc"
    OUTFILE_TEXT   = OUT_DIR+"/results/"+settings[0]+"/out"+settings[1]+".txt"

    #=======================================================
    #                   Fixed Parameters
    #=======================================================
    day_in_sec = 86400.           #  1 day in seconds
    Erad       = 6371e3           #  Earth's radius [m]
    EOmega     = 7.292e-5         #  Earth's rotation rate [s^-1]
    beta       = 2.*EOmega /Erad  #  Meridional gradient of the Coriolis coefficient [s^-1 m^-1]

    start_time = time.time()
    txtout = []

    #=========================================================
    #                    Read in Parameters
    #=========================================================

    spec = spec_from_file_location("params_module", PATH_INPUT_PARAM)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    
    P = module.get_params()

    N2,ω1,ω2,ω3 = P.N2, P.ω1, P.ω2, P.ω3
    n1_values,n2_values,n3_values = P.n1_values,P.n2_values,P.n3_values
    k1_values,k2_guesses = P.k1_values,P.k2_guesses
    m1_guess, m2_guesses = P.m1_guess, P.m2_guesses

    txtout.append("==>  Parameters are read from "+settings[0]+" and "+settings[1]+"  <==")
    print(txtout[-1])

    if not is_sum_zero(ω1,ω2,-ω3, verbose=False):
        raise ValueError("ω1 + ω2 = ω3 is not satisfied! Check params_XXX.py.")

    #    Non dimensionalize
    ω1,ω2,ω3, k1_values, m1_guess, k2_guesses, m2_guesses= \
        scaling(ω1, ω2, ω3, k1_values, m1_guess, k2_guesses, m2_guesses, beta,N2)

    #=========================================================
    #              Solve the Resonance Conditions
    #=========================================================

    # Build all tasks: each task is one (k1,m1,n1,n2,n3) combination
    tasks = []
    for k1 in k1_values:
        for n1 in n1_values:
            k1_checked = check_k_sign(n1, ω1, k1)
            m1 = compute_m1(ω1, k1_checked, n1, m1_guess)
            for n2 in n2_values:
                for n3 in n3_values:
                    if (n1 + n2 + n3) % 2 != 0:
                        tasks.append((k1_checked, m1, n1, n2, n3))

    #    Parameters in an array
    initial_points = np.array([[k2, m2] for k2 in k2_guesses for m2 in m2_guesses])

    # Define function to process one task
    def process_task(k1, m1, n1, n2, n3):
        solutions = solve_system(initial_points, k1, ω2, ω3, n2, n3, m1)
        N = solutions.shape[0]
        if N == 0:
            return None
        return {
            'k2': solutions[:,0],
            'm2': solutions[:,1],
            "k1": np.full(N, k1),
            "m1": np.full(N, m1),
            "n1": np.full(N, n1),
            "n2": np.full(N, n2),
            "n3": np.full(N, n3),
        }

    # Parallel execution with progress bar
    with tqdm_joblib(tqdm(desc="Solving systems", total=len(tasks))) as progress_bar:
        results = Parallel(n_jobs=62)(delayed(process_task)(*task) for task in tasks)

    # Collect results into lists
    k1_list, k2_list          = [], []
    m1_list, m2_list          = [], []
    n1_list, n2_list, n3_list = [], [], []

    for res in results:
        if res is not None:
            k2_list.extend(res['k2'])
            m2_list.extend(res['m2'])
            k1_list.extend(res['k1'])
            m1_list.extend(res['m1'])
            n1_list.extend(res['n1'])
            n2_list.extend(res['n2'])
            n3_list.extend(res['n3'])

    ## ===> Single CPU <=== Loop over all parameter combinations
    #for k1 in k1_values:
    #    for n1 in n1_values:
    #        m1 = compute_m1(ω1,k1,n1,m1_guess) #    Compute m1
    #        k1 = check_k_sign(n1,ω1,k1)
    #        for n2 in n2_values:
    #            for n3 in n3_values:
    #                if n1+n2+n3 % 2 != 0:    #   n1 + n2 + n3 = odd
    #                   solutions = solve_system(initial_points, k1, ω2, ω3, n2, n3, m1)
    #                   #
    #                   N = solutions.shape[0]
    #                   if N > 0:
    #                       # Add k2 and m2
    #                       k2_list.extend(solutions[:,0])
    #                       m2_list.extend(solutions[:,1])
    #                       # Add parameters with length N
    #                       k1_list.extend([k1]*N)
    #                       m1_list.extend([m1]*N)
    #                       n1_list.extend([n1]*N)
    #                       n2_list.extend([n2]*N)
    #                       n3_list.extend([n3]*N)

    #  Convert to NumPy arrays
    k1_out, k2_out         = np.array(k1_list), np.array(k2_list)
    m1_out, m2_out         = np.array(m1_list), np.array(m2_list)
    n1_out, n2_out, n3_out = np.array(n1_list), np.array(n2_list), np.array(n3_list)

    #  Compute k3 and m3 from the resonance conditions
    k3_out = k1_out + k2_out
    m3_out = m1_out + m2_out

    #  Rename
    k1, k2, k3 = k1_out, k2_out, k3_out
    m1, m2, m3 = m1_out, m2_out, m3_out
    n1, n2, n3 = n1_out, n2_out, n3_out

    #  Remove duplicates
    solutions = np.column_stack((k1,k2,k3,m1,m2,m3,n1,n2,n3))
    solutions = remove_duplicate_solutions(solutions)
    k1,k2,k3,m1,m2,m3,n1,n2,n3 = solutions.T

    #==========================================================================
    #             Check the solutions and Output to a text file
    #==========================================================================
    print("==>  Check the solutions and Output to a text file")

    #  Dimensionalize
    ω1,ω2,ω3, k1,k2,k3, m1,m2,m3 = unscaling(ω1,ω2,ω3, k1,k2,k3, m1,m2,m3, beta,N2)


    txtout.append('')
    txtout.append('==>  Fixed parameters  <==')
    txtout.append(f'                   N^2 = {N2:.3e} s^-2')
    omegas = np.array([   ω1,   ω2,   ω3]) *1e6
    Ts     = np.array([1./ω1,1./ω2,1./ω3]) *2.*np.pi /day_in_sec
    txtout.append("      ω1,    ω2,    ω3 = "+", ".join(f"{x:>7.2f}" for x in omegas)+"  x 10^-6 s^-1")
    txtout.append("   2π/ω1, 2π/ω2, 2π/ω3 = "+", ".join(f"{x:>7.2f}" for x in Ts )+"  days")
    txtout.append('')

    txtout.append('==>  Search range  <==')
    k1_range  = np.array([np.min(k1_values), np.max(k1_values)])
    Lx2_range = np.array([np.min(2.*np.pi/k2_guesses), np.max(2.*np.pi/k2_guesses)]) /1e3
    Lz2_range = np.array([np.min(2.*np.pi/m2_guesses), np.max(2.*np.pi/m2_guesses)])
    txtout.append("      k1 = "+" ~ ".join(f"{x:>7.2f}" for x in  k1_range)+" rad m^-1")
    txtout.append("   2π/k2 = "+" ~ ".join(f"{x:>7.2f}" for x in Lx2_range)+" km")
    txtout.append("   2π/m2 = "+" ~ ".join(f"{x:>7.2f}" for x in Lz2_range)+" m" )
    txtout.append('')

    for i in range(len(k2)):
        txtout.append(f'==>  Solution {i+1:2d}  <==')

        ωs  = np.array([   ω1   ,    ω2   ,    ω3   ])
        ks  = np.array([   k1[i],    k2[i],    k3[i]])
        ms  = np.array([   m1[i],    m2[i],    m3[i]])
        ns  = np.array([   n1[i],    n2[i],    n3[i]])
        Lxs = np.array([1./k1[i], 1./k2[i], 1./k3[i]]) *2.*np.pi
        Lzs = np.array([1./m1[i], 1./m2[i], 1./m3[i]]) *2.*np.pi

        ns = np.int32(ns)

        #----------------------------------
        #  Format each number, same width
        #----------------------------------
        ks_scaled  = ks *1e6
        ms_scaled  = ms *1e3
        Lxs_scaled = Lxs /1e3
        txtout.append("   n1, n2, n3 = "+", ".join(f"{x:>5d}" for x in ns))
        txtout.append("   k1, k2, k3 = "+", ".join(f"{x:>10.8f}" for x in ks_scaled )+"  x 10^-6 m^-1")
        txtout.append("   m1, m2, m3 = "+", ".join(f"{x:>10.8f}" for x in ms_scaled )+"  x 10^-3 m^-1")
        txtout.append('')
        txtout.append("   2π/k1, 2π/k2, 2π/k3 = "+", ".join(f"{x:>7.1f}" for x in Lxs_scaled)+"  km")
        txtout.append("   2π/m1, 2π/m2, 2π/m3 = "+", ".join(f"{x:>7.1f}" for x in Lzs       )+"  m")
        txtout.append('')

        #---------------------------------------------------
        #  Check if the resonance conditions are satisfied
        #---------------------------------------------------
        L_res_w = is_sum_zero(ωs[0],ωs[1], -1.*ωs[2], verbose=False)
        L_res_k = is_sum_zero(ks[0],ks[1], -1.*ks[2], verbose=False)
        L_res_m = is_sum_zero(ms[0],ms[1], -1.*ms[2], verbose=False)
        L_res_n = np.sum(ns) % 2 == 1

        #---------------------------------------------------
        #  Check if the dispersion relations are satisfied
        #---------------------------------------------------
        L_dis = np.zeros(3,dtype=bool)
        for j,(w,k,m,n) in enumerate(zip(ωs,ks,ms,ns)):
            term1,term2,term3,term4 = DipersionRelationTerms(w,k,m,n,N2,beta)
            L_dis[j] = is_sum_zero(term1,term2,term3,term4, verbose=False)

        #-----------------
        #  Print results
        #-----------------
        txtout.append('   Is ω1 + ω2 = ω3 satisfied? '+str(L_res_w))
        txtout.append('   Is k1 + k2 = k3 satisfied? '+str(L_res_k))
        txtout.append('   Is m1 + m2 = m3 satisfied? '+str(L_res_m))
        txtout.append('   Is n1 + n2 + n3 = odd satisfied? '+str(L_res_n))
        txtout.append('   Is the dispersion relation for mode 1 satisfied? '+str(L_dis[0]))
        txtout.append('   Is the dispersion relation for mode 2 satisfied? '+str(L_dis[1]))
        txtout.append('   Is the dispersion relation for mode 3 satisfied? '+str(L_dis[2]))
        txtout.append('')
        if not np.all([L_res_w,L_res_k,L_res_m,L_res_n,L_dis[0],L_dis[1],L_dis[2]]):
            txtout.append('   Some of the conditions are not met... This solution has been discarded!')
            txtout.append('')

        #-----------------------
        #  Save good solutions
        #-----------------------
        if np.all([L_res_w,L_res_k,L_res_m,L_res_n,L_dis[0],L_dis[1],L_dis[2]]):
            out_vars = np.array([k1[i],k2[i],k3[i], m1[i],m2[i],m3[i], n1[i],n2[i],n3[i]])
            try:
                out_to_netcdf = np.append(out_to_netcdf, out_vars)
            except NameError:  # if arr doesn't exist yet
                out_to_netcdf = out_vars

    #------------------------------------------------------------------------------------------
    #  Array is like out_to_netcdf[{# of solution},{parameters (k1,k2,k3,m1,m2,m3,n1,n2,n3)}]
    #------------------------------------------------------------------------------------------
    try:
        out_to_netcdf = out_to_netcdf.reshape(-1,9) #  collect every 9 elements
    except NameError:
        out_to_netcdf = np.zeros((1,9)) *np.nan   #  if no valid solution was obtained, create a dummy

    #==========================================================================
    #                           Output to a NetCDF file
    #==========================================================================
    import netCDF4
    print("==>  Output to a netcdf file")

    #---------------------
    #   collect results
    #---------------------
    nsol = len(out_to_netcdf[:,0])

    ks, ms, ns = np.zeros((3,3,nsol))

    for i in range(nsol):
        ks[:,i] = out_to_netcdf[i,0:3]
        ms[:,i] = out_to_netcdf[i,3:6]
        ns[:,i] = out_to_netcdf[i,6:9]
    ns = ns.astype(int)

    ωs = np.array([ω1,ω2,ω3])

    #---------------------
    #     write out
    #---------------------
    f1 = netCDF4.Dataset(OUTFILE_NETCDF,'w',format='NETCDF4')

    f1.createDimension('p',       1)
    f1.createDimension('triad',   3)
    f1.createDimension('solution',nsol)

    N2_1   = f1.createVariable('N2',   'float32', ('p'))
    beta_1 = f1.createVariable('beta', 'float32', ('p'))
    ωs_1   = f1.createVariable('omega', ωs.dtype, ('triad'))
    ks_1   = f1.createVariable('k',     ks.dtype, ('triad','solution'))
    ms_1   = f1.createVariable('m',     ms.dtype, ('triad','solution'))
    ns_1   = f1.createVariable('n',     ns.dtype, ('triad','solution'))

    setattr(f1,'Description','ω1, ω2, ω3, k1, n1, n2, n3 are prescribed')

    N2_1.setncattr(  'long_name','Buoyancy frequency squared')
    beta_1.setncattr('long_name','Meridional gradient of the Coriolis coefficient')
    ns_1.setncattr(  'long_name','Meridional modes')
    ωs_1.setncattr(  'long_name','Frequency')

    N2_1.setncattr(  'units','s^-2')
    beta_1.setncattr('units','s^-1 m^-1')
    ωs_1.setncattr(  'units','rad s^-1')
    ks_1.setncattr(  'units','rad m^-1')
    ms_1.setncattr(  'units','rad m^-1')

    N2_1[0]   = N2
    beta_1[0] = beta
    ωs_1[:]   = ωs[:]
    ks_1[:,:] = ks[:,:]
    ms_1[:,:] = ms[:,:]
    ns_1[:]   = ns[:,:]

    f1.close()

    #==========================================================================
    #                           Close the program
    #==========================================================================
    end_time = time.time()
    txtout.append(f"Elapsed time: {end_time - start_time:.3f} seconds")

    with open(OUTFILE_TEXT,'w') as fa:
        for text in txtout:
            fa.write(text)
            fa.write('\n')

# Entry point
if __name__ == "__main__":
    sys.exit(main())
