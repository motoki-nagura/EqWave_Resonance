import numpy as np
from scipy.optimize import root
from scipy.spatial import cKDTree

def F_system(vars, k1, ω2, ω3, n2, n3, m1): #{{{
    """
    Compute the system of equations F2,F3 given variables (k2,m2)
    """
    k2, m2 = vars
    m3 = m1 + m2
    k3 = k1 + k2

    # F2
    if   n2 == -1:  # Kelvin mode
        F2 = k2 - np.abs(m2) * ω2
    elif n2 ==  0:  # Yanai mode
        F2 = ω2 *k2 - np.abs(m2) *ω2**2 + 1.
    else:           # Rossby/IG modes
        F2 = ω2 * k2**2 + k2 - m2**2 * ω2**3 + np.abs(m2) * ω2 * (2*n2+1)

    # F3
    if   n3 == -1:  # Kelvin mode
        F3 = k3 - np.abs(m3) * ω3
    elif n3 ==  0:  # Yanai mode
        F3 = ω3 *k3 - np.abs(m3) *ω3**2 + 1.
    else:           # Rossby/IG modes
        F3 = ω3 * k3**2 + k3 - m3**2 * ω3**3 + np.abs(m3) * ω3 * (2*n3+1)

    return np.array([F2, F3])  # Only two equations to match two variables
#}}}

def filter_initial_points(initial_points, k1, ω2, ω3, n2, n3, m1, Rmax=100): #{{{
    """
    Filter initial points based on residual magnitude
    """
    filtered_points = []
    for pt in initial_points:
        res = F_system(pt, k1, ω2, ω3, n2, n3, m1)
        R = np.linalg.norm(res)
        if R <= Rmax:
            filtered_points.append(pt)
    return np.array(filtered_points)
#}}}

def remove_duplicate_solutions(solutions, tol=1e-5): #{{{
    """
    Remove duplicate solutions using numpy
    """
    tree = cKDTree(solutions)
    seen = np.zeros(len(solutions), dtype=bool)
    unique_idx = []

    for i in range(len(solutions)):
        if not seen[i]:
            # Find all rows close to solutions[i]
            neighbors = tree.query_ball_point(solutions[i], tol)
            seen[neighbors] = True
            unique_idx.append(i)

    return solutions[unique_idx]
#}}}

def solve_system(initial_points, k1, ω2, ω3, n2, n3, m1): #{{{
    """
    Solve the system for multiple initial points with filtering and duplicate removal
    """
    # Filter points based on residual
    filtered_points = filter_initial_points(initial_points, k1, ω2, ω3, n2, n3, m1)

    solutions = []
    for pt in filtered_points:
        sol = root(F_system, pt, args=(k1, ω2, ω3, n2, n3, m1), method="lm", tol=1e-20)
        #sol = root(F_system, pt, args=(k1, ω2, ω3, n2, n3, m1))
        if sol.success:
            solutions.append(sol.x)

    solutions = np.array(solutions)

#    if len(solutions) != 0:
#        solutions = remove_duplicate_solutions(solutions)

    return solutions
#}}}

def check_k_sign(n,ω,k): #{{{
    if   n == -1:    #  Kelvin
        if k != 0 and ω/k < 0:
            k = -1. *k    #  k flipped such that ω/k > 0 
    elif n >= 0:     #  Other modes
        if k != 0 and ω/k > 0:
            k = -1. *k    #  k flipped such that ω/k < 0 (considering Rossby)
    return k
#}}}

def scaling(ω1, ω2, ω3, k1_values, m1_guess, k2_guess, m2_guess, beta,N2): #{{{
    #    Make variables non dimensional
    #
    c0 = 1.0                         #   arbitrary constant
    ω_scale = np.sqrt(beta *c0)      #   Frequency scale
    k_scale = np.sqrt(beta /c0)      #   Zonal wavenumber scale
    m_scale = np.sqrt(N2) /c0        #   Vertical wavenumber scale

    ω1 /= ω_scale        #   non dimensional variables
    ω2 /= ω_scale
    ω3 /= ω_scale

    k1_values /= k_scale

    m1_guess /= m_scale

    k2_guess /= k_scale
    m2_guess /= m_scale

    return ω1, ω2, ω3, k1_values, m1_guess, k2_guess, m2_guess
#}}}

def unscaling(ω1,ω2,ω3, k1,k2,k3, m1,m2,m3, beta,N2): #{{{
    #
    #    Make variables dimensional
    #
    c0 = 1.0                         #   arbitrary constant
    ω_scale = np.sqrt(beta *c0)      #   Frequency scale
    k_scale = np.sqrt(beta /c0)      #   Zonal wavenumber scale
    m_scale = np.sqrt(N2) /c0        #   Vertical wavenumber scale

    ω1, ω2, ω3 = ω1*ω_scale, ω2*ω_scale, ω3*ω_scale
    k1, k2, k3 = k1*k_scale, k2*k_scale, k3*k_scale
    m1, m2, m3 = m1*m_scale, m2*m_scale, m3*m_scale

    return ω1,ω2,ω3, k1,k2,k3, m1,m2,m3
#}}}

def compute_m1(ω1,k1,n1,m1_guess): #{{{
    #-------------------------------------------------------------------------------
    #      Compute m1 by solving the dispersion relation for mode 1
    #            (m1 is assumed to be positive.)
    #-------------------------------------------------------------------------------
    def fun0(vars):
        m1 = vars
        # Dispersion relation for Yanai mode (non dimensional)
        F1 = ω1 *k1 - np.abs(m1) *ω1**2 + 1.
        return F1

    def fun1(vars):
        m1 = vars
        # Dispersion relation for Rossb/IG mode (non dimensional)
        F1 = ω1 *k1**2 + k1 - m1**2 *ω1**3 + np.abs(m1) *ω1 *(2.*n1+1.)
        return F1

    #  Compute m1
    if   n1 == -1:      #   Kelvin
        m1 = k1 /ω1

    elif n1 == 0:       #   Yanai
        x0 = m1_guess
        res = root(fun0, x0, method="lm", tol=1e-20)
        m1 = res.x[0]

    else:               #   Rossby/IG
        x0 = m1_guess
        res = root(fun1, x0, method="lm", tol=1e-20)
        m1 = res.x[0]

        #coeff = np.array([ ω1**2 /N2,
        #                   -1. *beta /np.sqrt(N2) *(2.*n1 + 1.),
        #                   -1. *beta *k1 /ω1 - k1**2] )
        #m = np.roots(coeff)
        #m = np.sort(np.abs(m))
        #m1 = m[0]

    return m1
#}}}

def is_sum_zero(*args, rtol=1e-5, verbose=False): #{{{
    """
    Check whether the sum of any number of variables is approximately zero
    (using absolute and relative tolerance). Optionally prints details.

    Parameters
    ----------
    *args : float or array-like
        Values to sum.
    rtol : float, optional
        Relative tolerance (default 1e-5)
    verbose : bool, optional
        If True, print diagnostic info when the condition fails.

    Returns
    -------
    bool or ndarray of bool
        True where |sum(args)| <= rtol * max(|args|).
    """
    import numpy as np
    arrs = np.array(args)
    total = np.sum(arrs, axis=0)
    ref = np.max(np.abs(arrs), axis=0)
    ok = np.abs(total) <= (rtol * ref)

    # Print diagnostic info if condition fails and verbose=True
    if verbose and not np.all(ok):
        print("Condition failed: |sum(args)| > rtol * max(|args|)")
        print("-" * 60)
        # Print each term
        for i, a in enumerate(args, 1):
            print(f"Term {i}: {a: .6e}")
        print(f"  Sum : {total: .6e}")
        print(f"  rtol: {rtol: .1e}, scale: {ref: .6e}")
        print(f"  Limit = {rtol * ref: .6e}")
        print("-" * 60)

    return ok
#}}}

def DipersionRelationTerms(ω,k,m,n,N2,beta): #{{{
    #------------------------------------------------------
    #    Return the terms of the dispersion relation
    #    *** Dimensional ***
    #------------------------------------------------------
    import numpy as np

    if   n == -1:             #  Kelvin
        term1 = k
        term2 = -1. *np.abs(m) *ω /np.sqrt(N2)
        term3 = 0e0
        term4 = 0e0

    elif n == 0:              #  Yanai
        term1 = ω *k
        term2 = -1. *np.abs(m) /np.sqrt(N2) *ω**2
        term3 = beta
        term4 = 0e0

    else:                     #  Rossby/IG
        term1 = ω *k**2
        term2 = beta *k
        term3 = -1. *m**2 *ω**3 /N2
        term4 = beta *np.abs(m) *ω /np.sqrt(N2) *(2*n+1)

    return term1,term2,term3,term4
#}}}

