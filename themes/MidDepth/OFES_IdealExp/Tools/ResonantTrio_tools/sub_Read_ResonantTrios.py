import numpy as np
import netCDF4

def Read_ResonantTrios(INFILE_RESONANCE):

    #--------------------------------------------------------------------------------------------#
    #                          Read in the resonance condition solutions
    #--------------------------------------------------------------------------------------------#

    print("IN <<< "+INFILE_RESONANCE)

    f1 = netCDF4.Dataset(INFILE_RESONANCE, "r")
    N2   = np.copy(f1.variables["N2"][0])
    β    = np.copy(f1.variables["beta"][0])
    ω    = np.copy(f1.variables["omega"])
    kwvn = np.copy(f1.variables["k"])
    mwvn = np.copy(f1.variables["m"])
    nmer = np.copy(f1.variables["n"])
    f1.close()

    #--------------------------------------------------------------------------------------------#
    #                        Find triad groups in terms of meridional modes
    #--------------------------------------------------------------------------------------------#
    #    Find the index to the identical triad in terms of nmer

    vals, idx, inv = np.unique(nmer.T, axis=0, return_index=True, return_inverse=True)

    num_groups = len(vals)  # Build list of groups
    groups = [[] for _ in range(num_groups)]

    for col_idx, group_id in enumerate(inv):
        groups[group_id].append(col_idx)

    # Get pattern for each group

    patterns = [vals[g] for g in range(len(vals))]

    #print("Groups:", groups)
    print("Patterns:", ", ".join(
        str(p.astype(int).tolist()).replace(" ", "")
        for p in patterns
    ))

    return  ω, mwvn,kwvn, nmer
