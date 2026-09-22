import numpy as np
from rich import traceback

def read_grid_dta(infile):

    # Open the file in binary read mode
    with open(infile, "rb") as f:

        # Read the three 8-byte integers (big-endian)
        a,b, i8_imt1, c, i8_jmt1, d, i8_km1 = np.fromfile(f, dtype='>i4', count=7)

        # Read the arrays
        a      = np.fromfile(f, dtype='>f8', count=1)
        r8_dxt = np.fromfile(f, dtype='>f8', count=i8_imt1)
        r8_dyt = np.fromfile(f, dtype='>f8', count=i8_jmt1)
        r8_dxu = np.fromfile(f, dtype='>f8', count=i8_imt1)
        r8_dyu = np.fromfile(f, dtype='>f8', count=i8_jmt1)
        r8_dzt = np.fromfile(f, dtype='>f8', count=i8_km1)
        r8_dzw = np.fromfile(f, dtype='>f8', count=i8_km1+1)  # 0:km in Fortran
        r8_xt  = np.fromfile(f, dtype='>f8', count=i8_imt1)
        r8_xu  = np.fromfile(f, dtype='>f8', count=i8_imt1)
        r8_yt  = np.fromfile(f, dtype='>f8', count=i8_jmt1)
        r8_yu  = np.fromfile(f, dtype='>f8', count=i8_jmt1)
        r8_zt  = np.fromfile(f, dtype='>f8', count=i8_km1)
        r8_zw  = np.fromfile(f, dtype='>f8', count=i8_km1)

    return r8_dxt, r8_dyt, r8_dxu, r8_dyu, r8_dzt, r8_dzw, r8_xt, r8_xu, r8_yt, r8_yu, r8_zt, r8_zw


def main():
    traceback.install()

    infile="/S/data01/G4006/y0330/data/ofes_exp/settings/box_bounded/topog/grid.dta.out"
    #infile="/S/data01/G4006/y0330/data/ofes_exp/settings/box_bnd_deep/topog/grid.dta.out"

    dxt, dyt, dxu, dyu, dzt, dzw, xt, xu, yt, yu, zt, zw = read_grid_dta(infile)

    # Print info
    print("read_grid")
    print('  len(xt), len(yt), len(zt) = ',len(xt), len(yt), len(zt))
    print('  dxt = ',  dxt[0],', ',  dxt[1],', ..., ',  dxt[-1])
    print('  dyt = ',  dyt[0],', ',  dyt[1],', ..., ',  dyt[-1])
    print('  dxu = ',  dxu[0],', ',  dxu[1],', ..., ',  dxu[-1])
    print('  dyu = ',  dyu[0],', ',  dyu[1],', ..., ',  dyu[-1])
    print('  dzt = ',  dzt[0],', ',  dzt[1],', ..., ',  dzt[-1])
    print('  dzw = ',  dzw[0],', ',  dzw[1],', ..., ',  dzw[-1])
    print('  xt  = ',  xt[ 0],', ',  xt[ 1],', ..., ',  xt[ -1])
    print('  xu  = ',  xu[ 0],', ',  xu[ 1],', ..., ',  xu[ -1])
    print('  yt  = ',  yt[ 0],', ',  yt[ 1],', ..., ',  yt[ -1])
    print('  yu  = ',  yu[ 0],', ',  yu[ 1],', ..., ',  yu[ -1])
    print('  zt  = ',  zt[ 0],', ',  zt[ 1],', ..., ',  zt[ -1])
    print('  zw  = ',  zw[ 0],', ',  zw[ 1],', ..., ',  zw[ -1])
    print()

if __name__ == "__main__":
    main()
