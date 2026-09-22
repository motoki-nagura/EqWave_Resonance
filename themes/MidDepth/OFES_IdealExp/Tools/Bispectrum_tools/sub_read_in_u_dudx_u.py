import numpy as np
import xarray as xr
from dataclasses import dataclass, field

@dataclass
class Input_Parameters:
    π:          float = 0.0
    deg2met:    float = 0.0
    day2sec:    float = 0.0
    EXP_NAME:   str   = ""
    dir_in:     str   = ""
    extension:  str   = ""
    VAR_NAMES:  list[str] = field(default_factory=list)
    z0:         float = 0.0
    z1:         float = 0.0
    x0:         float = 0.0
    x1:         float = 0.0
    block_size: int   = 0
    l0:         int   = 0
    l1:         int   = 0
    ω_range:    list[float] = field(default_factory=list)
    m_range:    list[float] = field(default_factory=list)
    k_range:    list[float] = field(default_factory=list)


def Read_In_U_DUDX_U(p):
    """
       Read in u3, dudx1, u2
    """
    if p.VAR_NAMES[0] == "w":
        in_Var1 = p.dir_in+"w_onU"+p.extension
    else:
        in_Var1 = p.dir_in+p.VAR_NAMES[0]+p.extension
    in_Var2 = p.dir_in+p.VAR_NAMES[1]+p.extension
    in_Var3 = p.dir_in+p.VAR_NAMES[2]+p.extension

    print('IN << '+in_Var1)
    print('IN << '+in_Var2)
    print('IN << '+in_Var3)

    isel_kw = dict( time=slice(p.l0, p.l1) )
    sel_kw  = dict(  lon=slice(p.x0, p.x1), lev=slice(p.z0, p.z1), )

    var1 = xr.open_dataset(in_Var1)[p.VAR_NAMES[0]].isel(**isel_kw).sel(**sel_kw)
    var2 = xr.open_dataset(in_Var2)[p.VAR_NAMES[1]].isel(**isel_kw).sel(**sel_kw)
    var3 = xr.open_dataset(in_Var3)[p.VAR_NAMES[2]].isel(**isel_kw).sel(**sel_kw)

    lon  = np.array(var1.lon)
    lev  = np.array(var1.lev)

    var1 = np.squeeze( np.copy(var1) )
    var2 = np.squeeze( np.copy(var2) )
    var3 = np.squeeze( np.copy(var3) )

    #   Remove NaN in the x direction
    i = np.argwhere( np.logical_and( ~np.isnan(var1[0,0,:]), \
                     np.logical_and( ~np.isnan(var2[0,0,:]),
                                     ~np.isnan(var3[0,0,:]) )) )[:,0]
    var1 = var1[:,:,i]
    var2 = var2[:,:,i]
    var3 = var3[:,:,i]
    lon  = lon[i]

    #   Depth coordinate (positive upward)
    lev  = np.flip( -lev )
    var1 = np.flip(var1, axis=1)
    var2 = np.flip(var2, axis=1)
    var3 = np.flip(var3, axis=1)

    #   Check NaN
    if np.any(np.isnan(var1)):
        raise ValueError('np.sum(np.isnan(var1)) = ',np.sum(np.isnan(var1)),
                         ',  var1.size = ',var1.size)
    if np.any(np.isnan(var2)):
        raise ValueError('np.sum(np.isnan(var2)) = ',np.sum(np.isnan(var2)),
                         ',  var2.size = ',var2.size)
    if np.any(np.isnan(var3)):
        raise ValueError('np.sum(np.isnan(var3)) = ',np.sum(np.isnan(var3)),
                         ',  var3.size = ',var3.size)

    #   Check shape
    if var1.shape == var2.shape == var3.shape:
        pass
    else:
        raise ValueError("var1, var2, and var3 are not in the same size.")

    #   Grid intervals
    dt = 86400.                       #  seconds
    dx = (lon[1]-lon[0]) *p.deg2met   #  meters
    dz = lev[1]-lev[0]                #  meters

    return var1, var2, var3, dt,dz,dx
