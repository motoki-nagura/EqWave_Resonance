"""
   Plot the maximum of the real part of bispectrum at each (k3, m3) grid
"""
from rich import traceback
import numpy as np
import sys, os
from dataclasses import dataclass

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import AutoMinorLocator
import matplotlib.patches as patches
from matplotlib.patches import PathPatch 
import matplotlib.colors as colors
from matplotlib.path import Path
from shapely.geometry import box, Polygon, MultiPolygon
from shapely.geometry.polygon import orient
from shapely.ops import unary_union

from common.new_page import (new_page,)
from common.sub_NetCDF_IO import (read_variable_from_netcdf,)
from themes.MidDepth.OFES_IdealExp.Tools.common.Overplot_DispersionRelation import (
        Overplot_DispersionRelation,)

π = np.pi
day2sec = 86400.   #  1 day in seconds
ω_180dy, ω_90dy, ω_60dy = 2.*π/(180.*day2sec), 2.*π/(90.*day2sec), 2.*π/(60.*day2sec)   #  [rad second^-1]


#==========================================================================================
#     Input Parameters
#==========================================================================================

#  Model run

EXP_NAME = { "std" :  "180dy+90dy",
             "sens" : "180dy",     }["std"]
#             "sens" : "180dy",     }["sens"]

#  Frequency set

FREQ_NAME = {"ABC":  "180-90-60day",
             "AAB":  "180-180-90day"}["ABC"]
#             "AAB":  "180-180-90day"}["AAB"]

#  Filter

FILTER = "3ptBox"
#FILTER = "NoSmth"

#  Component

COMPONENT = {"xyz": "xyz",  #  u ∂u/dx + v ∂u/∂y + w ∂u/∂z
             "xz":  "xz",   #  u ∂u/dx + w ∂u/∂z
             "x":   "x",    #  u ∂u/dx
             "y":   "y",    #  v ∂u/∂y
             "z":   "z",}   #  w ∂u/∂z
COMPONENT = COMPONENT["x"]

#==========================================================================================
#     Derived Parameters
#==========================================================================================

#  Input

match FREQ_NAME:
    case "180-90-60day":  tmp_label = "60dy"
    case "180-180-90day": tmp_label = "90dy"

#PATH_IN_LIST = \
#    "../k3m3_region/"+FREQ_NAME+"/results/x/out_list_"+tmp_label+"_Rossby_3ptBox.nc"

INFILE_EXT = { "NoSmth" : "NoSmth",
               "3ptBox" : "3ptBox",
             }[FILTER]

current_dir = os.getcwd()                        #  Get current directory path
INDIR_DATA = current_dir.replace("code", "data") #  Replace "code" with "data"

DIR_SUB = "/results/Exp_"+EXP_NAME+"/"+FREQ_NAME+"/"+COMPONENT+"/"

INFILE_NETCDF = INDIR_DATA+DIR_SUB+"out_"+INFILE_EXT+".nc"

#  Output

OUTFILE_FIG = "./"+DIR_SUB+"fig_"+INFILE_EXT+".pdf"

#  Set ω1, ω2, and k3-m3 region for Rossby

@dataclass(slots=True)
class BoxP:
    k3_inc1: list[float] | None = None  # float list
    m3_inc1: list[float] | None = None
    k3_exc1: list[float] | None = None
    m3_exc1: list[float] | None = None
    k3_exc2: list[float] | None = None
    m3_exc2: list[float] | None = None

boxP = BoxP()

if   FREQ_NAME == "180-90-60day":
    ω1, ω2, ttl_ω3 = ω_180dy, ω_90dy , "60 days"   #  ω3 is 2π/(60 days)
    boxP.k3_inc1, boxP.m3_inc1 = [-1.8e-5, 0.0e-5], [-0.4e-2, 0.4e-2] # Target region for Rossby

elif FREQ_NAME == "180-180-90day":
    ω1, ω2, ttl_ω3 = ω_180dy, ω_180dy, "90 days"   #  ω3 is 2π/(90 days)
    boxP.k3_inc1, boxP.m3_inc1 = [-2.0e-5, 0.0e-5], [-0.8e-2, 0.8e-2] # Target region for Rossby
    boxP.k3_exc1, boxP.m3_exc1 = [-0.5e-5, 0.0e-5], [ 0.5e-2, 0.8e-2] # Exclude this region
    boxP.k3_exc2, boxP.m3_exc2 = [-0.7e-5, 0.0e-5], [ 0.6e-2, 0.8e-2]


#==========================================================================================
def main():
    traceback.install()

    #==========================================================================================
    #  Read In
    #==========================================================================================
    BSR_k3m3, coords = read_variable_from_netcdf(
            INFILE_NETCDF, "BSR_k3m3", return_coords=True, )

    m_BSR = coords["m_BSR"]
    k_BSR = coords["k_BSR"]

    #==========================================================================================
    #  Maximum value on each (m3,k3) grid
    #==========================================================================================
    BSR_k3m3_max = np.nanmax(BSR_k3m3, axis=0)

    ##==========================================================================================
    ##              Read in Top N values list
    ##==========================================================================================
    #print("IN << ",PATH_IN_LIST)
    #f1 = netCDF4.Dataset(PATH_IN_LIST, 'r')
    #i_k1 = np.copy(f1.variables['index_k1'])
    #i_m1 = np.copy(f1.variables['index_m1'])
    #i_k2 = np.copy(f1.variables['index_k2'])
    #i_m2 = np.copy(f1.variables['index_m2'])
    #f1.close()

    #k3_top = k_BSR[i_k1] + k_BSR[i_k2]
    #m3_top = m_BSR[i_m1] + m_BSR[i_m2]

    #==========================================================================================
    #                    Figure
    #==========================================================================================

    ttl, ω, xlab, ylab = ttl_ω3, ω1+ω2, "$k_3$ (m$^{-1}$)", "$m_3$ (m$^{-1}$)"

    k_range, m_range = [-2e-5,2e-5], [-4e-2,4e-2]
    xcd,ycd, xax_range,yax_range = k_BSR,m_BSR, k_range,m_range

    xcd_ext = np.append(xcd, xcd[-1]+(xcd[-1]-xcd[-2]))

    nrows, ncols = 2, 2

    #----------

    print('OUT >>> '+OUTFILE_FIG)

    with PdfPages(OUTFILE_FIG) as pdf:

        #-------------------------------------------------------------------------------------
        #-------------------------------------------------------------------------------------
        fig, axes, iplot = new_page(nrows, ncols)

        #-------------------------------------------------------------------------------------
        #   Maximum of real bispectrum

        clab = 'Real Bispectrum'
        scale = np.nanmax(BSR_k3m3)
        arr   = BSR_k3m3_max/scale

        v0, cmap, ticks = 1.0, 'RdBu_r', [-1.0, -0.5, 0.5, 1.0]
        norm = colors.SymLogNorm( linthresh=1e0, linscale=1e0, vmin=-v0, vmax=v0)

        #------
        ax = axes[iplot]
        #
        pcm = ax.pcolormesh(xcd,ycd,arr, cmap=cmap, shading='auto', norm=norm)
        fig.colorbar(pcm, ax=ax, label=clab, location='bottom', ticks=ticks)
        #
        ax = set_ax(ax, "", xax_range, yax_range, xlab,ylab)
        #ax = set_ax(ax, "Maximum, "+ttl, xax_range, yax_range, xlab,ylab)
        Overplot_DispersionRelation(np.abs(ω), xcd_ext, ax)
        Draw_Box(ax, boxP)
        #ax.plot(k3_top,m3_top,'.',ms=1,c='w')    #    Mark top values with white dots
        #
        iplot += 1

        #-------------------------------------------------------------------------------------
        fig.delaxes(axes[1])
        fig.delaxes(axes[2])
        fig.delaxes(axes[3])

        plt.tight_layout()
        pdf.savefig()
        plt.close(fig)

        #-------------------------------------------------------------------------------------
        #-------------------------------------------------------------------------------------
        fig, axes, iplot = new_page(nrows, ncols)

        #-------------------------------------------------------------------------------------
        #   Maximum

        ax = axes[iplot]
        arr = BSR_k3m3_max/scale
        #
        pcm = ax.pcolormesh(xcd,ycd,arr, cmap=cmap, shading='auto', norm=norm)
        fig.colorbar(pcm, ax=ax, label=clab, location='bottom', ticks=ticks)
        #
        ax = set_ax(ax, "Maximum, "+ttl, xax_range, yax_range, xlab,ylab)
        Overplot_DispersionRelation(np.abs(ω), xcd_ext, ax)
        #
        iplot += 1

        #-------------------------------------------------------------------------------------
        #   Minimum

        ax = axes[iplot]
        arr = np.nanmin(BSR_k3m3, axis=0)/scale
        #
        pcm = ax.pcolormesh(xcd,ycd,arr, cmap=cmap, shading='auto', norm=norm)
        fig.colorbar(pcm, ax=ax, label=clab, location='bottom', ticks=ticks)
        #
        ax = set_ax(ax, "Minimum, "+ttl, xax_range, yax_range, xlab,ylab)
        Overplot_DispersionRelation(np.abs(ω), xcd_ext, ax)
        #
        iplot += 1

        #-------------------------------------------------------------------------------------
        #   Average

        ax = axes[iplot]
        arr = np.nanmean(BSR_k3m3, axis=0)/scale
        #
        pcm = ax.pcolormesh(xcd,ycd,arr, cmap=cmap, shading='auto', norm=norm)
        fig.colorbar(pcm, ax=ax, label=clab, location='bottom', ticks=ticks)
        #
        ax = set_ax(ax, "Average, "+ttl, xax_range, yax_range, xlab,ylab)
        Overplot_DispersionRelation(np.abs(ω), xcd_ext, ax)
        #
        iplot += 1

        #-------------------------------------------------------------------------------------
        #   Mean of positive

        v0 = 1e-2
        arr = BSR_k3m3/scale
        norm = colors.SymLogNorm( linthresh=1e0, linscale=1e0, vmin=-v0, vmax=v0)

        ax = axes[iplot]
        arr1 = np.nanmean(np.where(arr >= 0., arr, np.nan), axis=0)
        #
        pcm = ax.pcolormesh(xcd,ycd,arr1, cmap=cmap, shading='auto', norm=norm)
        fig.colorbar(pcm, ax=ax, label=clab, location='bottom')
        #
        ax = set_ax(ax, "Mean of Positive Values, "+ttl, xax_range, yax_range, xlab,ylab)
        Overplot_DispersionRelation(np.abs(ω), xcd_ext, ax)
        #
        iplot += 1

        #-------------------------------------------------------------------------------------
        plt.tight_layout()
        pdf.savefig()
        plt.close(fig)


def attach_Xaxis_ticks(label,ax):
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.xaxis.set_ticks(np.arange(-2.0e-5, 2.1e-5, 1.0e-5))


def attach_Yaxis_ticks(label,ax):
    ax.yaxis.set_ticks(np.arange(-0.04, 0.041, 0.02))
    ax.yaxis.set_minor_locator(AutoMinorLocator(4))


def set_ax(ax, ttl, xax_range, yax_range, xlab,ylab):
    ax.set_title(ttl)
    ax.set_xlim(xax_range)
    ax.set_ylim(yax_range)
    ax.axvline(x=0.,lw=0.5,c='k')
    ax.axhline(y=0.,lw=0.5,c='k')
    attach_Xaxis_ticks(xlab,ax)
    attach_Yaxis_ticks(ylab,ax)
    ax.set(xlabel=xlab, ylabel=ylab)
    return ax


def Add_Rectangle(ax,krange,mrange):
    rect = patches.Rectangle(
           (np.min(krange), np.min(mrange)),
           np.max(krange)-np.min(krange), np.max(mrange)-np.min(mrange),
           linestyle=':',linewidth=1, edgecolor='black', facecolor='none' )
    ax.add_patch(rect)
    return ax


def Draw_Box(ax, boxP):
    outer = box(boxP.k3_inc1[0], boxP.m3_inc1[0], boxP.k3_inc1[1], boxP.m3_inc1[1])

    if boxP.k3_exc1 == None:
        inner = None
    else:
        inner1 = box( boxP.k3_exc1[0], boxP.m3_exc1[0], boxP.k3_exc1[1], boxP.m3_exc1[1] )
        inner2 = box( boxP.k3_exc2[0], boxP.m3_exc2[0], boxP.k3_exc1[1], boxP.m3_exc2[1] )
        inner = [inner1, inner2]

    patches = shapely_difference_patch( outer, inner, facecolor="none", edgecolor="black",
                                        linewidth=1, linestyle=":" )
    for patch in patches:
        ax.add_patch(patch)


def shapely_difference_patch(
    outer,
    hole_polygons,
    **patch_kwargs
):
    """
    Subtract hole_polygons from outer, and then return PathPatch

    Parameters
    ----------
    outer : shapely.geometry.Polygon
        The outer polygon

    hole_polygons : list of shapely.geometry.Polygon
        Polygons to subtract from the outer polygon

    **patch_kwargs
        argument to give to PathPatch
        (facecolor, edgecolor, linewidth, etc.)

    Returns
    -------
    list[PathPatch]
        List of PathPatch to draw figures
    """

    if hole_polygons:
        holes_union = unary_union(hole_polygons)
        geometry = outer.difference(holes_union)
    else:
        geometry = outer

    if geometry.is_empty:
        return []

    if isinstance(geometry, Polygon):
        polygons = [geometry]
    elif isinstance(geometry, MultiPolygon):
        polygons = list(geometry.geoms)
    else:
        raise TypeError(
            f"Unsupported geometry type: {geometry.geom_type}"
        )

    patches = []

    for polygon in polygons:
        # outer polygon is anti-clockwise; holes are clockwise
        polygon = orient(polygon, sign=1.0)

        vertices = []
        codes = []

        def add_ring(ring):
            coords = np.asarray(ring.coords, dtype=float)

            if len(coords) < 4:
                return

            vertices.extend(coords)
            codes.extend(
                [Path.MOVETO]
                + [Path.LINETO] * (len(coords) - 2)
                + [Path.CLOSEPOLY]
            )

        add_ring(polygon.exterior)

        for interior in polygon.interiors:
            add_ring(interior)

        path = Path(
            np.asarray(vertices),
            np.asarray(codes)
        )

        patch = PathPatch(path, **patch_kwargs)
        patches.append(patch)

    return patches


#  Script entry point
if __name__ == "__main__":
    sys.exit(main())
