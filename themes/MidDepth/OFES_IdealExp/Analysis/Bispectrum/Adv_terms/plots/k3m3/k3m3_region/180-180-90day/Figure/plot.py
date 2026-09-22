"""
     Draw figures
"""
from rich import traceback
traceback.install()
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import AutoMinorLocator
import matplotlib.patches as patches
import matplotlib.colors as colors
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
import string
import sys, pathlib, os

MY_WORK_ROOT = str( pathlib.Path(os.environ["MY_WORK_ROOT"]) )

from themes.MidDepth.OFES_IdealExp.Tools.common.subs_ofes_idealexp import (
        overplot_DispersionRelation,)
from themes.MidDepth.OFES_IdealExp.Tools.Bispectrum_tools.subs_k3m3_region_TargetRange import (
        km_range_180_90_60, km_range_180_180_90,)
from themes.MidDepth.OFES_IdealExp.Tools.Bispectrum_tools.subs_process_BS_ResTrio import (
        Read_Processed_Bispectrum_NaNAvg_Sum, Process_ResonantTrios_For_Cluster,)

#==============================================================================
#       Fixed Parameters
#==============================================================================

π = np.pi
day2sec = 86400.              #  1 day in seconds
Erad    = 6371e3              #  Earth's radius [m]
deg2met = 2. *π *Erad /360.   #  1 deg in meters

#==============================================================================
#       Parameters for bispectrum
#==============================================================================

setting = "90dy_Rossby_3ptBox" 
COMPONENT = "x"

INPUT_BS_FILE = "../results/"+COMPONENT+"/out_"+setting+".nc"

k_range = [-2.0e-5, 2.0e-5]     #  Wave number range for figures
m_range = [-2.0e-2, 2.0e-2]

m_ticks = np.arange(-0.02, 0.021, 0.01)

ω_180dy, ω_90dy = 2.*π/(180.*day2sec), 2.*π/(90.*day2sec)   #  [rad second^-1]
ω1, ω2, ω3 = ω_180dy, ω_180dy, ω_90dy

#==============================================================================
#       Parameters for resonant triads
#==============================================================================

set_res_N2       = "N2_5e-6"
set_res_freq     = "180-180-90dy"
sets_res_n1       =("n1_1", "n1_1", "n1_3", )
sets_res_add     = ("", "_negative_m1", "", )
sets_res_cluster = ("C1", "C2", "C3",)

number_of_clusters = len(sets_res_cluster)

PATH_RES_OUTPUT = \
  MY_WORK_ROOT+"/data/themes/MidDepth/OFES_IdealExp/Analysis/ResonanceCondition/"+\
  "Numerical_ver2/"

#==============================================================================

def main():
    #vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
    #                   Read In
    #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    BSR_m1k1_avg, BSR_m2k2_avg, BSR_m3k3_sum, \
      kwvn, mwvn, scale, _, _, _, _ = \
      Read_Processed_Bispectrum_NaNAvg_Sum(
              INPUT_BS_FILE, cluster_considered=number_of_clusters)

    BSR_m1k1_avg[[0,1],:,:] = BSR_m1k1_avg[[1,0],:,:]    #  Swap Clusters 1 and 2
    BSR_m2k2_avg[[0,1],:,:] = BSR_m2k2_avg[[1,0],:,:]
    BSR_m3k3_sum[[0,1],:,:] = BSR_m3k3_sum[[1,0],:,:]

    k1_range_clus, m1_range_clus, k2_range_clus, m2_range_clus = \
            np.zeros((4,2,number_of_clusters))

    for iclus in range(number_of_clusters):
        k1_range_clus[:,iclus], m1_range_clus[:,iclus] = km_range_180_180_90("180", "$k_1$", iclus)
        k2_range_clus[:,iclus], m2_range_clus[:,iclus] = km_range_180_180_90("180", "$k_2$", iclus)

    km12_range_clus = {"k1" : k1_range_clus,  "m1" : m1_range_clus,
                       "k2" : k2_range_clus,  "m2" : m2_range_clus, }

    ω123, \
      m_filtered_C1, m_remaining_C1, k_filtered_C1, k_remaining_C1,  \
      m_filtered_C2, m_remaining_C2, k_filtered_C2, k_remaining_C2,  \
      m_filtered_C3, m_remaining_C3, k_filtered_C3, k_remaining_C3, = \
      Process_ResonantTrios_For_Cluster(
              set_res_N2, set_res_freq, sets_res_n1, sets_res_add, sets_res_cluster,
              km12_range_clus, PATH_RES_OUTPUT)

    #vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
    #                   Draw Figures
    #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    for i_cluster in range(number_of_clusters):

        OUTPUT_PDF = f"fig_C{i_cluster+1:d}.pdf"
        print('OUT >>> '+OUTPUT_PDF)
        with PdfPages(OUTPUT_PDF) as pdf:

            #  Set parameters
            clab, xcd,ycd, v0, cmap = \
                    "Real Bispectrum", kwvn,mwvn, 1.0, "RdBu_r"
            norm = colors.SymLogNorm(linthresh=1e0, linscale=1e0, vmin=-v0, vmax=v0)
            colorbar_ticks = [-1.0, -0.5, 0.5, 1.0]

            #  Set labels
            txt_labels = [f"({c})" for c in string.ascii_lowercase]
            icnt = 0

            #  Set panels
            fig = plt.figure(figsize=(8,6))
            gs = GridSpec(3, 3, height_ratios=[1, 0.05, 1], hspace=0.55, wspace=0.45)
            fig.subplots_adjust(left=0.10, right=0.95, bottom=0.10, top=0.95)

            axes = np.empty((2,3), dtype=object)

            #------------------------
            #     Plot bispectrum
            #------------------------
            for j, (ttl0, ω, xlab,ylab, arr) in enumerate(zip(
                    ("180","180","90"),
                    (ω1,ω2,ω3,),
                    ("$k_1$","$k_2$","$k_3$",),
                    ("$m_1$","$m_2$","$m_3$"),
                    (BSR_m1k1_avg[i_cluster,:,:]/scale,
                     BSR_m2k2_avg[i_cluster,:,:]/scale,
                     BSR_m3k3_sum[i_cluster,:,:]/scale),
                    )):
                #  Set ax
                i = 0
                axes[i,j] = fig.add_subplot(gs[i, j])
                ax = axes[i,j]

                #  Set title
                ttl = txt_labels[icnt]+"  "+f"Cluster {i_cluster+1}, "+ttl0+" days"
                icnt += 1

                #  Draw a figure
                pcm = ax.pcolormesh(xcd,ycd,arr, cmap=cmap, shading='auto', norm=norm)
                #
                ax = set_ax(ax, ttl, k_range, m_range,
                            xlab+r" (m$^{-1}$)",ylab+r" (m$^{-1}$)", m_ticks,
                            label_fontsize=8)
                #
                xcd1 = np.append(xcd, 2*xcd[-1]-xcd[-2])  #  linearly extend
                overplot_DispersionRelation(np.abs(ω), xcd1,ax, line_width=0.5)

                #  Specify k/m range by a green box
                if j in (0,1):
                    krng, mrng = km_range_180_180_90(ttl0, xlab, i_cluster) 
                    ax = Add_Rectangle(ax,krng,mrng, color='g')

            #  Color bar
            subgs = GridSpecFromSubplotSpec(
                1, 5, subplot_spec=gs[1, :], width_ratios=[3, 4, 4, 4, 3] )
            cax = fig.add_subplot(subgs[0, 1:4])
            #       Adjust the vertical position of the colorbar
            pos = cax.get_position()
            dy = 0.025  # Positive --> Moves up;  Negative --> Moves down
            cax.set_position([pos.x0, pos.y0 + dy, pos.width, pos.height])
            cb = fig.colorbar(pcm, cax=cax, orientation="horizontal", ticks=colorbar_ticks)
            cb.set_label(clab, fontsize=9)
            cb.ax.tick_params(labelsize=9)

            #------------------------------
            #     Plot resonant trios
            #------------------------------
            match i_cluster:
                case 0:
                    m_filtered   , m_remaining   , k_filtered   , k_remaining    = \
                    m_filtered_C1, m_remaining_C1, k_filtered_C1, k_remaining_C1
                case 1:
                    m_filtered   , m_remaining   , k_filtered   , k_remaining    = \
                    m_filtered_C2, m_remaining_C2, k_filtered_C2, k_remaining_C2
                case 2:
                    m_filtered   , m_remaining   , k_filtered   , k_remaining    = \
                    m_filtered_C3, m_remaining_C3, k_filtered_C3, k_remaining_C3

            for iω, ω in enumerate(ω123):
                match iω:
                    case 0:  ttl0, xlab,ylab = "180", "$k_1$", "$m_1$"
                    case 1:  ttl0, xlab,ylab = "180", "$k_2$", "$m_2$"
                    case 2:  ttl0, xlab,ylab =  "90", "$k_3$", "$m_3$"

                #  Set ax
                i = 1
                axes[i,iω] = fig.add_subplot(gs[i+1, iω])
                ax = axes[i,iω]

                k, m = k_remaining[iω,:], m_remaining[iω,:]
                ax.scatter(k, m, marker="o", c="0.7", s=6)  #    Resonant trios

                k, m = k_filtered[iω,:], m_filtered[iω,:]
                ax.scatter(k, m, marker="o", c='r', s=8)    #    Targetted resonant trios

                overplot_DispersionRelation(ω, kwvn, ax, line_width=0.5)
                                                            #    Dispersion Relation

                #    Axes
                ttl = txt_labels[icnt]+"  "+f"Cluster {i_cluster+1}, "+ttl0+" days"
                ax = set_ax(ax, ttl, k_range, m_range,
                            xlab+r" (m$^{-1}$)",ylab+r" (m$^{-1}$)", m_ticks,
                            label_fontsize=8)
                icnt += 1

                #   Additional specification of k/m range
                if iω in (0,1):
                    krng, mrng = km_range_180_180_90(ttl0, xlab, i_cluster) 
                    ax = Add_Rectangle(ax,krng,mrng, color='g')

                #if ttl0 == "90":
                #    i_cluster = 0
                #    k3_range, m3_range = km_range_180_90_60("90", i_cluster)
                #    ax = Add_Rectangle(ax,k3_range,m3_range, color='g', linestyle=':')

            #  Close the page
            pdf.savefig()
            plt.close(fig)

 
def set_ax(ax, ttl, xax_range, yax_range, xlab,ylab, yticks, label_fontsize=8):
    ax.set_title(ttl, fontsize=9)
    ax.set_xlim(xax_range)
    ax.set_ylim(yax_range)
    ax.axvline(x=0.,lw=0.5,c='k')
    ax.axhline(y=0.,lw=0.5,c='k')
    ax.set_xlabel(xlab, fontsize=label_fontsize)
    ax.set_ylabel(ylab, fontsize=label_fontsize)
    ax.yaxis.set_ticks(yticks)
    ax.tick_params(right=True,top=True,which='both')
    ax.xaxis.set_tick_params(labelsize=8)
    ax.yaxis.set_tick_params(labelsize=8)
    ax.xaxis.get_offset_text().set_fontsize(8)
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.yaxis.set_minor_locator(AutoMinorLocator(5))
    return ax


def Add_Rectangle(ax, x,y, linestyle='-', color='k'):
    if x == None:
        return ax

    xmin, xmax = np.min(x), np.max(x)
    ymin, ymax = np.min(y), np.max(y)
    rect = patches.Rectangle(
           (xmin, ymin), xmax-xmin, ymax-ymin,
           linestyle=linestyle, linewidth=1, edgecolor='black', facecolor='none',ec=color)
    ax.add_patch(rect)
    return ax


if __name__ == "__main__":
    sys.exit(main())
