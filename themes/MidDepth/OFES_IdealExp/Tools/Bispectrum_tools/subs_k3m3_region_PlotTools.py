import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from common.new_page import new_page

def manage_page(iplot,nrows,ncols,fig,pdf,axs):
    iplot += 1
    if iplot == nrows*ncols:
        plt.tight_layout()
        pdf.savefig()
        plt.close(fig)          #  Close
        fig, axs, iplot = new_page(nrows, ncols)   #  Open the new one
    return fig,pdf,axs, iplot


def attach_Xaxis_ticks(ax):
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.xaxis.set_ticks(np.arange(-2.0e-5, 2.1e-5, 1.0e-5))


def attach_Yaxis_ticks(ax):
    ax.yaxis.set_ticks(np.arange(-0.04, 0.041, 0.02))
    ax.yaxis.set_minor_locator(AutoMinorLocator(4))


def set_ax(ax, ttl, xax_range, yax_range, xlab,ylab, label_fontsize=None):
    ax.set_title(ttl, fontsize=9)
    ax.set_xlim(xax_range)
    ax.set_ylim(yax_range)
    ax.axvline(x=0.,lw=0.5,c='k')
    ax.axhline(y=0.,lw=0.5,c='k')
    attach_Xaxis_ticks(ax)
    attach_Yaxis_ticks(ax)
    ax.set_xlabel(xlab, fontsize=label_fontsize)
    ax.set_ylabel(ylab, fontsize=label_fontsize)
    return ax

