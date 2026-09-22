"""
    Supportive functions to draw figures
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, TwoSlopeNorm

def lon2x_180dy(x):
    return x-70.+30.
def x2lon_180dy(x):
    return x-30.+70.
def lon2x_90dy(x):
    return x-85.+45.
def x2lon_90dy(x):
    return x-45.+85.


def generate_axes():
    fig, axes = plt.subplots(nrows=6, ncols=1, figsize=(8,6))

    ax_ul = fig.add_axes([0.09, 0.73, 0.380, 0.22])  #  [left, bottom, width, height]
    ax_ml = fig.add_axes([0.09, 0.41, 0.380, 0.22])
    ax_bl = fig.add_axes([0.09, 0.05, 0.305, 0.22])
    ax_ur = fig.add_axes([0.57, 0.73, 0.380, 0.22])
    ax_mr = fig.add_axes([0.57, 0.41, 0.380, 0.22])
    ax_br = fig.add_axes([0.57, 0.05, 0.305, 0.22])

    for i in range(len(axes)):
        fig.delaxes(axes[i])
    axes = [ax_ul, ax_ml, ax_bl, ax_ur, ax_mr, ax_br]

    return fig, axes


def make_smooth_seismic(data, width=0.25, n_colors=256):
    """
    Create a seismic-like colormap with a smooth transition to pure white at zero.

    Parameters
    ----------
    data : array-like
        Data used to determine the normalization range.
    width : float, default=0.25
        Width of the region blended toward white, expressed as a fraction
        of the colormap half-range. Larger values produce a smoother and
        broader transition around zero.
    n_colors : int, default=256
        Number of colors used to construct the colormap.

    Returns
    -------
    cmap : ListedColormap
        Modified seismic colormap.
    norm : TwoSlopeNorm
        Normalization centered at zero.
    """
    # Sample colors from the original seismic colormap
    colors = plt.colormaps["seismic"](np.linspace(0, 1, n_colors))

    # Coordinates centered at zero
    x = np.linspace(-1, 1, n_colors)

    # Compute the blending weight toward white
    # The weight is 1 at zero and decreases linearly to 0
    w = np.clip(1.0 - np.abs(x) / width, 0.0, 1.0)

    # Blend the original RGB values with pure white
    colors[:, :3] = (
        colors[:, :3] * (1.0 - w[:, None])
        + w[:, None]
    )

    cmap = ListedColormap(colors, name="smooth_seismic")

    # Center the normalization at zero
    vmax = np.nanmax(np.abs(data))
    norm = TwoSlopeNorm(vmin=vmax*(-1.), vcenter=0.0, vmax=vmax)

    return cmap, norm
