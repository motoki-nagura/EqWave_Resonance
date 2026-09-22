import matplotlib.pyplot as plt

def new_page(nrows,ncols,figsize=(8,8),stack="row major"):
    iplot = 0
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=figsize)
    if nrows*ncols > 1:
        if   stack == "row major":
            axes = axes.flatten()
        elif stack == "column major":
            axes = axes.T.flatten()
    return fig, axes, iplot
