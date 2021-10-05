import matplotlib as mpl
from matplotlib import pyplot as plt

def _subplots(ax=None, **kwargs):
    if ax is None:
        fig, ax = plt.subplots(**kwargs)
    else:
        # Get the current figure, create a new one if it doesn't exist.
        fig = plt.gcf()
    return fig,ax

def _plot_profile(wavelength, data, ax=None, **kwargs):
    fig, ax = _subplots(ax)
    ax.plot(wavelength, data, **kwargs)
    
