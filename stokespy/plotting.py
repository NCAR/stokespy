import matplotlib as mpl
from matplotlib import pyplot as plt

def _subplots(ax=None, **kwargs):
    if ax is None:
        fig, ax = plt.subplots(**kwargs)
    else:
        fig = plt.gcf()
    return fig, ax

def _plot_profile(wavelength, data, ax=None, **kwargs):
    fig, ax = _subplots(ax)
    ax.plot(wavelength, data, **kwargs)
