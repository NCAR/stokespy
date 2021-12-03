import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.widgets import Slider, Button
import numpy as np

#import matplotlib
#matplotlib.use('TkAgg')

def _subplots(ax=None, **kwargs):
    if ax is None:
        fig, ax = plt.subplots(**kwargs)
    else:
        # Get the current figure, create a new one if it doesn't exist.
        fig = plt.gcf()
    return fig,ax

def _plot_profile(wavelengths, data, ax=None, **kwargs):
    fig, ax = _subplots(ax, nrows=1, ncols=1, figsize=[4, 4], dpi=100)
    fig.subplots_adjust(bottom=0.1, top=0.9, left=0.15, right=0.90, wspace=0.0, hspace=0.0)
    
    ax.plot(wavelengths, data, **kwargs)
    
    ax.set_title('Stokes X')
    ax.set_xlabel('Wavelength')
    
def _plot_image(data, ax=None, proj=None, meta=None, **kwargs):
    fig, ax = _subplots(ax, nrows=1, ncols=1, figsize=[5, 5], dpi=120, subplot_kw={'projection':proj})
    fig.subplots_adjust(bottom=0.25, top=0.85, left=0.05, right=0.90, wspace=0.0, hspace=0.0)
    
    ax.imshow(data, **kwargs)
    
    if meta['wav1'] is None:
        ax.set_title('Stokes ' + meta['stokes'] + '\n Single image \n $\lambda$ = ' + str(meta['wav0']))
    else:
        ax.set_title('Stokes ' + meta['stokes'] + '\n Summed over \n $\lambda\in$ [' + str(meta['wav0']) + ' ' + str(meta['wav1']) + ']')

def _plot_3d_cube(wavelengths, data, ax=None, proj=None, meta=None, init=0, **kwargs):
    fig, ax = _subplots(ax, nrows=1, ncols=1, figsize=[5, 5], dpi=120, subplot_kw={'projection':proj})
    fig.subplots_adjust(bottom=0.25, top=0.85, left=0.05, right=0.9, wspace=0.0, hspace=0.0)
    
    img_plot = ax.imshow(data[init,:,:], **kwargs)
    ax.set_title('Stokes ' + meta['stokes'] + '\n Single image \n $\lambda$ = ' + str(wavelengths[init]))
    
    # Make a horizontal slider to control the frequency.
    axcolor = 'lightgoldenrodyellow'
    wav_ax = plt.axes([0.25, 0.1, 0.55, 0.06], facecolor=axcolor)
    allowed_vals = np.arange(0, len(wavelengths))
    
    wav_slider = Slider(ax=wav_ax, label='Wavelength',
            valmin=0, valmax=len(wavelengths)-1,
            valinit=init, valstep=allowed_vals)

    def f(ixt):
        return data[ixt,:,:]
    
    # The function to be called anytime a slider's value changes
    def update(val):
        img_plot.set_data(f(val))
        ax.set_title('Stokes ' + meta['stokes'] + '\n Single image \n $\lambda$ = ' + str(wavelengths[val]))
        fig.canvas.draw_idle()

    # register the update function with the slider
    wav_slider.on_changed(update)

    plt.show()
    
    return wav_slider

