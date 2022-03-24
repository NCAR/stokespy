import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.widgets import Slider, Button
import numpy as np
import astropy.units as u

#import matplotlib
#matplotlib.use('TkAgg')

def _subplots(ax=None, **kwargs):
    if ax is None:
        fig, ax = plt.subplots(**kwargs)
    else:
        # Get the current figure, create a new one if it doesn't exist.
        fig = plt.gcf()
    return fig,ax

def _plot_profile(wavelengths, data, plot_u, ax=None, meta=None, **kwargs):
    if ax is None:
        fig, ax = _subplots(ax, nrows=1, ncols=1, figsize=[4, 4], dpi=100)
        fig.subplots_adjust(bottom=0.15, top=0.9, left=0.15, right=0.90, wspace=0.0, hspace=0.0)
    
    plot_wav = wavelengths.to(plot_u)
    ax.plot(plot_wav.value, data, **kwargs)
    
    ax.set_title('Stokes ' + meta['stokes'])
    
    x0_str = 'x0 = ' + str(round(meta['x0'].value, 1)) + ' ' + meta['x0'].unit.to_string()
    y0_str = 'y0 = ' + str(round(meta['y0'].value, 1)) + ' ' + meta['y0'].unit.to_string()
    ax.text(0.23, 0.95, x0_str, horizontalalignment='center', verticalalignment='center', transform=ax.transAxes) 
    ax.text(0.23, 0.9, y0_str, horizontalalignment='center', verticalalignment='center', transform=ax.transAxes) 
    
    ax.ticklabel_format(useOffset=False)
    ax.set_xlabel('Wavelength [' + plot_wav[0].unit.to_string() + ']')
    
    return ax

def _plot_context_all_profiles(wavelengths, data, context_img, plot_u, init=0, proj=None, ax=None, meta=None, **kwargs):
    """
    ax: a list with axis instances.
    """
    
    fig = plt.figure(constrained_layout=False)
    gs0 = fig.add_gridspec(1, 2, left=0.1, right=0.95, width_ratios=[2, 1],
                          hspace=0.2, wspace=0.2)

    gs00 = gs0[0].subgridspec(2, 1, wspace=0.0, hspace=0.2, 
                              width_ratios=[1], height_ratios=[5, 1])
    gs01 = gs0[1].subgridspec(4, 1, hspace=0.0, wspace=0.0)
    
    context_ax = []
    context_ax.append(fig.add_subplot(gs00[0,0], projection=proj))
    context_ax.append(fig.add_subplot(gs00[1,0]))
    
    # Resize the slider axis. 
    l, b, w, h = context_ax[1].get_position().bounds
    scale_w = 0.7
    scale_h = 0.7
    context_ax[1].set_position([l + w*(1-scale_w)/2, b + h*(1-scale_w)/2, 
                                w*scale_w, h*scale_w])
    
    stokes_ax = []
    for a in range(4):
        stokes_ax.append(fig.add_subplot(gs01[a]))
    
    #if ax is None:
    #    fig, ax = _subplots(ax, nrows=4, ncols=1, figsize=[3, 5], dpi=100)
    #    fig.subplots_adjust(bottom=0.1, top=0.95, left=0.15, right=0.90, wspace=0.0, hspace=0.0)
    
    # Plot the context image and selected point location.
    img_plot = context_ax[0].imshow(context_img[init,:,:], origin='lower')
    print(meta['x0_pix'], meta['y0_pix'])
    context_ax[0].plot(meta['x0_pix'], meta['y0_pix'], '+r')
    
    # Plot the Stokes vectors.
    plot_wav = wavelengths.to(plot_u)
    
    stokes_ax[0].plot(plot_wav.value, data[0,:], **kwargs)
    stokes_ax[1].plot(plot_wav.value, data[1,:], **kwargs)
    stokes_ax[2].plot(plot_wav.value, data[2,:], **kwargs)
    stokes_ax[3].plot(plot_wav.value, data[3,:], **kwargs)

    stokes_ax[0].set_title('I', y=0.75, fontweight=700, loc='center')
    stokes_ax[1].set_title('Q', y=0.75, fontweight=700, loc='center')
    stokes_ax[2].set_title('U', y=0.75, fontweight=700, loc='center')
    stokes_ax[3].set_title('V', y=0.75, fontweight=700, loc='center')

    # Plot the wavelength position shown in the context image.
    wav_positions = []
    for ax in stokes_ax:
        y0, y1 = ax.get_ylim()
        wav_positions.append(ax.plot([plot_wav[init].value, plot_wav[init].value], [y0,y1], '--k')[0])
    
    context_ax[0].set_xlabel('Helioprojective Longitude')
    context_ax[0].set_ylabel('Helioprojective Latitude')
    context_ax[0].set_title('Stokes ' + meta['stokes'] + '\n $\lambda$ = ' + str(np.round(plot_wav[init],2)))
    
    #title_string = '(x0, y0) = (' + str(round(meta['x0'].value, 1)) + ', ' + \
    #                str(round(meta['y0'].value, 1)) + ') ' + meta['x0'].unit.to_string()
    title_string = 'x0 = ' + str(round(meta['x0'].value, 1)) + ' arcsec \n' + \
                'y0 = ' + str(round(meta['y0'].value, 1)) + ' arcsec'
    #x0_str = 'x0 = ' + str(round(meta['x0'].value, 1)) + ' ' + meta['x0'].unit.to_string()
    #y0_str = 'y0 = ' + str(round(meta['y0'].value, 1)) + ' ' + meta['y0'].unit.to_string()
    #ax[0].text(0.2, 1.22, x0_str, horizontalalignment='left', verticalalignment='center', transform=ax[0].transAxes) 
    #ax[0].text(0.2, 1.1, y0_str, horizontalalignment='left', verticalalignment='center', transform=ax[0].transAxes) 
    
    stokes_ax[0].text(0.02, 1.03, title_string, fontsize=10, horizontalalignment='left', verticalalignment='bottom', transform=stokes_ax[0].transAxes)
    
    stokes_ax[3].ticklabel_format(useOffset=False)
    stokes_ax[3].set_xlabel('Wavelength [' + plot_wav[0].unit.to_string() + ']')
    
    # Make a horizontal slider to control the frequency.
    axcolor = 'lightgoldenrodyellow'
    #wav_ax = plt.axes([0.25, 0.1, 0.55, 0.06], facecolor=axcolor)
    allowed_vals = np.arange(0, len(plot_wav))
    
    wav_slider = Slider(ax=context_ax[1], label='Wavelength',
            valmin=0, valmax=len(plot_wav)-1,
            valinit=init, valstep=allowed_vals)

    def f(ixt):
        return context_img[ixt,:,:]
    
    # The function to be called anytime a slider's value changes
    def update(val):
        # Update the context image.
        img_plot.set_data(f(val))
        
        # Update the wavelength position.
        for wav_position in wav_positions:
            wav_position.set_xdata([plot_wav[val].value, plot_wav[val].value])
            
        # Update the title
        context_ax[0].set_title('Stokes ' + meta['stokes'] + '\n $\lambda$ = ' + str(np.round(plot_wav[val],2)))
        fig.canvas.draw_idle()

    # register the update function with the slider
    wav_slider.on_changed(update)

    plt.show()
    
    return wav_slider

def _plot_all_profiles(wavelengths, data, context_img, plot_u, init=0, proj=None, ax=None, meta=None, **kwargs):
    """
    ax: a list with axis instances.
    """
    
    ### Setup the axes in the figure ###
    fig = plt.figure()
    gs0 = fig.add_gridspec(1, 2, left=0.1, right=0.95, width_ratios=[2, 1],
                          hspace=0.2, wspace=0.2)

    gs00 = gs0[0].subgridspec(3, 2, wspace=0.05, hspace=0.05, 
                              width_ratios=[1,1], height_ratios=[4.5,4.5,1.5])
    gs01 = gs0[1].subgridspec(4, 1, hspace=0.0, wspace=0.0)

    # Axes are mapped as: [0,1,2,3] = [I,Q,V,U]
    # [0,0],[0,1],[1,0],[1,1]
    ax = {}
    stokes_list = [['I','Q'],['V','U']]
    for a in range(2):
        for b in range(2):
            ax[stokes_list[a][b]] = fig.add_subplot(gs00[a,b], projection=proj)
            
    slider_ax = fig.add_subplot(gs00[2,:])
    
    l, b, w, h = slider_ax.get_position().bounds
    scale_w = 0.7
    scale_h = 0.7
    slider_ax.set_position([l + w*(1-scale_w)/2, b - w*(1-scale_h)/5, 
                                w*scale_w, h*scale_h])
    
    stokes_ax = []
    for a in range(4):
        stokes_ax.append(fig.add_subplot(gs01[a]))
    
    # Title of the plot.
    plot_wav = wavelengths.to(plot_u)
    title_txt = ax['I'].text(1.0, 1.05, '$\lambda$ = ' + str(np.round(plot_wav[init],2)),
                    transform=ax['I'].transAxes, ha='center')
    
    # Plot the data at the initial wavelength.
    I_img = ax['I'].imshow(context_img[0,init,:,:], origin='lower', aspect='auto')
    Q_img = ax['Q'].imshow(context_img[1,init,:,:], origin='lower', aspect='auto')
    U_img = ax['U'].imshow(context_img[2,init,:,:], origin='lower', aspect='auto')
    V_img = ax['V'].imshow(context_img[3,init,:,:], origin='lower', aspect='auto')
    
    ax['I'].set_title('I', y=0.8, fontweight=700, loc='center')
    ax['Q'].set_title('Q', y=0.8, fontweight=700, loc='center')
    ax['U'].set_title('U', y=0.8, fontweight=700, loc='center')
    ax['V'].set_title('V', y=0.8, fontweight=700, loc='center')
    
    for i in ax:
        ax[i].tick_params(axis='both', direction='in')
        #ax[i].coords[0].set_axislabel('')
        #ax[i].coords[1].set_axislabel('')
            
    ax['I'].coords[0].set_ticklabel_visible(False)
    ax['Q'].coords[0].set_ticklabel_visible(False)
    ax['Q'].coords[1].set_ticklabel_visible(False)
    ax['U'].coords[1].set_ticklabel_visible(False)
    
    #ax['V'].set_xlabel('V')
    ax['V'].coords[0].set_axislabel(' ')
    ax['V'].coords[1].set_axislabel(' ')
    
    ax['U'].coords[0].set_axislabel(' ')
    ax['I'].coords[1].set_axislabel(' ')
    
    ax['I'].text(-0.3, 0.0, 'Helioprojective Latitude',
                 transform=ax['I'].transAxes, 
                 rotation = 'vertical', va='center')
    
    ax['V'].text(1.0, -0.22, 'Helioprojective Longitude',
                 transform=ax['V'].transAxes, 
                 rotation = 'horizontal', ha='center')
    
    # Plot a marker over the coordinates selected.
    for i in ax:
        ax[i].plot(meta['x0_pix'], meta['y0_pix'], '+r')
    
    # Plot the Stokes spectra at the coordinates.
    stokes_ax[0].plot(plot_wav.value, data[0,:], **kwargs)
    stokes_ax[1].plot(plot_wav.value, data[1,:], **kwargs)
    stokes_ax[2].plot(plot_wav.value, data[2,:], **kwargs)
    stokes_ax[3].plot(plot_wav.value, data[3,:], **kwargs)

    stokes_ax[0].set_title('I', y=0.75, fontweight=700, loc='center')
    stokes_ax[1].set_title('Q', y=0.75, fontweight=700, loc='center')
    stokes_ax[2].set_title('U', y=0.75, fontweight=700, loc='center')
    stokes_ax[3].set_title('V', y=0.75, fontweight=700, loc='center')

    # Plot the wavelength position shown in the context image.
    wav_positions = []
    for ax in stokes_ax:
        y0, y1 = ax.get_ylim()
        wav_positions.append(ax.plot([plot_wav[init].value, plot_wav[init].value], [y0,y1], '--k')[0])

    title_string = 'x0 = ' + str(round(meta['x0'].value, 1)) + ' arcsec \n' + \
                'y0 = ' + str(round(meta['y0'].value, 1)) + ' arcsec'
    
    stokes_ax[0].text(0.02, 1.03, title_string, fontsize=10,
                      horizontalalignment='left', verticalalignment='bottom',
                      transform=stokes_ax[0].transAxes)
    
    stokes_ax[3].ticklabel_format(useOffset=False)
    stokes_ax[3].set_xlabel('Wavelength [' + plot_wav[0].unit.to_string() + ']')
    plt.setp(stokes_ax[0].get_xticklabels(), visible=False)
    plt.setp(stokes_ax[1].get_xticklabels(), visible=False)
    plt.setp(stokes_ax[2].get_xticklabels(), visible=False)
    
    for i in stokes_ax:
        i.tick_params(axis='both', direction='in')
        
    # Make a horizontal slider to control the frequency.
    axcolor = 'lightgoldenrodyellow'
    allowed_vals = np.arange(0, len(plot_wav))
    
    wav_slider = Slider(ax=slider_ax, label='Wavelength',
            valmin=0, valmax=len(plot_wav)-1,
            valinit=init, valstep=allowed_vals)
        
    def f(ist, ixt):
        return context_img[ist,ixt,:,:]
    
    # The function to be called anytime a slider's value changes
    def update(val):
        # Update the context image.
        #img_plot.set_data(f(val))
        I_img.set_data(f(0,val))
        Q_img.set_data(f(1,val))
        U_img.set_data(f(2,val))
        V_img.set_data(f(3,val))
        
        # Update the wavelength position.
        for wav_position in wav_positions:
            wav_position.set_xdata([plot_wav[val].value, plot_wav[val].value])
            
        # Update the title
        title_txt.set_text('$\lambda$ = ' + str(np.round(plot_wav[val],2)))
        fig.canvas.draw_idle()

    # register the update function with the slider
    wav_slider.on_changed(update)

    plt.show()
    
    return wav_slider
    #return ax
    
def _plot_all_data(wavelengths, data, plot_u, init=0, proj=None, ax=None, meta=None, **kwargs):
    """
    Draw the default plot showing 2D images of the four Stokes parameters
    along with a way to change the wavelength.
    """
    
    #fig, ax = _subplots(ax, nrows=2, ncols=2, figsize=[5, 5], dpi=120, 
    #                    sharex=False, sharey=False, subplot_kw={'projection':proj})
    #fig.subplots_adjust(bottom=0.1, top=0.9, left=0.15, right=0.9, wspace=0.0, hspace=0.0)
    
    fig = plt.figure(figsize=[5,5], dpi=120)
    ax = fig.subplot_mosaic([["I", "Q"],["V", "U"]], 
                                 sharex=True, sharey=True, 
                                 subplot_kw={'projection':proj})
    fig.subplots_adjust(bottom=0.2, top=0.9, left=0.15, right=0.96, wspace=0.05, hspace=0.05)
    
    plot_wav = wavelengths.to(plot_u)
    
    # Plot the initial wavelength
    title_txt = ax['I'].text(1.0, 1.05, '$\lambda$ = ' + str(np.round(plot_wav[init],2)),
                    transform=ax['I'].transAxes, ha='center')
    
    # Plot the data at the initial wavelength.
    I_img = ax['I'].imshow(data[0,init,:,:], origin='lower', aspect='auto')
    Q_img = ax['Q'].imshow(data[1,init,:,:], origin='lower', aspect='auto')
    U_img = ax['U'].imshow(data[2,init,:,:], origin='lower', aspect='auto')
    V_img = ax['V'].imshow(data[3,init,:,:], origin='lower', aspect='auto')
    
    ax['I'].set_title('I', y=0.84, fontweight=700, loc='center')
    ax['Q'].set_title('Q', y=0.84, fontweight=700, loc='center')
    ax['U'].set_title('U', y=0.84, fontweight=700, loc='center')
    ax['V'].set_title('V', y=0.84, fontweight=700, loc='center')
    
    for i in ax:
        ax[i].tick_params(axis='both', direction='in')
        #ax[i].coords[0].set_axislabel('')
        #ax[i].coords[1].set_axislabel('')
            
    ax['I'].coords[0].set_ticklabel_visible(False)
    ax['Q'].coords[0].set_ticklabel_visible(False)
    ax['Q'].coords[1].set_ticklabel_visible(False)
    ax['U'].coords[1].set_ticklabel_visible(False)
    
    #ax['V'].set_xlabel('V')
    ax['V'].coords[0].set_axislabel(' ')
    ax['V'].coords[1].set_axislabel(' ')
    
    ax['U'].coords[0].set_axislabel(' ')
    ax['I'].coords[1].set_axislabel(' ')
    
    ax['I'].text(-0.34, 0.0, 'Helioprojective Latitude',
                 transform=ax['I'].transAxes, 
                 rotation = 'vertical', va='center')
    
    ax['V'].text(1.0, -0.22, 'Helioprojective Longitude',
                 transform=ax['V'].transAxes, 
                 rotation = 'horizontal', ha='center')
    
    # Make a horizontal slider to control the frequency.
    axcolor = 'lightgoldenrodyellow'
    wav_ax = plt.axes([0.25, 0.05, 0.55, 0.06], facecolor=axcolor)
    allowed_vals = np.arange(0, len(wavelengths))
    
    wav_slider = Slider(ax=wav_ax, label='Wavelength',
            valmin=0, valmax=len(wavelengths)-1,
            valinit=init, valstep=allowed_vals)
        
    def f(ist, ixt):
        return data[ist,ixt,:,:]
    
    # The function to be called anytime a slider's value changes
    def update(val):
        # Update the context image.
        #img_plot.set_data(f(val))
        I_img.set_data(f(0,val))
        Q_img.set_data(f(1,val))
        U_img.set_data(f(2,val))
        V_img.set_data(f(3,val))
        
        # Update the wavelength position.
        #for wav_position in wav_positions:
        #    wav_position.set_xdata([plot_wav[val].value, plot_wav[val].value])
            
        # Update the title
        title_txt.set_text('$\lambda$ = ' + str(np.round(plot_wav[val],2)))
        fig.canvas.draw_idle()

    # register the update function with the slider
    wav_slider.on_changed(update)

    plt.show()
    
    return wav_slider
    
def _plot_image(data, ax=None, proj=None, meta=None, plot_title=None, **kwargs):
    
    if ax is None:
        fig, ax = _subplots(ax, nrows=1, ncols=1, figsize=[5, 5], dpi=120, subplot_kw={'projection':proj})
        fig.subplots_adjust(bottom=0.25, top=0.85, left=0.05, right=0.90, wspace=0.0, hspace=0.0)

    ax.imshow(data, **kwargs)
    
    ax.set_title(plot_title, fontsize=10)
    ax.set_xlabel('Helioprojective Longitude')
    ax.set_ylabel('Helioprojective Latitude')
    
def _plot_3d_cube(wavelengths, data, plot_u, ax=None, proj=None, meta=None, init=0, **kwargs):
    fig, ax = _subplots(ax, nrows=1, ncols=1, figsize=[5, 5], dpi=120, subplot_kw={'projection':proj})
    fig.subplots_adjust(bottom=0.25, top=0.85, left=0.05, right=0.9, wspace=0.0, hspace=0.0)
    
    img_plot = ax.imshow(data[init,:,:], **kwargs)
    plot_wav = round(wavelengths[init].to(plot_u).value,3)
    ax.set_title('Stokes ' + meta['stokes'] + '\n $\lambda$ = ' + str(plot_wav) + ' ' + str(plot_u))
    
    ax.set_xlabel('Helioprojective Longitude')
    ax.set_ylabel('Helioprojective Latitude')
    
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
        plot_wav = round(wavelengths[val].to(plot_u).value,3)
        ax.set_title('Stokes ' + meta['stokes'] + '\n $\lambda$ = ' + str(plot_wav) + ' ' + str(plot_u))
        fig.canvas.draw_idle()

    # register the update function with the slider
    wav_slider.on_changed(update)

    plt.show()
    
    return wav_slider
