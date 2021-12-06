import numpy as np
import ndcube
import astropy.wcs
import astropy.units as u
from astropy.wcs.wcsapi import SlicedLowLevelWCS, HighLevelWCSWrapper

import matplotlib.pyplot as plt

from . import plotting
from matplotlib.widgets import Slider, Button

def make_def_wcs(naxis=3, ctype=None, cunit=None):
    """
    Function that generates a default wcs object.
    
    Parameters
    -----------
    naxis: `int`
        Number of axes that the NDCube will have.
    ctype: `tuple`
        Tuple of strings containing the axes types.
    cunit: `tuple`
        Tuple of strings containing the units for each axes. Must have the same number of elements as ctype.
    """
    wcs = astropy.wcs.WCS(naxis=naxis)
    wcs.wcs.ctype = ctype
    wcs.wcs.cunit = cunit
    wcs.wcs.set()
    return wcs


from . import plotting

class StokesParamCube(ndcube.ndcube.NDCubeBase):
    """Class representing a 2D map of a single Stokes profile with dimensions (wavelength, coord1, coord2)."""
    
    def __init__(self, data, wcs, **kwargs):
        
        # Init base NDCube with data and wcs
        super(StokesParamCube, self).__init__(data, wcs=wcs, **kwargs)
        
        self.n_spectral = None
        self._spectral_axis = None
        #self.stokes = self.meta['stokes']
        
        ## Define spectral_axis attribute from WCS
        # Find the axis index corresponding to the spectral axis
        # This is not necessary if the data structure passed has a known shape
        #ax_types = self.wcs.get_axis_types()
        #for i in range(len(ax_types)):
        #    if ax_types[i]['coordinate_type'] == 'spectral':
        #        iax = len(ax_types) - 1 - i 
        #self.n_spectral = self.data.shape[iax]
        if self.wcs.pixel_n_dim == 3:
            self.n_spectral = self.data.shape[0]
            self._spectral_axis = self._spectral_slice().array_index_to_world_values(np.arange(self.n_spectral))
        #print(self.wcs)
        
    def _spectral_slice(self):
        """Slice of the WCS containing only the spectral axis"""
        wcs_slice = [0] * self.wcs.pixel_n_dim
        wcs_slice[0] = slice(0, self.n_spectral)
        wcs_slice = SlicedLowLevelWCS(self.wcs.low_level_wcs, wcs_slice)
        return wcs_slice
    
    def _get_index(self,wavelength):
        """Test if a wavelength is inside the wavelength axis for the object and return the array index corresponding to that wavelength """
        
        # Test if the value is either an integer index or astropy quantity
        if isinstance(wavelength, int):
            ix = wavelength
            if (ix < 0) or (ix > self.n_spectral-1):
                ix = 0 if ix < 0 else (self.n_spectral-1)
                print('Warning: Wavelength selected outside of axis range: {} {}'.\
                        format(self._spectral_axis[0], self._spectral_axis[-1]))
                print('Defaulting to nearest wavelength at {}'.\
                        format(self._spectral_axis[ix]))
                
        elif isinstance(wavelength,u.Quantity):
            # Unit check the input wavelength
            wav = wavelength.to(self._spectral_slice().world_axis_units[0])
            ix = int(self._spectral_slice().world_to_array_index_values((wav.value,))[0])

            # Check that the wavelength is inside the cube spectral
            # axis range.
            if (ix < 0) or (ix > self.n_spectral-1):
                ix = 0 if ix < 0 else (self.n_spectral-1)
                print('Warning: Wavelength selected outside of axis range: {} {}'.\
                        format(self._spectral_axis[0], self._spectral_axis[-1]))
                print('Defaulting to nearest wavelength at {}'.\
                        format(self._spectral_axis[ix]))
        else:
            print('Warning: the quantity entered must be either an Astropy quantity with a unit of length or an integer index value.')
            pass
                
        return ix
    
    def plot(self, wavelength=None):
        
        if wavelength is None:
            # Default to the first wavelength.
            ix = 0
        else:
            ix = self._get_index(wavelength) 
        
        wav_slider = plotting._plot_3d_cube(self._spectral_axis, self.data, proj=self[ix,:,:].wcs, meta=self.meta, init=ix, origin='lower')
        
        return wav_slider
        
    def plot_old(self, wavelength=None, coord1=None, coord2=None):
        
        if wavelength is None:
            # Default to the first wavelength.
            ix = 0
        else:
            ix = self._get_index(wavelength) 
        
        # Create the plot window.
        image_display, ax = plt.subplots(nrows=1, ncols=1, figsize=[4, 4], dpi=120)
        image_display.subplots_adjust(bottom=0.25, top=0.85, left=0.05, right=0.9, wspace=0.0, hspace=0.0)
        
        img_plot = ax.imshow(self.data[ix,:,:], origin='lower')
        ax.set_title('Stokes ' + self.meta['stokes'] + '\n Single image \n $\lambda$ = ' + str(self._spectral_axis[ix]))
        
        # Make a horizontal slider to control the frequency.
        axcolor = 'lightgoldenrodyellow'
        wav_ax = plt.axes([0.25, 0.1, 0.55, 0.06], facecolor=axcolor)
        allowed_vals = []
        #for i in range(len(self._spectral_axis)):
        #    allowed_vals.append(self._spectral_axis[i])
        #allowed_vals = np.asarray(allowed_vals)
        
        allowed_vals = np.arange(0, self.n_spectral)
        
        #wav_slider = Slider(
        #    ax=wav_ax,
        #    label='Wavelength',
        #    valmin=self._spectral_axis[0],
        #    valmax=self._spectral_axis[-1],
        #    valinit=self._spectral_axis[ix],
        #    valstep=allowed_vals
        #)
        wav_slider = Slider(ax=wav_ax, label='Wavelength',
            valmin=0, valmax=self.n_spectral-1,
            valinit=ix,
            valstep=allowed_vals
        )
        
        def f(wav):
            #ixt=int(self._spectral_slice().world_to_array_index_values(wav))
            #ixt = self._get_index(wav)
            ixt = wav
            return self.data[ixt,:,:]
            
        # The function to be called anytime a slider's value changes
        def update(val):
            img_plot.set_data(f(val))
            ax.set_title('Stokes ' + self.meta['stokes'] + '\n Single image \n $\lambda$ = ' + str(self._spectral_axis[val]))
            #print(val)
            #ax.set_title('blah' + str(val))
            image_display.canvas.draw_idle()
        
        # register the update function with the slider
        wav_slider.on_changed(update)
        
        """Plot a slice of the Stokes parameter cube"""
        #print(f"TODO: implement {type(self)}.plot()")
        
        plt.show()

        return wav_slider
        
class StokesParamMap(ndcube.ndcube.NDCubeBase):
    """Class representing a 2D map of bandpass intensities of a single Stokes parameter 
    with dimensions (coord1, coord2).
    """

    def __init__(self, data, wcs, **kwargs):
            
        # Init base NDCube with data and wcs
        super().__init__(data, wcs=wcs, **kwargs)
    
        self.wav0 = None
        self.wav1 = None
    
        if 'wav0' in self.meta:
            self.wav0 = self.meta['wav0']
        if 'wav1' in self.meta:
            self.wav1 = self.meta['wav1']
        if 'stokes' in self.meta:
            self.stokes = self.meta['stokes']
    
    def average(self):
        "Return StokesProfile"
        pass
    
    def plot(self):
        """Plot a map of bandpass intensities"""
        
        plotting._plot_image(self.data, proj=self.wcs, meta=self.meta, origin='lower')
        
        #image_display, ax = plt.subplots(nrows=1, ncols=1, figsize=[5, 5], dpi=120, subplot_kw={'projection': self.wcs})
        #image_display.subplots_adjust(bottom=0.25, top=0.85, left=0.05, right=0.90, wspace=0.0, hspace=0.0)
        
        #img_plot = ax.imshow(self.data, origin='lower')
        
        #if self.wav1 is None:
        #    ax.set_title('Stokes ' + self.meta['stokes'] + '\n Single image \n $\lambda$ = ' + str(self.wav0))
        #else:
        #    ax.set_title('Stokes ' + self.meta['stokes'] + '\n Summed over \n $\lambda\in$ [' + str(self.wav0) + ' ' + str(self.wav1) + ']')
        #print(f"TODO: implement {type(self)}.plot()")

        
class StokesProfile(ndcube.ndcube.NDCubeBase):
    """Class representing a profile of a single Stokes parameter with dimensions (wavelength)
    """   
    def __init__(self, data, wcs, **kwargs):
            
        # Init base NDCube with data and wcs
        super().__init__(data, wcs=wcs, **kwargs)
        
        # Define spectral_axis attribute from WCS
        self.n_spectral = self.data.shape[0]
        self._spectral_axis = self._spectral_slice().array_index_to_world_values(np.arange(self.n_spectral))
        #print(self.wcs)
    
    def _spectral_slice(self):
        """Slice of the WCS containing only the spectral axis"""
        # Assume only a single spectral dimension.
        return self.wcs.low_level_wcs
        
    def plot(self):
    
        plotting._plot_profile(self._spectral_axis,self.data)
    
        #print(self.n_spectral)
        #print(self._spectral_axis)
        #n_spec = self.data.size
        
        # Create the plot window.
        #image_display, ax = plt.subplots(nrows=1, ncols=1, figsize=[4, 4], dpi=100)
        #image_display.subplots_adjust(bottom=0.1, top=0.9, left=0.05, right=0.90, wspace=0.0, hspace=0.0)
        
        #ax.plot(self._spectral_axis,self.data)
        
        # The object does not know which Stokes parameter it 
        # represents.
        #ax.set_title('Stokes X')
        #ax.set_xlabel('Wavelength')
        
        """Plot a Stokes profile"""
        #print(f"TODO: implement {type(self)}.plot()")
        
class StokesCube(ndcube.ndcube.NDCubeBase):
    """
    Class representing a 2D map of Stokes profiles with dimensions (stokes, wavelength, coord1, coord2).

    Parameters
    ----------
    data: `numpy.ndarray`
        The array holding the actual data in this object.  The array index order must be 
        (stokes, wavelength, coord1, coord2).

    wcs: `astropy.wcs.wcsapi.BaseLowLevelWCS`, `astropy.wcs.wcsapi.BaseHighLevelWCS`, optional
        The WCS object containing the axes' information.  If not provided, a WCS is constructed 
        using `wavelength_unit` and `coordinate_unit`, which default to pixels.

    stokes_params: `tuple` of `str`
        Tuple containing all or part of ('I', 'Q', 'U', 'V') defining the number and kind of 
        Stokes parameters available.

    normalize: `bool`
        Normalization for polarization parameters Q, U, V.  If `True`, then polarization parameters 
        are normalized by the intensity at each wavelength, e.g. Q => Q/I.  If a non-zero scalar is 
        provided then that will be used as the normalization, e.g. for a chosen continuum intensity.

    Additional kwargs are passed to the `NDCube` base class.
    """

    def __init__(self, data, wcs=None, stokes_params=('I', 'Q', 'U', 'V'), normalize=False, **kwargs):
        if wcs is None:
            # Define a default WCS where coordinates and wavelength axis are
            # in pixel units.  Note: cannot use "WAVE" ctype;
            # astropy.wcs.WCS enforces length units for that name
            wcs = make_def_wcs(naxis=4, ctype=["COORD2", "COORD1", "WAVEIX", "STOKES"], 
                               cunit=['pix', 'pix', 'pix', ''])

        # Init base NDCube with data and wcs
        super().__init__(data, wcs=wcs, **kwargs)
        self.normalize = normalize
        self.n_spectral = None
        self._spectral_axis = None
        self._stokes_axis = None
        
        if self.wcs.pixel_n_dim == 4:
            # Check and define Stokes axis.
            if len(stokes_params) != self.data.shape[0]:
                raise Exception(f"Data contains {self.data.shape[0]} Stokes parameters, "+
                            f"but {len(stokes_params)} parameters  were expected: {stokes_params}")
            self._stokes_axis = stokes_params
            # TODO: stokes index map for N params < 4; use below

            # Define spectral_axis attribute from WCS
            self.n_spectral = self.data.shape[1]
            self._spectral_axis = self._spectral_slice().array_index_to_world_values(np.arange(self.n_spectral))
        
    @property
    def stokes_axis(self):
        """The available Stokes parameters"""
        return self._stokes_axis

    @property
    def spectral_axis(self):
        """The spectral axis in physical units"""
        return self._spectral_axis

    def _spectral_slice(self):
        """Slice of the WCS containing only the spectral axis"""
        wcs_slice = [0] * self.wcs.pixel_n_dim
        wcs_slice[1] = slice(0, self.n_spectral)
        wcs_slice = SlicedLowLevelWCS(self.wcs.low_level_wcs, wcs_slice)
        return wcs_slice
    
    def coord1_axis(self, coord2):
        """The physical axis across the first spatial dimension"""
        # TODO: allow coord2 to be None assuming uniform coord1, return 1D array structure
        n_coord1 = self.data.shape[2]
        return self.wcs[0,0,:,coord2].array_index_to_world(np.arange(n_coord1))

    def coord2_axis(self, coord1):
        """The physical axis across the second spatial dimension"""
        # TODO: allow coord1 to be None assuming uniform coord2, return 1D array structure        
        n_coord2 = self.data.shape[3]
        return self.wcs[0,0,coord1,:].array_index_to_world(np.arange(n_coord2))

    def _stokes_slice(self, stokes_ix, normalize=False):
        """Return a 3D NDCube (wavelength, coord1, coord2) for a given Stokes parameter"""
        
        # Construct the WCS object for the smaller 3D cube.
        # This function should only called if the  
        d_sh = self.data.shape
        wcs_slice = [0] * self.wcs.pixel_n_dim
        wcs_slice[0] = stokes_ix
        wcs_slice[1] = slice(0, d_sh[1])
        wcs_slice[2] = slice(0, d_sh[2])
        wcs_slice[3] = slice(0, d_sh[3])
        #print(wcs_slice)
        wcs_slice = SlicedLowLevelWCS(self.wcs.low_level_wcs, wcs_slice)
        #newcube = StokesParamCube(self.data[stokes_ix,:,:,:], HighLevelWCSWrapper(wcs_slice), self._stokes_axis[stokes_ix])
        cube_meta = {'stokes': self._stokes_axis[stokes_ix]}
        newcube = StokesParamCube(self.data[stokes_ix,:,:,:], HighLevelWCSWrapper(wcs_slice), meta=cube_meta)
        
        # Normalize the spectra if wanted.
        if stokes_ix != 0:
            if self.normalize is True:
                # Normalize by I
                I = self._stokes_slice(0)
                newcube = StokesParamCube(newcube.data / I.data, newcube.wcs, self._stokes_axis[stokes_ix], meta=newcube.meta)
            elif self.normalize:
                # normalize by non-zero float
                # TODO: sanity check input
                newcube = StokesParamCube(newcube.data / self.normalize, newcube.wcs, self._stokes_axis[stokes_ix], meta=newcube.meta)
        
        #newcube.meta = {'stokes': self._stokes_axis[stokes_ix]}
        
        return newcube

    @property
    def I(self):
        """Intensity as a 3D NDCube (wavelength, coord1, coord2)"""
        return self._stokes_slice(0)

    @property
    def Q(self):
        """Linear polarization Q as a 3D NDCube (wavelength, coord1, coord2)"""
        return self._stokes_slice(1)
        
    @property
    def U(self):
        """Linear polarization U as a 3D NDCube (wavelength, coord1, coord2)"""
        return self._stokes_slice(2)
    
    @property
    def V(self):
        """Circular polarization as a 3D NDCube (wavelength, coord1, coord2)"""
        return self._stokes_slice(3)
    
    @property
    def P(self):
        """Total polarization P = sqrt(Q**2 + U**2 + V**2) a 3D NDCube (wavelength, coord1, coord2)"""
        Q = self.Q
        U = self.U
        V = self.V
        P = np.sqrt(Q.data**2 + U.data**2 + V.data**2)
        return StokesParamCube(P, Q.wcs, meta={'stokes': 'P'})
    
    @property
    def L(self):
        """Linear polarization L = sqrt(Q**2 + U**2) a 3D NDCube (wavelength, coord1, coord2)"""
        Q = self.Q
        U = self.U
        L = np.sqrt(Q.data**2 + U.data**2)
        return StokesParamCube(L, Q.wcs, meta={'stokes': 'L'})
    
    @property
    def theta(self):
        """Linear polarization angle theta = 0.5 arctan(U/Q) a 3D NDCube (wavelength, coord1, coord2)"""
        Q = self.Q
        U = self.U
        theta = 0.5 * np.arctan2(U.data, Q.data)
        return StokesParamCube(np.degrees(theta) * u.degree, Q.wcs, meta={'stokes': 'theta'})
    
    def _stokes_map(self, stokes_ix, wavelength, stop_wavelength=None):
        """Return a 2D NDCube (coord1, coord2) for a given Stokes parameter and wavelength selection"""        
        
        newcube = self._stokes_slice(stokes_ix)
        
        # Starting index.
        ix_0 = self._get_index(wavelength) 
        wav0 = self._spectral_axis[ix_0]
        print('ix_0, wav0 = ', ix_0, wav0)
        newcube.meta['wav0'] = wav0
        
        # Stopping index.
        ix_1 = None
        wav1 = None
        if stop_wavelength is not None:
            if self._get_index(stop_wavelength) < ix_0:
                print('Input stop_wavelength ahead of wavelength. Defaulting to the single wavelength map.')
            else:
                ix_1 = self._get_index(stop_wavelength)
                wav1 = self._spectral_axis[ix_1]
            print('ix_1, wav1 = ', ix_1, wav1)
        newcube.meta['wav1'] = wav1
        
        # Select the data to be included.
        if ix_1 is None:
            newcube_data = newcube.data[ix_0,:,:]
        else:
            # Sum over the selected wavelengths.
            newcube_data = np.sum(newcube.data[ix_0:ix_1+1,:,:], axis=0)
        
        newcube_wcs = newcube[ix_0,:,:].wcs
        newmap = StokesParamMap(newcube_data, newcube_wcs, meta=newcube.meta)
        return newmap
    
    def _get_index(self,wavelength):
        """Test if a wavelength is inside the wavelength axis for the object and return the array index corresponding to that wavelength """
        
        # Test if the value is either an integer index or astropy quantity
        if isinstance(wavelength, int):
            ix = wavelength
            if (ix < 0) or (ix > self.n_spectral-1):
                ix = 0 if ix < 0 else (self.n_spectral-1)
                print('Warning: Wavelength selected outside of axis range: {} {}'.\
                        format(self._spectral_axis[0], self._spectral_axis[-1]))
                print('Defaulting to nearest wavelength at {}'.\
                        format(self._spectral_axis[ix]))
                
        elif isinstance(wavelength,u.Quantity):
            # Unit check the input wavelength
            wav = wavelength.to(self._spectral_slice().world_axis_units[0])
            ix = int(self._spectral_slice().world_to_array_index_values((wav.value,))[0])

            # Check that the wavelength is inside the cube spectral
            # axis range.
            if (ix < 0) or (ix > self.n_spectral-1):
                ix = 0 if ix < 0 else (self.n_spectral-1)
                print('Warning: Wavelength selected outside of axis range: {} {}'.\
                        format(self._spectral_axis[0], self._spectral_axis[-1]))
                print('Defaulting to nearest wavelength at {}'.\
                        format(self._spectral_axis[ix]))
        else:
            print('Warning: the quantity entered must be either an Astropy quantity with a unit of length or an integer index value.')
            pass
                
        return ix
        
    def I_map(self, wavelength, stop_wavelength=None):
        """Intensity as a 2D NDCube (coord1, coord2)"""
        return self._stokes_map(0, wavelength, stop_wavelength=stop_wavelength)
    
    def Q_map(self, wavelength, stop_wavelength=None):
        """Linear polarization Q as a 2D NDCube (coord1, coord2)"""
        return self._stokes_map(1, wavelength, stop_wavelength=stop_wavelength)
    
    def U_map(self, wavelength, stop_wavelength=None):
        """Linear polarization U as a 2D NDCube (coord1, coord2)"""        
        return self._stokes_map(2, wavelength, stop_wavelength=stop_wavelength)
    
    def V_map(self, wavelength, stop_wavelength=None):
        """Circular polarization as a 2D NDCube (coord1, coord2)"""        
        return self._stokes_map(2, wavelength, stop_wavelength=stop_wavelength)
    
    def P_map(self, wavelength, stop_wavelength=None):
        """Total polarization P = sqrt(Q**2 + U**2 + V**2) as a 2D NDCube (coord1, coord2)"""
        Q = self.Q_map(wavelength, stop_wavelength=stop_wavelength)
        U = self.U_map(wavelength, stop_wavelength=stop_wavelength)
        V = self.V_map(wavelength, stop_wavelength=stop_wavelength)
        P = np.sqrt(Q.data**2 + U.data**2 + V.data**2)
        meta = Q.meta
        meta['stokes'] = 'P'
        return StokesParamMap(P, Q.wcs, meta=meta)
    
    def L_map(self, wavelength, stop_wavelength=None):
        """Linear polarization L = sqrt(Q**2 + U**2) as a 2D NDCube (coord1, coord2)"""
        Q = self.Q_map(wavelength, stop_wavelength=stop_wavelength)
        U = self.U_map(wavelength, stop_wavelength=stop_wavelength)
        L = np.sqrt(Q.data**2 + U.data**2)
        meta = Q.meta
        meta['stokes'] = 'L'
        return StokesParamMap(L, Q.wcs, meta=meta)
    
    def theta_map(self, wavelength, stop_wavelength=None):
        """Linear polarization angle theta = 0.5 arctan(U/Q) as a 2D NDCube (coord1, coord2)"""
        Q = self.Q_map(wavelength, stop_wavelength=stop_wavelength)
        U = self.U_map(wavelength, stop_wavelength=stop_wavelength)
        theta = np.arctan2(U.data, Q.data)
        meta = Q.meta
        meta['stokes'] = 'theta'
        return StokesParamMap(np.degrees(theta) * u.degree, Q.wcs, meta=meta)
    
    def _stokes_profile(self, stokes_ix, coord1, coord2):
        """Return a 1D NDCube (wavelength) for a given Stokes parameter and coordinate selection"""
        # TODO: allow to specify coords in physical units
        newcube = self._stokes_slice(stokes_ix)
        #newcube = self._stokes_slice(stokes_ix)[:,coord1,coord2]
        newcube = newcube[:,coord1,coord2] 
        return StokesProfile(newcube.data, newcube.wcs)
        
    def I_profile(self, coord1, coord2):
        """Intensity profile at a specific coordinate"""
        return self._stokes_profile(0, coord1, coord2)
    
    def Q_profile(self, coord1, coord2):
        """Linear polarization Q profile at a specific coordinate"""
        return self._stokes_profile(1, coord1, coord2)
    
    def U_profile(self, coord1, coord2):
        """Linear polarization U profile at a specific coordinate"""
        return self._stokes_profile(2, coord1, coord2)
    
    def V_profile(self, coord1, coord2):
        """Circular polarization profile at a specific coordinate"""
        return self._stokes_profile(3, coord1, coord2)
    
    def P_profile(self, coord1, coord2):
        """Total polarization P = sqrt(Q**2 + U**2 + V**2) profile at a specific coordinate"""
        Q = self.Q_profile(coord1, coord2)
        U = self.U_profile(coord1, coord2)
        V = self.V_profile(coord1, coord2)
        P = np.sqrt(Q.data**2 + U.data**2 + V.data**2)
        return StokesProfile(P, Q.wcs)
    
    def L_profile(self, coord1, coord2):
        """Linear polarization L = sqrt(Q**2 + U**2) profile at a specific coordinate"""        
        Q = self.Q_profile(coord1, coord2)
        U = self.U_profile(coord1, coord2)
        P = np.sqrt(Q.data**2 + U.data**2)
        return StokesProfile(P, Q.wcs)
    
    def theta_profile(self, coord1, coord2):
        """Linear polarization angle theta = 0.5 arctan(U/Q) profile at a specific coordinate"""
        Q = self.Q_profile(coord1, coord2)
        U = self.U_profile(coord1, coord2)
        theta = np.arctan2(U.data, Q.data)
        return StokesProfile(np.degrees(theta) * u.degree, Q.wcs)

    def plot(self, wavelength=None, coord1=None, coord2=None):
        """Create a four panel plot showing I,Q,U,V maps at a specific wavelength"""

        if (coord1 and coord2): 
            # Plot a four panel plot showing I,Q,U,V wavelengths at point(coord1, coord2)
        
             # Create the plot window.
            image_display, ax = plt.subplots(nrows=4, ncols=1, figsize=[2, 5], dpi=100)
            image_display.subplots_adjust(bottom=0.1, top=0.9, left=0.05, right=0.90, wspace=0.0, hspace=0.0)

            ax[0].plot(self.data[0,:,coord1,coord2])
            ax[1].plot(self.data[1,:,coord1,coord2])
            ax[2].plot(self.data[2,:,coord1,coord2])
            ax[3].plot(self.data[3,:,coord1,coord2])
            
            ax[0].set_title('I')
            ax[1].set_title('Q')
            ax[2].set_title('U')
            ax[3].set_title('V')
            
        elif (coord1 is None) and (coord2 is None):
            # Choose a default wavelength if none are provided.
            if wavelength is None:
                # Need a better way to select which wavelength to show.
                ix = int(self._spectral_axis.shape[0]/2)
            elif isinstance(wavelength,astropy.units.Quantity):
                ix = self._get_index(wavelength)
                # Test if the wavelength provided is a Quantity object.
                nwav = len(self._spectral_axis)
                ix = int(self._spectral_slice().world_to_array_index_values(wavelength))
                # Check that the selected value falls within the wavelength array.
                if (ix < 0) or (ix > nwav-1):
                    ix = 0 if ix < 0 else (nwav-1)
                    print('Warning: Wavelength selected outside of range: {} {}'.\
                          format(self._spectral_axis[0], self._spectral_axis[-1]))
                    print('Defaulting to nearest wavelength at {}'.\
                          format(self._spectral_axis[ix]))
            elif isinstance(wavelength, int):
                # Wavelength is array index.
                ix = wavelength
                    
            # Create the plot window.
            image_display, ax = plt.subplots(nrows=2, ncols=2, figsize=[8, 8], dpi=120)
            image_display.subplots_adjust(bottom=0.1, top=0.9, left=0.05, right=0.90, wspace=0.0, hspace=0.0)

            ax[0,0].imshow(self.data[0,ix,:,:], origin='lower', aspect='auto')
            ax[0,1].imshow(self.data[1,ix,:,:], origin='lower', aspect='auto')
            ax[1,0].imshow(self.data[2,ix,:,:], origin='lower', aspect='auto')
            ax[1,1].imshow(self.data[3,ix,:,:], origin='lower', aspect='auto')

            ax[0,0].set_title('I')
            ax[0,1].set_title('Q')
            ax[1,0].set_title('U')
            ax[1,1].set_title('V')
        
        """Plot all Stokes parameters"""
        print(f"TODO1: implement {type(self)}.plot()")
        
class MagVectorCube(ndcube.ndcube.NDCubeBase):
    """
    Class representing a 2D map of inverted magnetic field vectors.
    
    Parameters
    ----------
    data: `numpy.ndarray`
        The array holding the magnetic field data stored in the object. The array index order must be
        (magnetic, coord1, coord2) or should this be (magnetic, coord2, coord1)?
        
    wcs: `astropy.wcs.wcsapi.BaseLowLevelWCS`, `astropy.wcs.wcsapi.BaseHighLevelWCS`, optional
        The WCS object containing the axes' information.  If not provided, a WCS is constructed 
        using `wavelength_unit` and `coordinate_unit`, which default to pixels.
    
    magnetic_params: `tuple` or `str`
        Tuple containing all or part of the magnetic field components ('B', 'inclination', 'azimuth')
    
    """
    
    def __init__(self, data, wcs=None, magnetic_params=('B', 'inclination', 'azimuth'), **kwargs):
        if wcs is None:
            # Define a default WCS where coordinates are defined in pixel units.  
            wcs = make_def_wcs(naxis=3, ctype=["COORD2", "COORD1", "Parameter"],
                               cunit=['pix', 'pix', ''])

        # Init base NDCube with data and wcs
        super().__init__(data, wcs=wcs, **kwargs)

        # Check and define Stokes axis
        if len(magnetic_params) != self.data.shape[0]:
            raise Exception(f"Data contains {self.data.shape[0]} magnetic parameters, "+
                            f"but {magnetic_params} parameters ({len(magnetic_params)} were expected")
        self._magnetic_axis = magnetic_params
        
    @property
    def magnetic_axis(self):
        """The available magnetic parameters"""
        return self._magnetic_axis

    def _magnetic_map(self, magnetic_ix):
        """Return a 2D NDCube (coord1, coord2) for a given magnetic parameter"""        
        newcube = ndcube.NDCube(self.data, self.wcs)[magnetic_ix,:,:]
        return newcube
    
    @property
    def B(self):
        """Magnetic field strength as a 2D NDcube (coord1, coord2)"""
        return self._magnetic_map(0)

    @property
    def inclination(self):
        """Magnetic inclination as a 2D NDCube (coord1, coord2)"""
        return self._magnetic_map(1)
        
    @property
    def azimuth(self):
        """Magnetic azimuth as 2D NDCube (coord1, coord2)"""
        return self._magnetic_map(2)
