import numpy as np
import ndcube
import astropy.wcs
import astropy.units as u
import astropy.coordinates
from astropy.coordinates import SkyCoord, SpectralCoord
from astropy.wcs.wcsapi import SlicedLowLevelWCS, HighLevelWCSWrapper

import matplotlib.pyplot as plt

from . import plotting
from matplotlib.widgets import Slider, Button

def make_def_wcs(naxis=3, ctype=None, cunit=None):
    """
    Generate a default WCS object
    
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

class StokesParamCube(ndcube.ndcube.NDCube):
    """Class representing a 2D map of a single Stokes profile with dimensions (wavelength, coord1, coord2)."""
    
    def __init__(self, data, wcs, **kwargs):
        
        # Init base NDCube with data and wcs
        super(StokesParamCube, self).__init__(data, wcs=wcs, **kwargs)
        
        self.n_spectral = None
        self._spectral_axis = None
        
        if self.wcs.pixel_n_dim == 3:
            self.n_spectral = self.data.shape[0]
            self._spectral_axis = self._spectral_slice().array_index_to_world_values(np.arange(self.n_spectral)) * u.Quantity(1, self.wcs.world_axis_units[2])
        
    def _spectral_slice(self):
        """Slice of the WCS containing only the spectral axis"""
        wcs_slice = [0] * self.wcs.pixel_n_dim
        wcs_slice[0] = slice(0, self.n_spectral)
        wcs_slice = SlicedLowLevelWCS(self.wcs.low_level_wcs, wcs_slice)
        return wcs_slice
    
    def get_wav_ind(self,wavelength):
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
    
    def plot(self, plot_u=u.nm, wavelength=None, **kwargs):
        
        if wavelength is None:
            # Default to the first wavelength.
            ix = 0
        else:
            ix = self.get_wav_ind(wavelength) 
        
        wav_slider = plotting._plot_3d_cube(self._spectral_axis, self.data, plot_u, proj=self[ix,:,:].wcs, meta=self.meta, init=ix, origin='lower', **kwargs)
        
        return wav_slider


class StokesParamMap(ndcube.ndcube.NDCube):
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
        "Return StokesProfile over a pixel area."
        pass
    
    def plot(self, plot_u=u.nm, **kwargs):
        """Plot a map of bandpass intensities"""
        
        plot_wav0 = round(self.meta['wav0'].to(plot_u).value,3)
        
        plot_title = ''
        if 'inst' in self.meta.keys():
            plot_title += self.meta['inst'] + ' '
        
        if self.meta['wav1'] is None:
            plot_title += 'Stokes ' + self.meta['stokes'] + '\n $\lambda$ = ' + str(plot_wav0) + ' ' + str(plot_u)
        else:
            plot_wav1 = round(self.meta['wav1'].to(plot_u).value,3)
            plot_title += 'Stokes ' + self.meta['stokes'] + '\n $\lambda\in$ [' + str(plot_wav0) + ', ' + str(plot_wav1) + '] ' + str(plot_u)
        
        plotting._plot_image(self.data, proj=self.wcs, meta=self.meta, origin='lower', plot_title=plot_title, **kwargs)
        
        
class StokesProfile(ndcube.ndcube.NDCube):
    """Class representing a profile of a single Stokes parameter with dimensions (wavelength)
    """   
    def __init__(self, data, wcs, **kwargs):
            
        # Init base NDCube with data and wcs
        super().__init__(data, wcs=wcs, **kwargs)
        
        # Define spectral_axis attribute from WCS
        self.n_spectral = self.data.shape[0]
        self._spectral_axis = self._spectral_slice().array_index_to_world_values(np.arange(self.n_spectral)) * u.Quantity(1, self.wcs.world_axis_units[0])
        #print(self.wcs)
    
    def _spectral_slice(self):
        """Slice of the WCS containing only the spectral axis"""
        # Assume only a single spectral dimension.
        return self.wcs.low_level_wcs
        
    def plot(self, plot_u=u.nm, **kwargs):
        """ Single panel plot showing the dispersed spectrum"""
        plotting._plot_profile(self._spectral_axis, self.data, plot_u, meta=self.meta, **kwargs)
    
    
class StokesCube(ndcube.ndcube.NDCube):
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
        self._frame = None
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
            self._spectral_axis = self._spectral_slice().array_index_to_world_values(np.arange(self.n_spectral)) * u.Quantity(1, self.wcs.world_axis_units[2])
            
            # Define the observer frame if it exists.
            tmp_pixel = np.zeros(self.wcs.pixel_n_dim, dtype='int')
            tmp_world = self.wcs.pixel_to_world(*tmp_pixel)
            if hasattr(tmp_world[0], 'frame'):
                self.meta['frame'] = tmp_world[0].frame
                self._frame = self.meta['frame']
            
    @property
    def stokes_axis(self):
        """The available Stokes parameters"""
        return self._stokes_axis

    @property
    def spectral_axis(self):
        """The spectral axis in physical units"""
        return self._spectral_axis

    @property
    def frame(self):
        "The observed frame for the dat if it exists"
        return self._frame
    
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
        return self[0,0,:,coord2].wcs.array_index_to_world(np.arange(n_coord1))

    def coord2_axis(self, coord1):
        """The physical axis across the second spatial dimension"""
        # TODO: allow coord1 to be None assuming uniform coord2, return 1D array structure        
        n_coord2 = self.data.shape[3]
        return self[0,0,coord1,:].wcs.array_index_to_world(np.arange(n_coord2))

    ##############################
    ####### Stokes Slices ########
    ##############################
    
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
        #cube_meta = {'stokes': self._stokes_axis[stokes_ix]}
        
        cube_meta = self.meta.copy()
        cube_meta['stokes'] = self._stokes_axis[stokes_ix]
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

    ###########################
    ####### Stokes Maps #######
    ###########################
    
    def _stokes_map(self, stokes_ix, wavelength, stop_wavelength=None):
        """Return a 2D NDCube (coord1, coord2) for a given Stokes parameter and wavelength selection"""        
        
        newcube = self._stokes_slice(stokes_ix)
        
        # Starting index.
        ix_0 = self.get_wav_ind(wavelength) 
        wav0 = self._spectral_axis[ix_0]
        print('ix_0, wav0 = ', ix_0, wav0)
        newcube.meta['wav0'] = wav0
        
        # Stopping index.
        ix_1 = None
        wav1 = None
        if stop_wavelength is not None:
            if self.get_wav_ind(stop_wavelength) < ix_0:
                print('Input stop_wavelength ahead of wavelength. Defaulting to the single wavelength map.')
            else:
                ix_1 = self.get_wav_ind(stop_wavelength)
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
    
    def get_wav_ind(self,wavelength):
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
        return self._stokes_map(3, wavelength, stop_wavelength=stop_wavelength)
    
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
    
    ##############################
    ####### Stokes Profile #######
    ##############################
    
    def _stokes_profile(self, stokes_ix, coords):
        """Return a 1D NDCube (wavelength) for a given Stokes parameter and coordinate selection"""
        
        # Tranform input coordinates into a SkyCoord object.
        coords, coords_pix = self.get_spatial_ind(coords)
        
        newcube = self._stokes_slice(stokes_ix)
        
        newcube = newcube[:, coords_pix[0],coords_pix[1]]
        newcube.meta['x0_pix'] = coords_pix[1]
        newcube.meta['y0_pix'] = coords_pix[0]
        newcube.meta['x0'] = self.coord2_axis(0)[coords_pix[1]].Tx
        newcube.meta['y0'] = self.coord1_axis(0)[coords_pix[0]].Ty
        
        return StokesProfile(newcube.data, newcube.wcs, meta=newcube.meta)
    
    def get_spatial_ind(self,coords):
        """Test if a set of coordinates fit inside the dimensions of the 2D images."""
        
        # TODO: allow to specify coords in physical units
        if (isinstance(coords, list) or isinstance(coords, tuple)) and (len(coords) == 2):
            if isinstance(coords[0], u.Quantity) and isinstance(coords[1], u.Quantity):
                coords = SkyCoord(Tx = coords[0], Ty= coords[1], frame=self.frame)
        elif isinstance(coords, SkyCoord):
            if coords.frame.__class__ is astropy.coordinates.builtin_frames.icrs.ICRS:
                """This is the default frame"""
                Tx = coords.ra.to(self.wcs.world_axis_units[0])
                Ty = coords.dec.to(self.wcs.world_axis_units[1])
                coords = SkyCoord(Tx = Tx, Ty = Ty, frame = self.frame)
        else:
            print('Invalid coordinate type.')
            return None
        
        coords_pix = self[0,0,:,:].wcs.world_to_array_index(coords)
        
        return coords, coords_pix
    
    def I_profile(self, coords):
        """Intensity profile at a specific coordinate"""
        return self._stokes_profile(0, coords)
    
    def Q_profile(self, coords):
        """Linear polarization Q profile at a specific coordinate"""
        return self._stokes_profile(1, coords)
    
    def U_profile(self, coords):
        """Linear polarization U profile at a specific coordinate"""
        return self._stokes_profile(2, coords)
    
    def V_profile(self, coords):
        """Circular polarization profile at a specific coordinate"""
        return self._stokes_profile(3, coords)
    
    def P_profile(self, coords):
        """Total polarization P = sqrt(Q**2 + U**2 + V**2) profile at a specific coordinate"""
        Q = self.Q_profile(coords)
        U = self.U_profile(coords)
        V = self.V_profile(coords)
        P = np.sqrt(Q.data**2 + U.data**2 + V.data**2)
        meta = Q.meta
        meta['stokes'] = 'P'
        return StokesProfile(P, Q.wcs, meta=meta)
    
    def L_profile(self, coords):
        """Linear polarization L = sqrt(Q**2 + U**2) profile at a specific coordinate"""        
        Q = self.Q_profile(coords)
        U = self.U_profile(coords)
        L = np.sqrt(Q.data**2 + U.data**2)
        meta = Q.meta
        meta['stokes'] = 'L'
        return StokesProfile(L, Q.wcs, meta=meta)
    
    def theta_profile(self, coords):
        """Linear polarization angle theta = 0.5 arctan(U/Q) profile at a specific coordinate"""
        Q = self.Q_profile(coords)
        U = self.U_profile(coords)
        theta = np.arctan2(U.data, Q.data)
        meta=Q.meta
        meta['stokes'] = 'theta'
        return StokesProfile(np.degrees(theta) * u.degree, Q.wcs)

    ###############################################
    ##### Plotting functionality for the cube #####
    ###############################################
    
    def plot(self, wavelength=None, coords=None, context=None, plot_u=u.nm, **kwargs):
        """Create a four panel plot showing I,Q,U,V maps at a specific wavelength"""
        
        if (coords is not None): 
            # Tranform input coordinates into a SkyCoord object.
            coords, coords_pix = self.get_spatial_ind(coords)
            
            plt_meta = self.meta.copy()
            plt_meta['x0_pix'] = coords_pix[1]
            plt_meta['y0_pix'] = coords_pix[0]
            plt_meta['x0'] = self.coord2_axis(0)[coords_pix[1]].Tx
            plt_meta['y0'] = self.coord1_axis(0)[coords_pix[0]].Ty
            
            if context is None:
                return plotting._plot_all_profiles(self._spectral_axis,
                                            self.data[:,:,coords_pix[0], coords_pix[1]], 
                                            self.data, plot_u, meta=plt_meta, 
                                            proj=self[0,0,:,:].wcs, **kwargs)
            else:
                context_ind = self._stokes_axis.index(context)
                plt_meta['stokes'] = context
                return plotting._plot_context_all_profiles(self._spectral_axis,
                                                        self.data[:,:,coords_pix[0],coords_pix[1]], 
                                                        self.data[context_ind,:,:,:], plot_u,
                                                        proj=self[0,0,:,:].wcs, meta=plt_meta,
                                                        **kwargs)
        elif coords is None:
            # Default is to plot all four Stokes parameters.
            return plotting._plot_all_data(self._spectral_axis, self.data, plot_u, proj=self[0,0,:,:].wcs, meta=self.meta, **kwargs)
    
class MagVectorMap(ndcube.ndcube.NDCube):
    """Class representing a 2D map of one magnetic field component.
    with dimensions (coord1, coord2).
    """

    def __init__(self, data, wcs, **kwargs):
            
        # Init base NDCube with data and wcs
        super().__init__(data, wcs=wcs, **kwargs)
        
    def plot(self, **kwargs):
        """Plot a map of bandpass intensities"""
        # Set title in keyword arguments
        plot_title = 'Magnetic parameter: ' + self.meta['magnetic_param']
        
        plotting._plot_image(self.data, proj=self.wcs, meta=self.meta, plot_title=plot_title, origin='lower', **kwargs)
               

class MagVectorCube(ndcube.ndcube.NDCube):
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

        if self.wcs.pixel_n_dim == 3:
            # Check and define Stokes axis
            #if len(magnetic_params) != self.data.shape[0]:
                #raise Exception(f"Data contains {self.data.shape[0]} magnetic parameters, " + f"but {magnetic_params} parameters ({len(magnetic_params)} were expected") 
            self._magnetic_axis = magnetic_params
        
    @property
    def magnetic_axis(self):
        """The available magnetic parameters"""
        return self._magnetic_axis

    def _magnetic_map(self, magnetic_ix):
        """Return a 2D NDCube (coord1, coord2) for a given magnetic parameter""" 
        
        newcube_data = self.data[magnetic_ix,:,:]
        newcube_wcs = self.wcs[magnetic_ix,:,:]
        new_meta = {'magnetic_param': self._magnetic_axis[magnetic_ix]}
        newmap = MagVectorMap(newcube_data, newcube_wcs, meta=new_meta)
        
        #newcube = ndcube.NDCube(self.data, self.wcs)[magnetic_ix,:,:]
        return newmap
    
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

    def coord1_axis(self, coord2):
        """The physical axis across the first spatial dimension"""
        # TODO: allow coord2 to be None assuming uniform coord1, return 1D array structure
        n_coord1 = self.data.shape[1]
        return self.wcs[0,:,coord2].array_index_to_world(np.arange(n_coord1))

    def coord2_axis(self, coord1):
        """The physical axis across the second spatial dimension"""
        # TODO: allow coord1 to be None assuming uniform coord2, return 1D array structure        
        n_coord2 = self.data.shape[2]
        return self.wcs[0,coord1,:].array_index_to_world(np.arange(n_coord2))
