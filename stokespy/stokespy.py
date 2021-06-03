import numpy as np
import ndcube
import astropy.wcs
import astropy.units as u

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


class StokesParamCube(ndcube.ndcube.NDCubeBase):
    """Class representing a 2D map of a single Stokes profile with dimensions (wavelength, coord1, coord2)."""
    def plot(self, wavelength=None, coord1=None, coord2=None):
        """Plot a slice of the Stokes parameter cube"""
        print(f"TODO: implement {type(self)}.plot()")

class StokesParamMap(ndcube.ndcube.NDCubeBase):
    """Class representing a 2D map of bandpass intensities of a single Stokes parameter 
    with dimensions (coord1, coord2).
    """    
    def plot(self):
        """Plot a map of bandpass intensities"""
        print(f"TODO: implement {type(self)}.plot()")

class StokesProfile(ndcube.ndcube.NDCubeBase):
    """Class representing a profile of a single Stokes parameter with dimensions (wavelength)
    """
    def plot(self):
        """Plot a Stokes profile"""
        print(f"TODO: implement {type(self)}.plot()")

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

        # Check and define Stokes axis
        if len(stokes_params) != self.data.shape[0]:
            raise Exception(f"Data contains {self.data.shape[0]} Stokes parameters, "+
                            f"but {stokes_params} parameters ({len(stokes_params)} were expected")
        self._stokes_axis = stokes_params
        # TODO: stokes index map for N params < 4; use below

        # Define spectral_axis attribute from WCS
        n_spectral = self.data.shape[1]
        self._spectral_axis = self.wcs[0,:,0,0].array_index_to_world(np.arange(n_spectral))

    @property
    def stokes_axis(self):
        """The available Stokes parameters"""
        return self._stokes_axis

    @property
    def spectral_axis(self):
        """The spectral axis in physical units"""
        return self._spectral_axis

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
        newcube = StokesParamCube(self.data, self.wcs)[stokes_ix]
        if stokes_ix != 0:
            if self.normalize is True:
                # Normalize by I
                I = self._get_stokes_slice(0)
                newcube = ndcube.NDCube(newcube.data / I.data, newcube.wcs)
            elif self.normalize:
                # normalize by non-zero float
                # TODO: sanity check input
                newcube = ndcube.NDCube(newcube.data / normalize, newcube.wcs)
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
        return StokesParamCube(P, Q.wcs)
    
    @property
    def L(self):
        """Linear polarization L = sqrt(Q**2 + U**2) a 3D NDCube (wavelength, coord1, coord2)"""
        Q = self.Q
        U = self.U
        L = np.sqrt(Q.data**2 + U.data**2)
        return StokesParamCube(L, Q.wcs)
    
    @property
    def theta(self):
        """Linear polarization angle theta = 0.5 arctan(U/Q) a 3D NDCube (wavelength, coord1, coord2)"""
        Q = self.Q
        U = self.U
        theta = 0.5 * np.arctan2(U.data, Q.data)
        return StokesParamCube(np.degrees(theta) * u.degree, Q.wcs)
    
    def _spectral_slice(self):
        """Slice of the WCS containing only the spectral axis"""
        n_spectral = self.data.shape[1]
        wcs_slice = [0] * self.wcs.naxis
        wcs_slice[1] = slice(0, n_spectral)
        wcs_slice = self.wcs.slice(wcs_slice)
        return wcs_slice
    
    def _stokes_map(self, stokes_ix, wavelength, stop_wavelength=None):
        """Return a 2D NDCube (coord1, coord2) for a given Stokes parameter and wavelength selection"""        
        newcube = self._stokes_slice(stokes_ix)
        spectral_wcs = self._spectral_slice()
        ix = int(spectral_wcs.world_to_array_index_values(wavelength))
        # TODO: understand how the binning is done better...
        newcube = newcube[ix]
        return StokesParamMap(newcube.data, newcube.wcs)
    
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
        return ndcube.NDCube(P, Q.wcs)
    
    def L_map(self, wavelength, stop_wavelength=None):
        """Linear polarization L = sqrt(Q**2 + U**2) as a 2D NDCube (coord1, coord2)"""
        Q = self.Q_map(wavelength, stop_wavelength=stop_wavelength)
        U = self.U_map(wavelength, stop_wavelength=stop_wavelength)
        L = np.sqrt(Q.data**2 + U.data**2)
        return ndcube.NDCube(L, Q.wcs)
    
    def theta_map(self, wavelength, stop_wavelength=None):
        """Linear polarization angle theta = 0.5 arctan(U/Q) as a 2D NDCube (coord1, coord2)"""
        Q = self.Q_map(wavelength, stop_wavelength=stop_wavelength)
        U = self.U_map(wavelength, stop_wavelength=stop_wavelength)
        theta = np.arctan2(U.data, Q.data)
        return ndcube.NDCube(np.degrees(theta) * u.degree, Q.wcs)
    
    def _stokes_profile(self, stokes_ix, coord1, coord2):
        """Return a 1D NDCube (wavelength) for a given Stokes parameter and coordinate selection"""
        # TODO: allow to specify coords in physical units
        newcube = self._stokes_slice(stokes_ix)
        newcube = newcube[:, coord1, coord2]
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
        """Plot all Stokes parameters"""
        print(f"TODO: implement {type(self)}.plot()")

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
