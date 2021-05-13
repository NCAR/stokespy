import numpy as np
import ndcube
import astropy.units as u


def _lookup_axis_index(wcs, name):
    axis_types = wcs.world_axis_object_components # FITS order
    for ix, ax in enumerate(axis_types):
        if ax[0] == name:
            return(wcs.naxis - 1) - ix # convert to array_index order

class StokesCube(ndcube.NDCube):
    def __init__(self, data, wcs=None, spectral_axis=None, stokes_axis=None,
                 stokes_params=('I', 'Q', 'U', 'V'), normalize=False, **kwargs):
        
        # If not specified, attempt to find wavelength axis in WCS
        if spectral_axis is None and wcs is not None:
            spectral_axis = _lookup_axis_index(wcs, 'spectral')
        
        # If not specified, attempt to find stokes axis in WCS
        if stokes_axis is None and wcs is not None:
            stokes_axis = _lookup_axis_index(wcs, 'stokes')
                
        print("XXX spectral_axis:", spectral_axis, "stokes_axis:", stokes_axis)
        # Init base NDCube with data and wcs
        super().__init__(data, wcs=wcs, **kwargs)
        self._spectral_axis_ix = spectral_axis
        self._stokes_axis_ix = stokes_axis
        self.normalize = normalize
        
        # Define spectral_axis attribute from WCS
        if spectral_axis is not None:
            n_spectral = self.data.shape[spectral_axis]
            wcs_slice = [0] * self.wcs.naxis
            wcs_slice[spectral_axis] = slice(0, n_spectral)
            wcs_slice = self.wcs.slice(wcs_slice)
            self.spectral_axis = wcs_slice.array_index_to_world(np.arange(n_spectral))
            
        # Define Stokes parameter attributes
        if stokes_axis is not None:
            self.stokes_axis = stokes_axis
                
    def _full_slice(self):
        full_slice = [slice(stop) for stop in self.dimensions.value.astype(int)]
        return full_slice

    def _stokes_slice(self, stokes_ix, normalize=False):
        cube_slice = self._full_slice()
        cube_slice[self.stokes_axis] = stokes_ix
        newcube = ndcube.NDCube(self.data, self.wcs)[cube_slice]
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
        return self._stokes_slice(0)

    @property
    def Q(self):
        return self._stokes_slice(1)
        
    @property
    def U(self):
        return self._stokes_slice(2)
    
    @property
    def V(self):
        return self._stokes_slice(3)
    
    @property
    def P(self):
        """Total polarization P = sqrt(Q**2 + U**2 + V**2)"""
        Q = self.Q
        U = self.U
        V = self.V
        P = np.sqrt(Q.data**2 + U.data**2 + V.data**2)
        return ndcube.NDCube(P, Q.wcs)
    
    @property
    def L(self):
        """Linear polarization L = sqrt(Q**2 + U**2)"""
        Q = self.Q
        U = self.U
        L = np.sqrt(Q.data**2 + U.data**2)
        return ndcube.NDCube(L, Q.wcs)
    
    @property
    def theta(self):
        """Linear polarization angle theta = 0.5 arctan(U/Q)"""
        Q = self.Q
        U = self.U
        theta = 0.5 * np.arctan2(U.data, Q.data)
        return ndcube.NDCube(np.degrees(theta) * u.degree, Q.wcs)
    
    def _spectral_slice(self):
        spectral_axis = self._spectral_axis_ix
        n_spectral = self.data.shape[spectral_axis]
        wcs_slice = [0] * self.wcs.naxis
        wcs_slice[spectral_axis] = slice(0, n_spectral)
        wcs_slice = self.wcs.slice(wcs_slice)
        return wcs_slice
    
    def _stokes_map(self, stokes_ix, wavelength, stop_wavelength=None):
        newcube = self._stokes_slice(stokes_ix)
        spectral_wcs = self._spectral_slice()
        ix = int(spectral_wcs.world_to_array_index_values(wavelength))
        # TODO: understand how the binning is done better...
        # TODO: assumes wavelength axis is first in newcube
        # Should I just require an axis order from user?
        return newcube[ix]
    
    def I_map(self, wavelength, stop_wavelength=None):
        return self._stokes_map(0, wavelength, stop_wavelength=stop_wavelength)
    
    def Q_map(self, wavelength, stop_wavelength=None):
        return self._stokes_map(1, wavelength, stop_wavelength=stop_wavelength)
    
    def U_map(self, wavelength, stop_wavelength=None):
        return self._stokes_map(2, wavelength, stop_wavelength=stop_wavelength)
    
    def V_map(self, wavelength, stop_wavelength=None):
        return self._stokes_map(2, wavelength, stop_wavelength=stop_wavelength)
    
    def P_map(self, wavelength, stop_wavelength=None):
        """Total polarization P = sqrt(Q**2 + U**2 + V**2)"""
        Q = self.Q_map(wavelength, stop_wavelength=stop_wavelength)
        U = self.U_map(wavelength, stop_wavelength=stop_wavelength)
        V = self.V_map(wavelength, stop_wavelength=stop_wavelength)
        P = np.sqrt(Q.data**2 + U.data**2 + V.data**2)
        return ndcube.NDCube(P, Q.wcs)
    
    def L_map(self, wavelength, stop_wavelength=None):
        """Linear polarization L = sqrt(Q**2 + U**2)"""
        Q = self.Q_map(wavelength, stop_wavelength=stop_wavelength)
        U = self.U_map(wavelength, stop_wavelength=stop_wavelength)
        L = np.sqrt(Q.data**2 + U.data**2)
        return ndcube.NDCube(L, Q.wcs)
    
    def theta_map(self, wavelength, stop_wavelength=None):
        """Linear polarization angle theta = 0.5 arctan(U/Q)"""
        Q = self.Q_map(wavelength, stop_wavelength=stop_wavelength)
        U = self.U_map(wavelength, stop_wavelength=stop_wavelength)
        theta = np.arctan2(U.data, Q.data)
        return ndcube.NDCube(np.degrees(theta) * u.degree, Q.wcs)
    
    def _stokes_profile(self, stokes_ix, coord1, coord2):
        newcube = self._stokes_slice(stokes_ix)
        # TODO: assumes wavelength axis is first in newcube, then coord1, coord2
        # Should I just require an axis order from user?
        return newcube[:, coord1, coord2]
    
    def I_profile(self, coord1, coord2):
        return self._stokes_profile(0, coord1, coord2)
    
    def Q_profile(self, coord1, coord2):
        return self._stokes_profile(1, coord1, coord2)
    
    def U_profile(self, coord1, coord2):
        return self._stokes_profile(2, coord1, coord2)
    
    def V_profile(self, coord1, coord2):
        return self._stokes_profile(3, coord1, coord2)
    
    def P_profile(self, coord1, coord2):
        Q = self.Q_profile(coord1, coord2)
        U = self.U_profile(coord1, coord2)
        V = self.V_profile(coord1, coord2)
        P = np.sqrt(Q.data**2 + U.data**2 + V.data**2)
        return ndcube.NDCube(P, Q.wcs)
    
    def L_profile(self, coord1, coord2):
        Q = self.Q_profile(coord1, coord2)
        U = self.U_profile(coord1, coord2)
        P = np.sqrt(Q.data**2 + U.data**2)
        return ndcube.NDCube(P, Q.wcs)
    
    def theta_profile(self, coord1, coord2):
        """Linear polarization angle theta = 0.5 arctan(U/Q)"""
        Q = self.Q_profile(coord1, coord2)
        U = self.U_profile(coord1, coord2)
        theta = np.arctan2(U.data, Q.data)
        return ndcube.NDCube(np.degrees(theta) * u.degree, Q.wcs)
