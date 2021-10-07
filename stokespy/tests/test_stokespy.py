import pytest
import astropy.units as u
import numpy as np
import stokespy

def generate_data():
    # Order: stokes, wavelength, x, y
    shape = (4, 10, 20, 30)
    data = np.arange(np.product(shape))
    return data.reshape(shape)

def test_init_nowcs():
    data = generate_data()
    stokescube = stokespy.StokesCube(data)
    assert stokescube.data is data
