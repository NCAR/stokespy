Introduction to StokesPy
========================

StokesPy provides several classes to represent different slices of
spectropolarimetric data and its inversion products.  For
spectropolarimetric data, the highest dimensionality object is the
`~stokespy.StokesCube`, which contains dimensions of (stokes, wavelength,
coord1, coord2), i.e. a map on the Sun of full Stokes profiles.
StokesPy data objects uses the `ndcube` project's `ndcube.NDCube` object as a
base class.  In StokesPy, as in ndcube, data and coordinates are
contained in a single object, allowing array-like slicing operations
to modify both the data and coordinate systems simultanously.  This
object can be constructed by providing a properly ordered data array
and World Coordinate System (WCS) as follows:

.. code-block:: python
   
  >>> stokes = StokesCube(data, wcs)

Where ``data`` should have dimensions (stokes, wavelength, coord1,
coord2) and ``wcs`` should describe those dimensions in physical
(world) coordinates.  The provided ``wcs`` should be
`APE-14 compliant
<https://docs.astropy.org/en/stable/wcs/wcsapi.html>`_, which in
practice generally means being a properly constructed
`astropy.wcs.WCS` object or a `gwcs.wcs.WCS`. The interested user
should consult the following references to better understand the
astropy WCS and how it is used by ndcube, StokesPy's base class:

* `NDCube Introduction <https://docs.sunpy.org/projects/ndcube/en/latest/introduction.html>`_
* `WCS APE-14 <https://docs.astropy.org/en/stable/wcs/wcsapi.html>`_   
* `AstroPy WCS <https://docs.astropy.org/en/stable/wcs/index.html>`_
* `gWCS WCS <https://gwcs.readthedocs.io/en/stable>`_  

However, one of the goals of the StokesPy project is to take advantage of the
WCS while insulating the user from the complexities of setting up a
valid WCS.  This is accomplished through custom-built instrument
loader functions. 

.. code-block:: python
  
  >>> stokes = stokespy.instload.load_HMI_stokes()

In this way, the user can immediately begin working with StokesPy
objects in their analysis.

A comprehensive
`example notebook <https://github.com/NCAR/stokespy_notebooks>`_
has been developed demonstrating loading StokesPy objects from HMI and
Hinode SOT/SP data files, using the slicing functions, and plotting
the objects.
