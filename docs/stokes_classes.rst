Stokes Classes
==============

For spectropolarimetric data, the highest dimensionality object is the
`~stokespy.StokesCube`, which contains dimensions of (stokes, wavelength,
coord1, coord2), i.e. a map on the Sun of full Stokes profiles.  In a
Jupyter notebook environment this multidimensional array can be fully
visualized by calling the ``plot()`` method:

.. code-block:: python

  >>> import stokespy
  >>> stokes = stokespy.instload.load_HinodeSP_stokes(...)
  >>> stokes.plot()

.. image:: images/StokesCube_plot.png
  :alt: Plot of a StokesCube object

Here the wavelength slider allows one to scan the full Stokes data at
each wavelength point.

Convenince methods for slicing give lower-dimensionality objects that
can also be visualized, for example, for a map of only Stokes V:

.. code-block:: python

  >>> stokes.V.plot()

.. image:: images/StokesParamCube_plot.png
  :alt: Plot of a StokesParamCube object


Here again the wavelength slider allows one to scan the data across
wavelength.

A map at a single wavelength point can also be obtained:

.. code-block:: python

  >>> import astropy.units as u
  >>> stokes.V_map(630.142 * u.nm).plot()

.. image:: images/StokesParamMap_plot1.png
  :alt: Plot of a StokesParamMap object from a physical wavelength

Notice how the wavelength was provided in physical units instead of an
index.  The `~stokespy.StokesCube` world coordinate system (WCS) is
used to map physical (world) coordinates to array index, finding the
nearest valid index.  In this way more general analysis code can be
written that supports arrays of various sizes, or sourced from
different instruments.  If an integer is provided, it will be
interpreted as an index on the wavelength axis:

.. code-block:: python

  >>> stokes.V_map(25).plot()

.. image:: images/StokesParamMap_plot2.png
  :alt: Plot of a StokesParamMap object from a wavelength index


In an analagous way, a line profile from any of the Stokes parameters
can be obtained:

.. code-block:: python

  >>> coord = SkyCoord(Tx = 148 * u.arcsec, Ty = 67 * u.arcsec, frame = stokes.meta['frame'])
  >>> stokes.V_profile(coord).plot()

.. image:: images/StokesProfile_plot.png
  :alt: Plot of a StokesProfile

In this example the Stokes V profile from disk center is plotted.
Again, the advantage StokesPy provides is that this line of code will
provide the disk center profile for any valid `~stokespy.StokesCube`
object, regardless of its dimensionality or instrument of origin.

See the `~stokespy.StokesCube` API reference for the full set of
slicing access methods.  Methods are also provided for common
polarimetric calculations, such as the total polarization with
``stokes.P``.
