Stokes Classes
==============

For spectropolarimetric data, the highest dimensionality object is the
`~stokespy.StokesCube`, which contains dimensions of (stokes, wavelength,
coord1, coord2), i.e. a map on the Sun of full Stokes profiles.  In a
Jupyter notebook environment this multidimensional array can be fully
visualized by calling the ``plot()`` method:

.. code-block:: python

  >>> import stokespy
  >>> stokes = stokespy.instload.load_hinode_stokes(...)
  >>> stokes.plot()

TODO: plot example

Here the wavelength slider allows one to scan the full Stokes data at
each wavelength point.

Convenince methods for slicing give lower-dimensionality objects that
can also be visualized, for example, for a map of only Stokes V:

.. code-block:: python

  >>> stokes.V.plot()

TODO: plot example

Here again the wavelength slider allows one to scan the data across
wavelength.

A map at a single wavelength point can also be obtained:

.. code-block:: python

  >>> import astropy.units as u
  >>> stokes.V_map(630.2 * u.nm).plot()

TODO: plot example

Notice how the wavelength was provided in physical units instead of an
index.  The `~stokespy.StokesCube` world coordinate system (WCS) is
used to map physical (world) coordinates to array index, finding the
nearest valid index.  In this way more general analysis code can be
written that supports arrays of various sizes, or sourced from
different instruments.  If an integer is provided, it will be
interpreted as an index on the wavelength axis:

.. code-block:: python

  >>> stokes.V_map(30).plot()

TODO: plot example

In an analagous way, a line profile from any of the Stokes parameters
can be obtained:

.. code-block:: python

  >>> stokes.V_profile(0 * u.arcsec, 0 * u.arcsec).plot()

TODO: plot example

In this example the Stokes V profile from disk center is plotted.
Again, the advantage StokesPy provides is that this line of code will
provide the disk center profile for any valid `~stokespy.StokesCube`
object, regardless of its dimensionality or instrument of origin.

See the `~stokespy.StokesCube` API reference for the full set of
slicing access methods.  Methods are also provided for common
polarimetric calculations, such as the total polarization with
``stokes.P``.
