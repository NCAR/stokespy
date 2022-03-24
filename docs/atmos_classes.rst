Atmosphere Classes
==================

StokesPy provides classes to represent inversions of
spectropolarimetric data in a set of "atmosphere" classes.  At the
moment only `~stokespy.MagVectorCube`, representing a map of the
vector magnetic field, has been developed.  In the future classes
representing the 3D atmosphere (including depth) will be developed.

.. code-block:: python

  >>> import stokespy
  >>> field = stokespy.instload.load_HinodeSP_magvec(...)
  >>> field.B.plot()

.. image:: images/MagVectorMap_plot.png
  :alt: Plot of a MagVectorMap object
