0.5.0 (2022-03-24)
==================

This is the first release of StokesPy!

New Features
------------
- Added ``StokesCube`` class for representing maps of Stokes profiles
- Added ``MagVectorCube`` class for representing map of vector magnetic field
- Added plotting functionality for ``StokesCube``, ``MagVectorCube`` and
  the derived slice classes. Plots in a jupyter widget environment contain
  interactive sliders when necessary.
- Added HMI data downloader/loader functions ``load_HMI_stokes()`` and ``load_HMI_magvec()``
- Added Hinode SOT/SP data loader functions ``load_HinodeSP_stokes()`` and ``load_HinodeSP_magvec()``
- Added installation and usage documentation
