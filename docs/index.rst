StokesPy Documentation
----------------------

StokesPy is a (future) Sunpy-affiliated packaged for the analysis of
polarized light from the Sun and the inferred solar atmosphere and
magnetic field (inversion products).  StokesPy aims to simplify
importing and visualizing common datasets of polarized light and
inversion products, as well as provide common operations general to
all such datasets.

For polarized spectra, the highest dimensionality data object is the
`StokesCube`, which represents a map (two spatial dimensions) of
profiles (wavelength) of the 4 Stokes parameters (IQUV), making a
seven-dimensional object.  StokesPy is based on the ndcube package
which provides support for coordinate-aware slicing of this
multidimensional data.  StokesPy implements dimensionality-aware
subclasses to represent and visualize maps of a given Stokes parameter
and wavelength, slits, and profiles at a given pixel coordinate.

For inversion products, StokesPy currently implements an object to represent
vector magnetograms, MagVectorCube.  Future versions of StokesPy will provide
representations of slab atmospheres and stratified atmospheres, and
provide tools to aid in the comparison of synthesized polarized
spectra from model atmospheres to observational data.

.. toctree::
   :maxdepth: 1

   introduction
   installation
   stokes_classes
   atmos_classes
   api
   whatsnew/index

.. _getting_help:

Getting Help
============

If you are unsure of the general features of the ndcube package, your
first stop should be this document which we strive to make as
comprehensive, yet understandable, as possible.

If you would like to know more about how to use a specific function,
class, etc., documentation for each object is in the :ref:`api` of
this guide.  If you cannot find the answer to your issue there, would
like to provide feedback, or would like help contributing code to
StokesPy, you can get help directly by mailing XXX.


Reporting Bugs
--------------

If you run into unexpected behavior or a bug please report it.  All
bugs are raised and stored on our `issue tracker`_.  If you are not
sure if your problem is a bug, you can either send a message to the
chat room or open an issue.  Ideally, we would like a short code
example so we can run into the bug on our own machines.

.. _issue tracker: https://github.com/NCAR/stokespy/issues


Using StokesPy in your project
------------------------------

The vision of StokesPy is to provide a common interface and toolset
for the various spectropolarimetric datasets and inversion codes in
the community.  In addition to the data containers, it provides a
repository of instrument/inversion code-specific data loading
routines.  If you are interested in adding data loaders for your
dataset/inversion code to StokesPy, please get in touch with the
developers to begin a discussion.

.. _contributing:

Contributing to StokesPy
========================

We are always enthusiastic to welcome new users and developers who
want to enhance StokesPy You can contribute in several ways, from
providing feedback, reporting bugs, contributing code, and reviewing
pull requests.  There is a role for almost any level of engagement.

If you would like to contribute code, it is strongly recommended that
you first discuss your aims with the StokesPy community.  We strive to
be an open and welcoming community for developers of all experience
levels.  Discussing your ideas before you start can give you new
insights that will make your development easier, lead to a better end
product, and reduce the chances of your work being regretfully
rejected because of an issue you weren't aware of, e.g. the
functionality already exists elsewhere.

The best guide we provide comes from the `sunpy newcomers guide
<https://docs.sunpy.org/en/latest/dev_guide/contents/newcomers.html>`__.
