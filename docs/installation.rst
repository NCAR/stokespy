.. _installation:

===================
Installing StokesPy
===================

StokesPy requires Python 3.7+, sunpy 3.0+, ndcube 2.0+, astropy 4.2+, numpy and matplotlib.

Installing the Stable Version
-----------------------------

There are two options for installing the stable version of SunPy.

The first is via the miniconda distribution using the conda-forge channel.
(The anaconda distribution can also be used but because miniconda is more lightweight we recommend it.)
For more information on installing the miniconda or anaconda distribution, see the `miniconda website`_.

.. code-block:: console

    $ conda config --add channels conda-forge
    $ conda config --set channel_priority strict
    $ conda install stokespy

To update stokespy do:

.. code-block:: console

    $ conda update stokespy

The second option for installing the stable version of StokesPy is via pip.

.. code-block:: console

        $ pip install stokespy

Then to update StokesPy do:

.. code-block:: console

        $ pip install stokespy --upgrade

.. _dev_install:

Installing the Development Version
----------------------------------

The stable version of StokesPy will be reliable.
However, if you value getting the latest updates over reliability, or want to contribute to the development of StokesPy, you will need to install the development version via `stokespy GitHub repository`_.
Let's step through how to do this using miniconda.
For information on installing the miniconda distribution, see the `miniconda website`_.

First, create a conda environment on your local machine to hold the StokesPy bleeding edge version.
Using a new environment allows you to keep your root environment for stable package releases.
Let's call the new conda environment ``stokespy-dev``.
Type the following into a terminal:

.. code-block:: console

    $ conda config --append channels conda-forge
    $ conda create -n stokespy-dev pip

Be sure to activate the environment, i.e. switch into it.
In Linux or MacOS, type:

.. code-block:: console

    $ conda activate stokespy-dev

In Windows, type:

.. code-block:: console

        > conda stokespy-dev

Next clone the StokesPy repo from GitHub to a new directory.
Let's call it stokespy-git.

.. note::

    If you want to develop StokesPy, you should fork the repository and then clone your fork here and not the main StokesPy repository.

.. code-block:: console

    $ git clone https://github.com/NCAR/stokespy.git stokespy-git

To install, change into the new directory and run the install script.

.. code-block:: console

        $ cd stokespy-git
        $ pip install -e .[dev]

Voila!
The StokesPy development version is now installed!
Be sure you get the latest updates by regularly doing:

.. code-block:: console

    $ git pull origin main

.. _miniconda website: https://docs.conda.io/en/latest/miniconda.html
.. _stokespy GitHub repository: https://github.com/NCAR/stokespy
