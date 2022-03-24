.. _installation:

===================
Installing StokesPy
===================

StokesPy requires Python 3.7+, sunpy 3.0+, ndcube 2.0+, astropy 4.2+, numpy and matplotlib.

Installing the Stable Version
-----------------------------

StokesPy can be installed using pip:

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
