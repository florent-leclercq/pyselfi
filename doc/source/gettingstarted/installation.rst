Installation
============

This page contains instructions for installing the **pySELFI** code. Two python packages are installed: :mod:`pyselfi` and :mod:`pyselfi_examples`.

Minimal installation using pip
------------------------------

pySELFI requires Python 3.8 or greater. It can be installed easily from the command line using `pip <http://www.pip-installer.org/>`_:

.. code-block:: console

    pip install pyselfi

Installing from pip gives access to the latest stable version of the code, but not to the features under development, which are only available from the `GitHub repository <https://github.com/florent-leclercq/pyselfi>`_.

Dependencies
~~~~~~~~~~~~

pySELFI depends on several other Python packages, which have their own dependencies: ``numpy``, ``scipy``, ``h5py``, ``matplotlib``, ``pathlib``, and ``sklearn``. They are listed in ``pyproject.toml`` and should be automatically installed by pip, if necessary.

Installation from sources
-------------------------

The sources for pySELFI can be downloaded from the `GitHub repository <https://github.com/florent-leclercq/pyselfi>`_. To download and install pySELFI from sources, clone the repository and install from there:

.. code-block:: console

    git clone git@github.com:florent-leclercq/pyselfi.git
    cd pyselfi
    pip install .

This will install pySELFI along with its dependencies. Note that the dot at the end of the last command means the current folder.

Dependencies and data for pyselfi_examples
------------------------------------------

Additional dependencies
~~~~~~~~~~~~~~~~~~~~~~~

* Running the `GRF and Simbelmynë examples <../examples/power_spectrum_inference.html#examples-with-blackboxes-using-simbelmyne>`_ additionally requires the ``pysbmy`` package, as provided by the `Simbelmynë <http://simbelmyne.florent-leclercq.eu>`_ code, version 0.2.0 or greater. Please refer to the `Simbelmynë documentation <https://simbelmyne.readthedocs.io/en/latest/>`_ for installation instructions. For the moment, the installation of ``pysbmy`` is not done automatically. Please make sure there is no version mismatch between the dependencies of Simbelmynë and the dependencies of pySELFI.
* From version 2.0, pySELFI includes an example using a Lotka-Volterra simulator and subsequent inference using likelihood-free rejection sampling. Running this example (found in ``src/pyselfi_examples/lotkavolterra/``) additionally requires `lotkavolterra_simulator <https://github.com/florent-leclercq/lotkavolterra_simulator>`_ (version 1.0) and `ELFI <https://github.com/elfi-dev/elfi>`_ (version 0.8.3). They are installed automatically when specifying ``lotkavolterra_example`` as an "extra" when installing the package, using

    .. code-block:: console

        pip install pyselfi[lotkavolterra_example]

Simulation data
~~~~~~~~~~~~~~~

Due to their size, the raw simulation data used in the pyselfi examples are not stored in the git repository. They are publicly available from the Aquila Consortium website at https://cloud.aquila-consortium.org/s/H4jeyqdtq2dTMCa. The files should be placed in the relevant directory (e.g. ``pyselfi/src/pyselfi_examples/grf/sims/``).

Developer installation from sources
-----------------------------------

If you plan to develop pySELFI (and contributions are welcome!), the best solution is to fork the repository on GitHub, clone your fork, and install it in editable mode:

.. code-block:: console

    git clone git@github.com:your-username/pyselfi.git
    cd pyselfi
    git remote add distant git@github.com:florent-leclercq/pyselfi.git
    pip install -e .

or

.. code-block:: console

    pip install -e .[lotkavolterra_example]

The flag ``-e`` will make pip install the project in editable mode.
