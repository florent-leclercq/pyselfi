Installation
============

This page contains instructions for installing the **pySELFI** code. Two python packages are installed: :mod:`pyselfi` and :mod:`pyselfi_examples`.

Minimal installation using pip
------------------------------

pySELFI requires Python 3.5 or greater. It can be installed easily from the command line using `pip <http://www.pip-installer.org/>`_:

.. code-block:: console

    pip install pyselfi

Installing from pip gives access to the latest stable version of the code, but not to the features under development, which are only available from the `git repository`__.

Dependencies
++++++++++++++

pySELFI depends on several other Python packages, which have their own dependencies: ``numpy``, ``scipy``, ``h5py``, ``pathlib``, and ``sklearn``. They should be automatically installed by pip, if necessary.

__ 

Installation from sources
-------------------------

The sources for pySELFI can be downloaded from the `GitHub repository <https://github.com/florent-leclercq/pyselfi>`_. To download and install pySELFI from sources, clone the repository and install from there:

.. code-block:: console

    git clone git@github.com:florent-leclercq/pyselfi.git
    cd pyselfi
    pip install .

This will install pySELFI along with its dependencies. Note that the dot at the end of the last command means the current folder.

Additional dependency: pysbmy
+++++++++++++++++++++++++++++

Running the `examples <examples.html>`_ additionally requires the ``pysbmy`` package, as provided by the `Simbelmynë <http://simbelmyne.florent-leclercq.eu>`_ code, version 0.2.0 or greater. Please refer to the `Simbelmynë documentation <https://simbelmyne.readthedocs.io/en/latest/>`_ for installation instructions. Please make sure there is no version mismatch between the dependencies of Simbelmynë and the dependencies of pySELFI.


Developer installation from sources
-----------------------------------

If you plan to develop pySELFI (and contributions are welcome!), the best solution is to fork the repository on GitHub, clone your fork, and install it in editable mode:

.. code-block:: console

    git clone git@github.com:your-username/pyselfi.git
    cd pyselfi
    git remote add distant git@github.com:florent-leclercq/pyselfi.git
    pip install -e .

The flag ``-e`` will make pip install the project in editable mode.
