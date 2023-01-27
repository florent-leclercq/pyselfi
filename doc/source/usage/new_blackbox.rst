Specifying a new blackbox simulator
===================================

One of the key features of the SELFI algorithm is its *likelihood-free* nature, which implies its ability to deal with arbitrary simulators as data models, i.e.  blackboxes.

In pySELFI, new blackbox simulators can be defined by the user in a new python module, as a class named :class:`blackbox`.

The constructor for :obj:`blackbox` objects shall initialize :math:`P` (the dimension of summaries space, only mandatory attribute), as well as any other attribute. Its prototype shall therefore be:

.. code-block:: python

    def __init__(self, P, **kwargs):
        self.P=P
        # Other initializations

.. currentmodule:: pyselfi.pool

The only requirement for a blackbox class to work with the pySELFI core infrastructure is a method named :meth:`compute_pool` with the following prototype and returning a :obj:`pool` object:

.. code-block:: python

        def compute_pool(self, theta, d, pool_fname, N):
            from pyselfi.pool import pool
            p=pool(pool_fname,N)
            # Compute the simulations going into the pool
            return p

The arguments (handed by pySELFI) are **theta** (a :math:`S`-dimensional vector of input parameters), **d** (the expansion index: :math:`0` for the expansion point and :math:`1` to :math:`S` for all possible directions in parameter space), the pool filename **pool_fname**, and the required number of blackbox realizations **N**.

.. currentmodule:: pyselfi.selfi

Instances of the new :class:`blackbox` class can then be used to define and manipulate :obj:`selfi` objects, as described in `this page <../gettingstarted/basic_usage.html>`__.

Working examples of :class:`blackbox` classes can be found in ``src/pyselfi_examples/grf/model/blackbox_GRF.py``, ``src/pyselfi_examples/simbelmyne/model/blackbox_SBMY.py``, or ``src/pyselfi_examples/lotkavolterra/model/blackbox_LV.py``.

Blackbox class template
-----------------------

A more extensive template for a new :class:`blackbox` class is provided below and downloadable here: :download:`new_blackbox.py <https://github.com/florent-leclercq/pyselfi/tree/master/src/pyselfi_examples/templates/new_blackbox.py>`.

.. literalinclude :: ../../../src/pyselfi_examples/templates/new_blackbox.py
