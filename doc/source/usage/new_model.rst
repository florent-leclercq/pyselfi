Specifying a new statistical model
===================================

.. _citepLeclercqetal2019: https://arxiv.org/abs/1902.10149

.. |citepLeclercqetal2019| replace:: Leclercq *et al*. 2019


pySELFI makes it possible to include new statistical models, including in particular new (Gaussian) priors and new methods (such as prior optimization, described in |citepLeclercqetal2019|_) for the core architecture.

The code for new statistical models can go into new python modules, placed for example in a subdirectory of ``pyselfi/`` such as ``pyselfi/power_spectrum/``, which is provided within the code.

Adding a new prior can be done by defining a new class :class:`prior`. This class shall have a method :meth:`compute` (not taking any argument except **self**). The only requirement is that after :meth:`compute` is called by the pySELFI core infrastructure, three attributes shall have been defined (by the constructor, by methods called by :meth:`compute`, or by :meth:`compute` itself):

* :obj:`self.mean`: an array of size :math:`S` giving the prior mean (recall that in the SELFI algorithm, the prior mean shall also be the expansion point :math:`\boldsymbol{\theta}_0`, so consistency is required here),
* :obj:`self.covariance`: an array of size :math:`S \times S` giving the prior covariance matrix :math:`\textbf{S}`,
* :obj:`self.inv_covariance`: an array of size :math:`S \times S` giving the prior inverse covariance matrix :math:`\textbf{S}^{-1}`,

.. currentmodule:: pyselfi.selfi

The :class:`prior` class shall also have a method :meth:`save` (with arguments **self, fname**) and a class method :meth:`load` (with arguments **cls, fname**) to save to and load prior data from the pySELFI hdf5 output file (in particular the three mandatory attributes, see above). Instances of the new :class:`prior` class can then be used to define and manipulate :obj:`selfi` objects, as described in `this page <../gettingstarted/basic_usage.html>`__. An example can be found in ``pyselfi/power_spectrum/prior.py``.

Adding new features to pySELFI can simply be done by defining a child of the :class:`selfi` class and using it instead of the parent class. An example can be found in ``pyselfi/power_spectrum/selfi.py``.

.. currentmodule:: pyselfi

Prior class template
--------------------

A template for a new :class:`prior` class is provided below and downloadable here: :download:`new_blackbox.py <https://github.com/florent-leclercq/pyselfi/tree/master/pyselfi_examples/templates/new_prior.py>`.

.. literalinclude :: ../../../pyselfi_examples/templates/new_prior.py


Selfi class template
--------------------

A template for a new :class:`child_selfi` class is provided below and downloadable here: :download:`new_selfi.py <https://github.com/florent-leclercq/pyselfi/tree/master/pyselfi_examples/templates/new_selfi.py>`.

.. literalinclude :: ../../../pyselfi_examples/templates/new_selfi.py
