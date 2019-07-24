Basic usage
============

.. _Leclercqetal2019: https://arxiv.org/abs/1902.10149

.. |Leclercqetal2019| replace:: Leclercq *et al*. (2019)

.. _citepLeclercqetal2019: https://arxiv.org/abs/1902.10149

.. |citepLeclercqetal2019| replace:: Leclercq *et al*. 2019

pySELFI requires 3 basic ingredients:

* a **blackbox simulator**, which is used as the data model. It is represented as a :obj:`blackbox` object. See `this page <../usage/new_blackbox.html>`__ for instructions on how to define a new blackbox.
* a **statistical model**, which includes in particular the mean and covariance matrix of the (Gaussian) prior. A statistical model is characterized by a :obj:`prior` object, and (if necessary) an instance of a child class of the :class:`pyselfi.selfi` class. See `this page <../usage/new_model.html>`__ for instructions on how to define a new model.
* the **core infrastructure**, which computes the effective likelihood and the effective posterior. It should not be necessary to modify this part of the code, which is located in the root of the directory ``pyselfi/``.

The basic usage then consists of 5 steps:

1. `Setting up the prior, blackbox, and observations`_
2. `Computing the prior`_
3. `Running the necessary simulations`_
4. `Computing the effective likelihood`_
5. `Computing the effective posterior`_

Setting up the prior, blackbox, and observations
------------------------------------------------

The first step is to setup

* a :obj:`prior` object (see `this page <power_spectrum_inference.html>`__ for how to use the power spectrum prior defined in |citepLeclercqetal2019|_, or `this page <../usage/new_model.html>`__ for how to setup a new prior - if the model has :math:`S` parameters, the prior should have dimension :math:`S` too),
* a :obj:`blackbox` object (see `this page <examples.html>`__ for how to use `Simbelmynë <http://simbelmyne.florent-leclercq.eu>`_ as a blackbox, or `this page <../usage/new_blackbox.html>`__ for how to setup a new blackbox),
* and the observed summaries vector **phi_obs** (a vector of dimension :math:`P` if the blackbox returns vectors of dimension :math:`P`).

After that, define the :obj:`selfi` object using:

.. code-block:: python

    from pyselfi.selfi import selfi
    selfi=selfi(fname,pool_prefix,pool_suffix,prior,blackbox,theta_0,Ne,Ns,Delta_theta,phi_obs)

or an instance of a child class, such as:

.. code-block:: python

    from pyselfi.power_spectrum.selfi import power_spectrum_selfi
    selfi=power_spectrum_selfi(fname,pool_prefix,pool_suffix,prior,blackbox,
                               theta_0,Ne,Ns,Delta_theta,phi_obs)

**fname** corresponds to the name of an hdf5 file where results will be saved (and from which they can be loaded if they have been precomputed). **pool_prefix** and **pool_suffix** are the prefix and suffix of simulation pool filenames, so that a pool is identified by ``pool_prefix{d}poolsufffix`` where :math:`d` runs from :math:`0` to :math:`S` (:math:`0` for the expansion point and :math:`1` to :math:`S` for all possible directions in parameter space). **theta_0** is the expansion point (a vector of size :math:`S`). **Ne** is the number of required simulations at the expansion point (:math:`N_0` in the notations of |citepLeclercqetal2019|_), **Ns** the number of required simulations at each of the expansion points (:math:`N_s` in the notations of |citepLeclercqetal2019|_), and **Delta_theta** the step size for finite differencing (:math:`h` in the notations of |citepLeclercqetal2019|_).

For more details, see the `API documentation <../api/selfi.html#pyselfi.selfi.selfi>`_.

Computing the prior
--------------------

The second step is to compute and save the prior (mean, covariance matrix :math:`\textbf{S}` and inverse of the covariance matrix :math:`\textbf{S}^{-1}`). This is achieved using:

.. code-block:: python

    selfi.compute_prior()
    selfi.save_prior()

If the prior has already been computed and saved to the pySELFI output file, it can instead simply be loaded using:

.. code-block:: python

    selfi.load_prior()

The prior mean and covariance matrix can then be accessed using respectively **selfi.prior.mean** and **selfi.prior.covariance**.
    
Running the necessary simulations
---------------------------------

The third step is to run or load the necessary simulations, it is typically the expensive step. This is achieved using:

.. code-block:: python

    selfi.run_simulations()

Since pySELFI uses the concept of *simulation pool*, it will first try to load as many simulations as possible from the specified pool files, then run any additional simulation required. Therefore, if all simulations have been precomputed, this step only consists in loading the results.

Computing the effective likelihood
----------------------------------

The fourth step is to compute and save the effective likelihood (covariance matrix of the summaries at the expansion point :math:`\textbf{C}_0` and its inverse :math:`\textbf{C}_0^{-1}`, and gradient of the average blackbox :math:`\boldsymbol{\nabla}\textbf{f}_0`). This is achieved using:

.. code-block:: python

    selfi.compute_likelihood()
    selfi.save_likelihood()

This step assumes that the necessary simulations are accessible.
    
If the effective likelihood has already been computed and saved to the pySELFI output file, it can instead simply be loaded using:

.. code-block:: python

    selfi.load_likelihood()
    
This covariance matrix at the expansion point and the gradient of the blackbox are then accessed using respectively **selfi.likelihood.C_0** and **selfi.likelihood.grad_f**.

Computing the effective posterior
---------------------------------

The fifth step is to compute and save the effective posterior (mean :math:`\boldsymbol{\gamma}` and covariance matrix :math:`\boldsymbol{\Gamma}`). This is achieved using:

.. code-block:: python

    selfi.compute_posterior()
    selfi.save_posterior()

This step assumes that the prior and effective likelihood have been computed.

If the effective likelihood has already been computed and saved to the pySELFI output file, it can instead simply be loaded using:

.. code-block:: python

    selfi.load_posterior()

The effective posterior mean and covariance matrix can then be accessed using respectively **selfi.posterior.mean** and **selfi.posterior.covariance**: these are the end products of pySELFI.
