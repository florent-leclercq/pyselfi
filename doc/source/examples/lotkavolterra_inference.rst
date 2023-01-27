Lotka-Volterra system inference
===============================

.. _Leclercqetal2019: https://arxiv.org/abs/1902.10149

.. |Leclercqetal2019| replace:: Leclercq *et al*. (2019)

.. _citepLeclercqetal2019: https://arxiv.org/abs/1902.10149

.. |citepLeclercqetal2019| replace:: Leclercq *et al*. 2019

.. _Leclercq2022: https://arxiv.org/abs/2209.11057

.. |Leclercq2022| replace:: Leclercq (2022)

.. _citepLeclercq2022: https://arxiv.org/abs/2209.11057

.. |citepLeclercq2022| replace:: Leclercq 2022

.. currentmodule:: pyselfi

The Lotka-Volterra model
------------------------

The Lotka-Volterra model (found in ``src/pyselfi/lotkavolterra/``) defines a :obj:`prior` (:class:`lotkavolterra.prior`) and a :obj:`selfi` child class (:class:`lotkavolterra.selfi`). The :class:`lotkavolterra.selfi` class allows optimization of prior hyperparameters (see sections II.C and II.E in |citepLeclercqetal2019|_, respectively).

The prior for time-dependent smooth functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The class :class:`lotkavolterra.prior` introduces a prior for time-dependent smooth functions. Let :math:`x(t_i)` be the number function of prey and :math:`y(t_i)` be the number function of predators at timesteps :math:`t_i`, with initially :math:`x_0 \equiv x(t_0)` and :math:`y_0 \equiv y(t_0)`.

The prior on :math:`\left\lbrace \left\lbrace x(t_i) \right\rbrace, \left\lbrace y(t_i) \right\rbrace \right\rbrace` is a Gaussian with mean :math:`\boldsymbol{\theta}_0` (a fiducial population number function) and covariance matrix :math:`\textbf{S} = \alpha_\mathrm{norm}^2 \textbf{V} \circ \textbf{K}` (where :math:`\circ` is the Hadamard product), such that:
:math:`\textbf{K} = \begin{pmatrix} \textbf{K}_x & \textbf{0} \\ \textbf{0} & \textbf{K}_y \end{pmatrix}`, :math:`(\textbf{K}_z)_{ij} = -\dfrac{1}{2} \left( \dfrac{t_i-t_j}{t_\mathrm{smooth}} \right)^2` for :math:`z \in \left\lbrace x,y \right\rbrace`, :math:`\textbf{V} = \begin{pmatrix} x_0 \textbf{u}\textbf{u}^\intercal & \textbf{0} \\ \textbf{0} & y_0 \textbf{u}\textbf{u}^\intercal \end{pmatrix}`, :math:`(\textbf{u})_i = 1 + \dfrac{t_i}{t_\mathrm{chaos}}`.

The prior is characterized by 3 hyperparameters:

    - :math:`\alpha_\mathrm{norm}` is the overall amplitude,
    - :math:`t_\mathrm{smooth}` is a typical time characterizing the smoothness of the population functions,
    - :math:`t_\mathrm{chaos}` is a typical time characterizing the chaotic behaviour of the system.

It is used like this:

.. code-block:: python

    from pyselfi.lotkavolterra.prior import lotkavolterra_prior
    alpha_norm = 0.02
    t_smooth = 1.6
    t_chaos = 8.2
    prior=lotkavolterra_prior(t_s,X0,Y0,theta_0,alpha_norm,t_smooth,t_chaos)

**t_s** is a vector containing the support wavenumbers, **X0** and **Y0** are the initial conditions, and **theta_0** is the desired prior mean (which should correspond to the expansion point for the effective likelihood). **t_s** and **theta_0** shall be both of size :math:`S`.

The SELFI model
~~~~~~~~~~~~~~~

The model is used like this:

.. code-block:: python

    from pyselfi.lotkavolterra.selfi import lotkavolterra_selfi
    selfi=lotkavolterra_selfi(fname,pool_prefix,pool_suffix,prior,blackbox,
                              theta_0,Ne,Ns,Delta_theta,phi_obs,LV)

where **LV** is a :class:`LVsimulator` (introduced by `lotkavolterra_simulator <https://github.com/florent-leclercq/lotkavolterra_simulator>`_).

The main computations (calculation of the effective likelihood) are done using:

.. code-block:: python

    selfi.run_simulations()
    selfi.compute_likelihood()
    selfi.save_likelihood()

Once the effective likelihood has been characterized, the prior can be later modified without having to rerun any simulation, for example:

.. code-block:: python

    alpha_norm = 0.05
    t_smooth = 2.1
    selfi.prior.alpha_norm=alpha_norm
    selfi.prior.t_smooth=t_smooth
    selfi.compute_prior()
    selfi.save_prior()
    selfi.compute_posterior()
    selfi.save_posterior()

The prior can be optimized using a set of fiducial simulations **theta_fid**:

.. code-block:: python

    selfi.optimize_prior(theta_fid, x0,
                         alpha_norm_min, alpha_norm_max,
                         t_smooth_min, t_smooth_max,
                         t_chaos_min, t_chaos_max)
    selfi.save_prior()


Check for model misspecification and score compression
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Mahalanobis distance between the posterior mean and the prior can be used to check for model misspecification (see |citepLeclercq2022|_). It is obtained using

.. code-block:: python

    Mahalanobis = selfi.Mahalanobis_ms_check()

A score compressor for top-level parameters can be obtained without having to rerun any simulation. It is obtained using

.. code-block:: python

    selfi.compute_grad_f_omega(t, t_s, tres, Delta_omega)
    selfi.compute_Fisher_matrix()
    omega_tilde = selfi.score_compression(phi_obs)

For more details see the example notebook (`selfi_LV.ipynb <https://nbviewer.jupyter.org/github/florent-leclercq/pyselfi/blob/master/src/pyselfi_examples/lotkavolterra/selfi_LV.ipynb>`__) and check the `API documentation <../api/inference_api.html#lotka-volterra-system>`__ (:class:`pyselfi.lotkavolterra.prior`, :class:`pyselfi.lotkavolterra.selfi`).


Inference of top-level parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Fisher-Rao distance between two points in top-level parameter space can be obtained using:

.. code-block:: python

    Fdist = selfi.FisherRao_distance(omega_tilde1, omega_tilde2)

It can be used for example, within likelihood-free rejection sampling (Approximate Bayesian Computation, ABC). See the example notebook `ABC_LV.ipynb <https://nbviewer.jupyter.org/github/florent-leclercq/pyselfi/blob/master/src/pyselfi_examples/lotkavolterra/ABC_LV.ipynb>`__.

Blackbox using lotkavolterra_simulator
--------------------------------------

.. currentmodule:: pyselfi_examples.lotkavolterra.model

A blackbox using the `lotkavolterra_simulator <https://github.com/florent-leclercq/lotkavolterra_simulator>`_ code is available in the package :mod:`pyselfi_examples` installed with the code. It is defined as an instance of the class :class:`blackbox_LV.blackbox` (see |citepLeclercq2022|_ for details):

.. code-block:: python

    from pyselfi_examples.lotkavolterra.model.blackbox_LV import blackbox
    # correct data model
    bbA=blackbox(P, X0=X0, Y0=Y0, seed=seed, fixnoise=fixnoise, model=0,
                 Xefficiency=Xefficiency0, Yefficiency=Yefficiency0,
                 threshold=threshold, mask=mask, obsP=obsP,
                 obsQ=obsQ, obsR=obsR0, obsS=obsS, obsT=obsT)
    # misspecified data model
    bbB=blackbox(P, X0=X0, Y0=Y0, seed=seed, fixnoise=fixnoise, model=1,
                 mask=mask, obsR=obsR1)


See the `API documentation <../api/blackboxes/lotkavolterra.html>`__ for the description of arguments, and the file ``src/pyselfi_examples/lotkavolterra/model/setup_LV.py`` for an example.

.. currentmodule:: pyselfi_examples.lotkavolterra.model.blackbox_LV.blackbox

Synthetic observations to be used in the inference can be generated using the method :meth:`make_data`:

.. code-block:: python

    phi_obs = bbA.make_data(theta_true, seed=seed)

The :obj:`blackbox` object also contains the method :meth:`compute_pool` with the necessary prototype (see `this page <../usage/new_blackbox.html>`__), as well as a method :meth:`evaluate` to generate individual realizations from any vector of input parameters :math:`\boldsymbol{\theta}`.

A full worked-out example using the Lotka-Volterra blackbox is available in `this notebook <https://nbviewer.jupyter.org/github/florent-leclercq/pyselfi/blob/master/src/pyselfi_examples/lotkavolterra/selfi_LV.ipynb>`__.
