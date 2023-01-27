Primordial power spectrum inference
===================================

.. _Leclercqetal2019: https://arxiv.org/abs/1902.10149

.. |Leclercqetal2019| replace:: Leclercq *et al*. (2019)

.. _citepLeclercqetal2019: https://arxiv.org/abs/1902.10149

.. |citepLeclercqetal2019| replace:: Leclercq *et al*. 2019

.. currentmodule:: pyselfi

The primordial power spectrum model
-----------------------------------

The primordial power spectrum inference model (found in ``src/pyselfi/power_spectrum/``) defines a :obj:`prior` (:class:`power_spectrum.prior`) and a :obj:`selfi` child class (:class:`power_spectrum.selfi`). The :class:`power_spectrum.selfi` class allows optimization of prior hyperparameters (see sections II.C and II.E in |citepLeclercqetal2019|_, respectively).

The prior is characterized by 3 hyperparameters (:math:`\theta_\mathrm{norm}`, :math:`k_\mathrm{corr}` and :math:`\alpha_\mathrm{cv}`), and is used like this:

.. code-block:: python

    from pyselfi.power_spectrum.prior import power_spectrum_prior
    theta_norm=0.20
    k_corr=0.015
    alpha_cv=0.0008847681471875314
    prior=power_spectrum_prior(k_s,theta_0,theta_norm,k_corr,alpha_cv)

**k_s** is a vector containing the support wavenumbers, and **theta_0** is the desired prior mean (which should correspond to the expansion point for the effective likelihood). They shall be both of size :math:`S`.

The model is initialised like this:

.. code-block:: python

    from pyselfi.power_spectrum.selfi import power_spectrum_selfi
    selfi=power_spectrum_selfi(fname,pool_prefix,pool_suffix,prior,blackbox,
                               theta_0,Ne,Ns,Delta_theta,phi_obs)

The main computations (calculation of the effective likelihood) are done using:

.. code-block:: python

    selfi.run_simulations()
    selfi.compute_likelihood()
    selfi.save_likelihood()

Once the effective likelihood has been characterized, the prior can be later modified without having to rerun any simulation, for example:

.. code-block:: python

    theta_norm=0.05
    k_corr=0.010
    selfi.prior.theta_norm=theta_norm
    selfi.prior.k_corr=k_corr
    selfi.compute_prior()
    selfi.save_prior()
    selfi.compute_posterior()
    selfi.save_posterior()

The prior can be optimized following the procedure described in |Leclercqetal2019|_ section II.E, using:

.. code-block:: python

    selfi.optimize_prior(theta_fiducial,k_opt_min,k_opt_max,x0=x0,
                    theta_norm_min=theta_norm_min,theta_norm_max=theta_norm_max,
                    theta_norm_mean=theta_norm_mean, theta_norm_std=theta_norm_std,
                    k_corr_min=k_corr_min,k_corr_max=k_corr_max,
                    k_corr_mean=k_corr_mean,k_corr_std=k_corr_std)
    selfi.save_prior()

For more details see the example notebooks (`Gaussian random field <https://nbviewer.jupyter.org/github/florent-leclercq/pyselfi/blob/master/src/pyselfi_examples/grf/selfi_GRF.ipynb>`__, `Galaxy survey <https://nbviewer.jupyter.org/github/florent-leclercq/pyselfi/blob/master/src/pyselfi_examples/simbelmyne/selfi_SBMY.ipynb>`__) and check the `API documentation <../api/inference_api.html#primordial-power-spectrum>`__ (:class:`pyselfi.power_spectrum.prior`, :class:`pyselfi.power_spectrum.selfi`).

Blackboxes using Simbelmynë
---------------------------

Two blackboxes using the `Simbelmynë code <http://simbelmyne.florent-leclercq.eu>`__ are available in the package :mod:`pyselfi_examples` installed with the code:

* `GRF <#gaussian-random-fields>`__ - A blackbox to generate Gaussian random fields and evaluate their power spectrum,
* `SBMY <#synthetic-galaxy-observations>`__ - A blackbox to generate synthetic galaxy observations and evaluate their power spectrum.

Gaussian random fields
~~~~~~~~~~~~~~~~~~~~~~

.. currentmodule:: pyselfi_examples.grf.model

The blackbox is defined as an instance of the class :class:`blackbox_GRF.blackbox`:

.. code-block:: python

    from pyselfi_examples.grf.model.blackbox_GRF import blackbox
    blackbox=blackbox(P,theta2P=theta2P,k_s=k_s,G_sim=G_sim,P_ss=P_ss,
                      corner0=corner0,corner1=corner1,corner2=corner2,a_min=a_min,
                      noise_std=noise_std,seedphases=seedphases,fixphases=fixphases,
                      seednoise=seednoise,fixnoise=fixnoise,
                      save_frequency=save_frequency)


See the `API documentation <../api/blackboxes/grf.html#module-pyselfi_examples.grf.model.blackbox_GRF>`__ for the description of arguments, and the file ``src/pyselfi_examples/grf/model/setup_GRF.py`` for an example.

.. currentmodule:: pyselfi_examples.grf.model.blackbox_GRF.blackbox

Synthetic observations to be used in the inference can be generated using the method :meth:`make_data`:

.. code-block:: python

    phi_obs=blackbox.make_data(cosmo_obs,fname_input_power.h5,seedphasesobs,seednoiseobs,force)

The :obj:`blackbox` object also contains the method :meth:`compute_pool` with the necessary prototype (see `this page <../usage/new_blackbox.html>`__), as well as a method :meth:`evaluate` to generate individual realizations from any vector of input parameters :math:`\boldsymbol{\theta}`.

A full worked-out example using the GRF blackbox is available in `this notebook <https://nbviewer.jupyter.org/github/florent-leclercq/pyselfi/blob/master/src/pyselfi_examples/grf/selfi_GRF.ipynb>`__.

Synthetic galaxy observations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. currentmodule:: pyselfi_examples.simbelmyne.model

The blackbox is defined as an instance of the class :class:`blackbox_SBMY.blackbox`:

.. code-block:: python

    from pyselfi_examples.simbelmyne.model.blackbox_SBMY import blackbox
    blackbox=blackbox(P,theta2P=theta2P,k_s=k_s,G_sim=G_sim,G_ss=G_ss,P_ss=P_ss,
                      corner0=corner0,corner1=corner1,corner2=corner2,Np0=Np0,Npm0=Npm0,
                      fdir=fdir,fsimdir=fdir,fname_inputsurveygeometry=fname_inputsurveygeometry,
                      b_cut=b_cut,bright_apparent_magnitude_cut=bright_apparent_magnitude_cut,
                      faint_apparent_magnitude_cut=faint_apparent_magnitude_cut,
                      bright_absolute_magnitude_cut=bright_absolute_magnitude_cut,
                      faint_absolute_magnitude_cut=faint_absolute_magnitude_cut,
                      Mstar=Mstar,alpha=alpha,save_frequency=save_frequency)

See the `API documentation <../api/blackboxes/simbelmyne.html#module-pyselfi_examples.simbelmyne.model.blackbox_SBMY>`__ for the description of arguments, and the file ``src/pyselfi_examples/simbelmyne/model/setup_SBMY.py`` for an example.

.. currentmodule:: pyselfi_examples.simbelmyne.model.blackbox_SBMY.blackbox

A survey geometry file is required. A mock survey geometry file can be generated using the method :meth:`make_survey_geometry`:

.. code-block:: python

    blackbox.make_survey_geometry(N_CAT,cosmo,force)

Synthetic observations (with an index :math:`i`) to be used in the inference can be generated using the method :meth:`make_data`:

.. code-block:: python

    phi_obs=blackbox.make_data(cosmo_obs,i,
                               force_powerspectrum,force_parfiles,
                               force_sim,force_mock,force_cosmo)

The :obj:`blackbox` object also contains the method :meth:`compute_pool` with the necessary prototype (see `this page <../usage/new_blackbox.html>`__), as well as a method :meth:`evaluate` to generate individual realizations from any vector of input parameters :math:`\boldsymbol{\theta}`.

A full worked-out example using the SBMY blackbox is available in `this notebook <https://nbviewer.jupyter.org/github/florent-leclercq/pyselfi/blob/master/src/pyselfi_examples/simbelmyne/selfi_SBMY.ipynb>`__.
