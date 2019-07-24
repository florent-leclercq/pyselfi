Primordial power spectrum inference
===================================

.. _Leclercqetal2019: https://arxiv.org/abs/1902.10149

.. |Leclercqetal2019| replace:: Leclercq *et al*. (2019)

.. _citepLeclercqetal2019: https://arxiv.org/abs/1902.10149

.. |citepLeclercqetal2019| replace:: Leclercq *et al*. 2019

.. currentmodule:: pyselfi

The primordial power spectrum inference model defines a :obj:`prior` (:class:`power_spectrum.prior`) and a :obj:`selfi` child class (:class:`power_spectrum.selfi`), which allows optimization of prior hyperparameters (see sections II.C and II.E in |citepLeclercqetal2019|_, respectively).

The prior is characterized by 3 hyperparameters (:math:`\theta_\mathrm{norm}`, :math:`k_\mathrm{corr}` and :math:`\alpha_\mathrm{cv}`), and is used like this:

.. code-block:: python

    from pyselfi.power_spectrum.prior import power_spectrum_prior
    theta_norm=0.20
    k_corr=0.015
    alpha_cv=0.0008847681471875314
    prior=power_spectrum_prior(k_s,theta_0,theta_norm,k_corr,alpha_cv)

**k_s** is a vector containing the support wavenumbers, and **theta_0** is the desired prior mean (which should correspond to the expansion point for the effective likelihood). They shall be both of size :math:`S`.

The model is used like this:

.. code-block:: python

    from pyselfi.power_spectrum.selfi import power_spectrum_selfi
    selfi=power_spectrum_selfi(fname,pool_prefix,pool_suffix,prior,blackbox,
                               theta_0,Ne,Ns,Delta_theta,phi_obs)

Once the effective likelihood has been characterized, the prior can be later modified without having to rerun any simulation, for example:

.. code-block:: python
    
    theta_norm=0.05
    k_corr=0.010
    selfi.prior.theta_norm=theta_norm
    selfi.prior.k_corr=k_corr
    selfi.compute_prior()
    selfi.compute_posterior()
    
The prior can be optimized following the procedure described in |Leclercqetal2019|_ section II.E, using:

.. code-block:: python

    selfi.optimize_prior(theta_fiducial,k_opt_min,k_opt_max,x0=x0,
                    theta_norm_min=theta_norm_min,theta_norm_max=theta_norm_max,
                    theta_norm_mean=theta_norm_mean, theta_norm_std=theta_norm_std,
                    k_corr_min=k_corr_min,k_corr_max=k_corr_max,
                    k_corr_mean=k_corr_mean,k_corr_std=k_corr_std)
    selfi.save_prior()
    selfi.save_posterior()
    
For more details see the example notebooks (`Gaussian random field <https://nbviewer.jupyter.org/github/florent-leclercq/pyselfi/blob/master/pyselfi_examples/grf/selfi_GRF.ipynb>`__, `Galaxy survey <https://nbviewer.jupyter.org/github/florent-leclercq/pyselfi/blob/master/pyselfi_examples/simbelmyne/selfi_SBMY.ipynb>`__) and check the `API documentation <../api/inference_api.html#primordial-power-spectrum>`__ (:class:`pyselfi.power_spectrum.prior`, :class:`pyselfi.power_spectrum.selfi`).



