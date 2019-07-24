Examples with blackboxes using Simbelmynë
=========================================

Two blackboxes using the `Simbelmynë code <http://simbelmyne.florent-leclercq.eu>`__ are available in the package :mod:`pyselfi_examples` installed with the code:

* `GRF <#gaussian-random-fields>`__ - A blackbox to generate Gaussian random fields and evaluate their power spectrum,
* `SBMY <#synthetic-galaxy-observations>`__ - A blackbox to generate synthetic galaxy observations and evaluate their power spectrum.

Gaussian random fields
----------------------

.. currentmodule:: pyselfi_examples.grf.model

The blackbox is defined as an instance of the class :class:`blackbox_GRF.blackbox`:

.. code-block:: python

    from pyselfi_examples.grf.model.blackbox_GRF import blackbox
    blackbox=blackbox(P,theta2P=theta2P,k_s=k_s,G_sim=G_sim,P_ss=P_ss,
                      corner0=corner0,corner1=corner1,corner2=corner2,a_min=a_min,
                      noise_std=noise_std,seedphases=seedphases,fixphases=fixphases,
                      seednoise=seednoise,fixnoise=fixnoise,
                      save_frequency=save_frequency)


See the `API documentation <../api/blackboxes/grf.html#module-pyselfi_examples.grf.model.blackbox_GRF>`__ for the description of arguments, and the file ``pyselfi_examples/grf/model/setup_GRF.py`` for an example.

.. currentmodule:: pyselfi_examples.grf.model.blackbox_GRF.blackbox

Synthetic observations to be used in the inference can be generated using the method :meth:`make_data`:

.. code-block:: python

    phi_obs=blackbox.make_data(cosmo_obs,fname_input_power.h5,seedphasesobs,seednoiseobs,force)

The :obj:`blackbox` object also contains the method :meth:`compute_pool` with the necessary prototype (see `this page <../usage/new_blackbox.html>`__), as well as a method :meth:`evaluate` to generate individual realizations from any vector of input parameters :math:`\boldsymbol{\theta}`.

A full worked-out example using the GRF blackbox is available in `this notebook <https://nbviewer.jupyter.org/github/florent-leclercq/pyselfi/blob/master/pyselfi_examples/grf/selfi_GRF.ipynb>`__.

Synthetic galaxy observations
-----------------------------

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

See the `API documentation <../api/blackboxes/simbelmyne.html#module-pyselfi_examples.simbelmyne.model.blackbox_SBMY>`__ for the description of arguments, and the file ``pyselfi_examples/simbelmyne/model/setup_SBMY.py`` for an example.

.. currentmodule:: pyselfi_examples.simbelmyne.model.blackbox_SBMY.blackbox

A survey geometry file is required. A mock survey geometry file can be generated using the method :meth:`make_survey_geometry`:

.. code-block:: python

    blackbox.make_survey_geometry(N_CAT,cosmo,force)

Synthetic observations (with an index **i**) to be used in the inference can be generated using the method :meth:`make_data`:

.. code-block:: python

    phi_obs=blackbox.make_data(cosmo_obs,i,
                               force_powerspectrum,force_parfiles,
                               force_sim,force_mock,force_cosmo)

The :obj:`blackbox` object also contains the method :meth:`compute_pool` with the necessary prototype (see `this page <../usage/new_blackbox.html>`__), as well as a method :meth:`evaluate` to generate individual realizations from any vector of input parameters :math:`\boldsymbol{\theta}`.

A full worked-out example using the SBMY blackbox is available in `this notebook <https://nbviewer.jupyter.org/github/florent-leclercq/pyselfi/blob/master/pyselfi_examples/simbelmyne/selfi_SBMY.ipynb>`__.
