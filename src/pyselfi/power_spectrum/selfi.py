#!/usr/bin/env python
#-------------------------------------------------------------------------------------
# pySELFI v2.0 -- src/pyselfi/power_spectrum/selfi.py
# Copyright (C) 2019-2023 Florent Leclercq.
#
# This file is part of the pySELFI distribution
# (https://github.com/florent-leclercq/pyselfi/)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# The text of the license is located in the root directory of the source package.
#-------------------------------------------------------------------------------------

"""SELFI for power spectrum inference: main class.

.. _Leclercqetal2019: https://arxiv.org/abs/1902.10149

.. |Leclercqetal2019| replace:: Leclercq *et al*. (2019)
"""

__author__  = "Florent Leclercq"
__version__ = "2.0"
__date__    = "2018-2023"
__license__ = "GPLv3"

from pyselfi import selfi

class power_spectrum_selfi(selfi):
    """Main class for power spectrum inference with SELFI. Child of the main class of the SELFI code.
    """

    #################
    # Prior optimizer
    #################

    def logprior_hyperparameters(self,theta_norm,theta_norm_mean,theta_norm_std,k_corr,k_corr_mean,k_corr_std):
        """Returns the log pdf of the Gaussian prior for hyperparameters theta_norm and k_corr.

        Parameters
        ----------
        theta_norm : double
            theta_norm value where to evaluate the prior
        theta_norm_mean : double
            mean of the Gaussian prior on theta_norm
        theta_norm_std : double
            standard deviations of the Gaussian prior on theta_norm
        k_corr : double
            k_corr value where to evaluate the prior
        k_corr_mean : double
            mean of the Gaussian prior on k_corr
        k_corr_std : double
            standard deviations of the Gaussian prior on k_corr

        Returns
        -------
        logpdf : double
            log pdf of the Gaussian prior at the evaluation point

        """
        from scipy.stats import norm
        return norm.logpdf(theta_norm,theta_norm_mean,theta_norm_std)+norm.logpdf(k_corr,k_corr_mean,k_corr_std)

    def loglikelihood_hyperparams(self,theta_fiducial,Nbin_opt_min,Nbin_opt_max,theta_norm,k_corr):
        """Returns the log pdf of the likelihood for hyperparameters theta_norm and k_corr. See equation (27) in |Leclercqetal2019|_.

        Parameters
        ----------
        theta_fiducial : array, double, dimension=S
            fiducial point in parameter space to optimize the prior
        Nbin_opt_min : int
            minimum index in k_s considered for optimization
        Nbin_opt_max : int
            maximum index in k_s considered for optimization
        theta_norm : double
            theta_norm value where to evaluate the likelihood
        k_corr : double
            k_corr value where to evaluate the likelihood

        Returns
        -------
        logpdf : double
            log pdf of the likelihood at the evaluation point

        """
        self.prior.theta_norm=theta_norm
        self.prior.k_corr=k_corr
        self.compute_prior()
        self.compute_posterior()
        theta_mean, theta_covariance, theta_icov = self.restrict_posterior(Nbin_opt_min,Nbin_opt_max)
        theta_fiducial=theta_fiducial[Nbin_opt_min:Nbin_opt_max]
        return self.posterior.logpdf(theta_fiducial, theta_mean, theta_covariance, theta_icov)

    def logposterior_hyperparameters(self,theta_fiducial,Nbin_opt_min,Nbin_opt_max,theta_norm,k_corr,theta_norm_mean=0.2,theta_norm_std=0.3,k_corr_mean=0.025,k_corr_std=0.015):
        """Returns the log pdf of the unnormalized posterior for hyperparameters theta_norm and k_corr (logprior_hyperparameters + loglikelihood_hyperparams).

        Parameters
        ----------
        theta_fiducial : array, double, dimension=S
            fiducial point in parameter space to optimize the prior
        Nbin_opt_min : int
            minimum index in k_s considered for optimization
        Nbin_opt_max : int
            maximum index in k_s considered for optimization
        theta_norm : double
         theta_norm value where to evaluate the likelihood
        k_corr : double
            k_corr value where to evaluate the likelihood
        theta_norm_mean : double, optional, default=0.2
            mean of the Gaussian prior on theta_norm
        theta_norm_std : double, optional, default=0.3
            standard deviations of the Gaussian prior on theta_norm
        k_corr_mean : double, optional, default=0.025
            mean of the Gaussian prior on k_corr
        k_corr_std : double, optional, default=0.015
            standard deviations of the Gaussian prior on k_corr

        Returns
        -------
        logpdf : double
            log pdf of the unnormalized posterior at the evaluation point

        """
        logprior = self.logprior_hyperparameters(theta_norm,theta_norm_mean,theta_norm_std,k_corr,k_corr_mean,k_corr_std)
        loglikelihood = self.loglikelihood_hyperparams(theta_fiducial,Nbin_opt_min,Nbin_opt_max,theta_norm,k_corr)
        logposterior = logprior + loglikelihood
        return logposterior

    def optimize_prior(self, theta_fiducial, k_opt_min, k_opt_max, x0=None,
                       theta_norm_min=1e-5, theta_norm_max=5., theta_norm_mean=0.2,
                       theta_norm_std=0.3, k_corr_min=0.005, k_corr_max=0.05,
                       k_corr_mean=0.025, k_corr_std=0.015, method='L-BFGS-B', options={'maxiter':50, 'ftol':1e-20, 'gtol':1e-20, 'eps':1e-6, 'disp':True}):
        """Optimizes the SELFI prior. See section II.E. in |Leclercqetal2019|_.

        Parameters
        ----------
        theta_fiducial : array, double, dimension=S
            fiducial point in parameter space to optimize the prior
        k_opt_min : double
            minimum wavenumber to be used in prior optimization
        k_opt_max : double
            maximum wavenumber to be used in prior optimization
        x0 : array, double, dimension=2, optional, default=(prior.theta_norm,prior.k_corr))
            starting point in parameter space
        theta_norm_min : double, optional, default=1e-5
            boundary in parameter space
        theta_norm_max : double, optional, default=5.
            boundary in parameter space
        theta_norm_mean : double, optional, default=0.2
            mean of the Gaussian prior on theta_norm
        theta_norm_std : double, optional, default=0.3
            standard deviations of the Gaussian prior on theta_norm
        k_corr_min : double, optional, default=0.005
            boundary in parameter space
        k_corr_max : double, optional, default=0.05
            boundary in parameter space
        k_corr_mean : double, optional, default=0.025
            mean of the Gaussian prior on k_corr
        k_corr_std : double, optional, default=0.015
            standard deviations of the Gaussian prior on k_corr
        method : :obj:`str`, *optional, default='L-BFGS-B'*
            optimization method. See documentation of scipy.optimize.minimize
        options : dictionary, optional, default={'maxiter':50, 'ftol':1e-20, 'gtol':1e-20, 'eps':1e-6, 'disp':True}
            optimization options. See documentation of scipy.optimize.minimize

        """
        Nbin_opt_min=self.prior.Nbin_min(k_opt_min)
        Nbin_opt_max=self.prior.Nbin_max(k_opt_max)
        from scipy.optimize import minimize
        def potential(x):
            import numpy as np
            theta_norm, k_corr = x
            return -self.logposterior_hyperparameters(theta_fiducial,Nbin_opt_min,Nbin_opt_max,theta_norm,k_corr,theta_norm_mean,theta_norm_std,k_corr_mean,k_corr_std)

        x0=x0 or [self.prior.theta_norm,self.prior.k_corr]
        res = minimize(potential, x0, method=method,
               options=options, bounds=[[max(1e-10,theta_norm_min),theta_norm_max],[max(1e-10,k_corr_min),k_corr_max]])
        print(res)

        self.prior.theta_norm=res.x[0]
        self.prior.k_corr=res.x[1]
        self.compute_prior()
        self.compute_posterior()

    def restrict_prior(self,Nbin_min,Nbin_max):
        """Restricts the SELFI prior to some range of scales.

        Parameters
        ----------
        Nbin_min : int
            minimal index in k_s to be used
        Nbin_max : int
            maximal index in k_s to be used

        Returns
        -------
        theta_mean : array, double, dimension=S'
            restricted prior mean
        theta_covariance : array, double, dimension=(S',S')
            restricted prior covariance
        theta_icov : array, double, dimension=(S',S')
            restricted inverse prior covariance where S'=Nbin_max-Nbin_min

        """
        from pyselfi.utils import regular_inv
        theta_mean=self.prior.mean[Nbin_min:Nbin_max]
        theta_covariance=self.prior.covariance[Nbin_min:Nbin_max,Nbin_min:Nbin_max]
        theta_icov=regular_inv(theta_covariance)
        return theta_mean, theta_covariance, theta_icov

    def restrict_posterior(self,Nbin_min,Nbin_max):
        """Restricts the SELFI posterior to some range of scales.

        Parameters
        ----------
        Nbin_min : int
            minimal index in k_s to be used
        Nbin_max : int
            maximal index in k_s to be used

        Returns
        -------
        theta_mean : array, double, dimension=S'
            restricted posterior mean
        theta_covariance : array, double, dimension=(S',S')
            restricted posterior covariance
        theta_icov : array, double, dimension=(S',S')
            restricted inverse posterior covariance where S'=Nbin_max-Nbin_min

        """
        from pyselfi.utils import regular_inv
        theta_mean=self.posterior.mean[Nbin_min:Nbin_max]
        theta_covariance=self.posterior.covariance[Nbin_min:Nbin_max,Nbin_min:Nbin_max]
        theta_icov=regular_inv(theta_covariance)
        return theta_mean, theta_covariance, theta_icov
# end class(power_spectrum_selfi)
