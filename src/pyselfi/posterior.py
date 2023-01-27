#!/usr/bin/env python
#-------------------------------------------------------------------------------------
# pySELFI v2.0 -- src/pyselfi/posterior.py
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

"""Routines related to the SELFI posterior.

.. _Leclercqetal2019: https://arxiv.org/abs/1902.10149

.. |Leclercqetal2019| replace:: Leclercq *et al*. (2019)
"""

__author__  = "Florent Leclercq"
__version__ = "2.0"
__date__    = "2018-2023"
__license__ = "GPLv3"

class posterior(object):
    """This class represents the SELFI posterior. See equations (25) and (26) in |Leclercqetal2019|_ for expressions.

    Attributes
    ----------
    theta_0 : array, double, dimension=S
        prior mean and expansion point
    prior_covariance : array, double, dimension=(S,S)
        prior covariance matrix
    inv_prior_covariance : array, double, dimension=(S,S)
        inverse prior covariance matrix
    f_0 : array, double, dimension=P
        mean blackbox at the expansion point
    C_0 : array, double, dimension=(P,P)
        covariance matrix of summaries at the expansion point
    inv_C_0 : array, double, dimension=(P,P)
        inverse covariance matrix of summaries at the expansion point
    grad_f : array, double, dimension=(S,P)
        gradient of the blackbox at the expansion point
    phi_obs : array, double, dimension=P
        observed summaries vector

    """

    # Initialization
    def __init__(self,theta_0,prior_covariance,inv_prior_covariance,f_0,C_0,inv_C_0,grad_f,phi_obs):
        """Initializes a posterior object.
        """
        self.theta_0=theta_0
        self.prior_covariance=prior_covariance
        self.inv_prior_covariance=inv_prior_covariance
        self.f_0=f_0
        self.C_0=C_0
        self.inv_C_0=inv_C_0
        self.grad_f=grad_f
        self.phi_obs=phi_obs
        self.EPS_inv_covariance=1e-7
        self.EPS_inv_N=1e-7
        self.EPS_SplusN=1e-7
        self.EPS_Gamma=1e-7
        self.EPS_residual=1e-3

    def _compute_inv_N(self):
        """Computes the inverse of the 'noise' matrix.
        """
        import numpy as np
        inv_C_0=self.inv_C_0
        grad_f=self.grad_f

        self.inv_N = inv_N = grad_f.T.dot(inv_C_0).dot(grad_f)

    def _compute_inverse_posterior_covariance(self):
        """Computes the inverse posterior covariance.
        """
        inv_prior_covariance=self.inv_prior_covariance
        self._compute_inv_N()
        inv_N=self.inv_N

        self.inv_covariance = inv_N + inv_prior_covariance

    def _get_posterior_covariance(self):
        """Gets the posterior covariance. See equation (26) in |Leclercqetal2019|_.

        Returns
        -------
        Gamma : array, double, dimension=(S,S)
            posterior covariance matrix

        """
        import numpy as np
        from pyselfi.utils import regular_inv
        self._compute_inverse_posterior_covariance()

        inv_covariance=self.inv_covariance
        Gamma = regular_inv(inv_covariance, self.EPS_inv_covariance, self.EPS_residual)

        return Gamma

    def _get_posterior_covariance_alt(self):
        """Gets the posterior covariance. See equation (26) in |Leclercqetal2019|_. Alternative algebra: can be used if numerically more stable.

        Returns
        -------
        Gamma : array, double, dimension=(S,S)
            posterior covariance matrix

        """
        from pyselfi.utils import regular_inv
        S=self.prior_covariance
        self._compute_inv_N()
        inv_N=self.inv_N

        N = regular_inv(inv_N, self.EPS_inv_N, self.EPS_residual)
        Gamma = S-S.dot(regular_inv(S+N, self.EPS_SplusN, self.EPS_residual)).dot(S)
        self.inv_covariance = regular_inv(Gamma, self.EPS_Gamma, self.EPS_residual)
        return Gamma

    def _get_posterior_mean(self):
        """Gets the posterior mean. See equation (25) in |Leclercqetal2019|_.

        Returns
        -------
        gamma : array, double, dimension=S
            posterior mean

        """
        import scipy.linalg as sla
        theta_0=self.theta_0
        f_0=self.f_0
        inv_C_0=self.inv_C_0
        grad_f=self.grad_f
        phi_obs=self.phi_obs
        inv_covariance=self.inv_covariance

        j = grad_f.T.dot(inv_C_0).dot(phi_obs-f_0)
        theta_rec = sla.solve(inv_covariance, j)
        gamma = theta_0 + theta_rec
        return gamma

    def _get_posterior_mean_alt(self):
        """Gets the posterior mean. See equation (25) in |Leclercqetal2019|_. Alternative algebra: can be used if numerically more stable.

        Returns
        -------
        gamma : array, double, dimension=S
            posterior mean

        """
        theta_0=self.theta_0
        f_0=self.f_0
        inv_C_0=self.inv_C_0
        grad_f=self.grad_f
        phi_obs=self.phi_obs
        Gamma=self.covariance

        j = grad_f.T.dot(inv_C_0).dot(phi_obs-f_0)
        theta_rec = Gamma.dot(j)
        gamma = theta_0 + theta_rec
        return gamma

    def compute(self):
        """Computes the posterior (mean and covariance matrix).
        """
        self.covariance = self._get_posterior_covariance()
        self.mean = self._get_posterior_mean()

    def logpdf(self,theta,theta_mean,theta_covariance,theta_icov):
        """Returns the log posterior probability at a given point in parameter space. See equation (24) in |Leclercqetal2019|_.

        Parameters
        ----------
        theta : array, double, dimension=S
            evaluation point in parameter space
        theta_mean : array, double, dimension=S
            posterior mean
        theta_covariance : array, double, dimension=(S,S)
            posterior covariance
        theta_icov : array, double, dimension=(S,S)
            inverse posterior covariance

        Returns
        -------
        logpdf : double
            log posterior probability

        """
        import numpy as np
        diff = theta-theta_mean
        return -diff.dot(theta_icov).dot(diff)/2. - np.linalg.slogdet(2*np.pi*theta_covariance)[1]/2.

    def save(self,fname):
        """Saves the posterior to an output file.

        Parameters
        ----------
        fname : :obj:`str`
            output filename

        """
        import h5py as h5
        from ctypes import c_double
        from pyselfi.utils import PrintMessage, save_replace_dataset
        PrintMessage(3, "Writing posterior in data file '{}'...".format(fname))

        with h5.File(fname, 'r+') as hf:
            save_replace_dataset(hf, '/posterior/theta_0', data=self.theta_0, maxshape=(None,), dtype=c_double)
            save_replace_dataset(hf, '/posterior/f_0', data=self.f_0, maxshape=(None,), dtype=c_double)
            save_replace_dataset(hf, '/posterior/C_0', data=self.C_0, maxshape=(None, None), dtype=c_double)
            save_replace_dataset(hf, '/posterior/inv_C_0', data=self.inv_C_0, maxshape=(None, None), dtype=c_double)
            save_replace_dataset(hf, '/posterior/grad_f', data=self.grad_f, maxshape=(None, None), dtype=c_double)
            save_replace_dataset(hf, '/posterior/phi_obs', data=self.phi_obs, maxshape=(None,), dtype=c_double)
            save_replace_dataset(hf, '/posterior/prior_covariance', data=self.prior_covariance, maxshape=(None, None), dtype=c_double)
            save_replace_dataset(hf, '/posterior/inv_prior_covariance', data=self.inv_prior_covariance, maxshape=(None, None), dtype=c_double)
            save_replace_dataset(hf, '/posterior/inv_covariance', data=self.inv_covariance, maxshape=(None, None), dtype=c_double)
            save_replace_dataset(hf, '/posterior/mean', data=self.mean, maxshape=(None,), dtype=c_double)
            save_replace_dataset(hf, '/posterior/covariance', data=self.covariance, maxshape=(None, None), dtype=c_double)

        PrintMessage(3, "Writing posterior in data file '{}' done.".format(fname))

    @classmethod
    def load(cls,fname):
        """Loads the posterior from an input file.

        Parameters
        ----------
        fname : :obj:`str`
            input filename

        Returns
        -------
        posterior : :obj:`posterior`
            loaded posterior object

        """
        import h5py as h5
        import numpy as np
        from ctypes import c_double
        from pyselfi.utils import PrintMessage
        PrintMessage(3, "Reading posterior in data file '{}'...".format(fname))

        with h5.File(fname,'r') as hf:
            theta_0=np.array(hf.get('/posterior/theta_0'), dtype=c_double)
            prior_covariance=np.array(hf.get('/posterior/prior_covariance'), dtype=c_double)
            inv_prior_covariance=np.array(hf.get('/posterior/inv_prior_covariance'), dtype=c_double)
            f_0=np.array(hf.get('/posterior/f_0'), dtype=c_double)
            C_0=np.array(hf.get('/posterior/C_0'), dtype=c_double)
            inv_C_0=np.array(hf.get('/posterior/inv_C_0'), dtype=c_double)
            grad_f=np.array(hf.get('/posterior/grad_f'), dtype=c_double)
            phi_obs=np.array(hf.get('/posterior/phi_obs'), dtype=c_double)
            posterior=cls(theta_0,prior_covariance,inv_prior_covariance,f_0,C_0,inv_C_0,grad_f,phi_obs)
            posterior.inv_covariance=np.array(hf.get('/posterior/inv_covariance'), dtype=c_double)
            posterior.mean=np.array(hf.get('/posterior/mean'), dtype=c_double)
            posterior.covariance = np.array(hf.get('/posterior/covariance'), dtype=c_double)

        PrintMessage(3, "Reading posterior in data file '{}' done.".format(fname))
        return posterior
#end class(posterior)
