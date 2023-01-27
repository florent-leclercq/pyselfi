#!/usr/bin/env python
#-------------------------------------------------------------------------------------
# pySELFI v2.0 -- src/pyselfi/power_spectrum/prior.py
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

"""Routines related to the SELFI power spectrum prior.

.. _Leclercqetal2019: https://arxiv.org/abs/1902.10149

.. |Leclercqetal2019| replace:: Leclercq *et al*. (2019)
"""

__author__  = "Florent Leclercq"
__version__ = "2.0"
__date__    = "2018-2023"
__license__ = "GPLv3"

class power_spectrum_prior(object):
    """This class represents the SELFI power spectrum prior. See equations (20)-(23) in |Leclercqetal2019|_ for expressions.

    Attributes
    ----------
    k_s : array, double, dimension=S
        array of support wavenumbers
    theta_0 : array, double, dimension=S
        prior mean
    theta_norm : double
        overall amplitude of the prior covariance matrix
    k_corr : double
        wavenumber of the length scale of prior correlations
    alpha_cv : double
        strength of cosmic variance
    log_kmodes : bool, optional, default=False
        take RBF in log(k) instead of k?

    """

    # Initialization
    def __init__(self,k_s,theta_0,theta_norm,k_corr,alpha_cv,log_kmodes=False):
        """Initializes a prior object.
        """
        self.k_s=k_s
        self.mean=theta_0
        self.theta_norm=theta_norm
        self.k_corr=k_corr
        self.alpha_cv=alpha_cv
        self.log_kmodes=log_kmodes
        self.EPS_K=1e-7
        self.EPS_residual=1e-3

    @property
    def gamma(self):
        """double: Defined by :math:`\gamma \equiv 1/(2 k_\mathrm{corr}^2)`
        """
        return 1/(2*self.k_corr**2)

    @property
    def gamma_log(self):
        """double: Defined by :math:`\gamma_\mathrm{log} \equiv 1/(2 \log k_\mathrm{corr}^2)`
        """
        import numpy as np
        return 1/(2*np.log(self.k_corr)**2)

    def Nbin_min(self,k_min):
        """Finds the index of the minimal wavenumber given k_min.

        Parameters
        ----------
        k_min : double
            minimal wavenumber

        Returns
        -------
        Nbin_min : int
            minimal index such that k_s[Nbin_min] >= k_min

        """
        import numpy as np
        return np.where(self.k_s>=k_min)[0].min()

    def Nbin_max(self,k_max):
        """Finds the index of the maximal wavenumber given k_max.

        Parameters
        ----------
        k_max : double
            maximal wavenumber

        Returns
        -------
        Nbin_max : int
            maximal index such that k_s[Nbin_max] <= k_max

        """
        import numpy as np
        return np.where(self.k_s<=k_max)[0].max()+1

    def _get_rbf(self):
        """Gets the radial basis function (RBF) part of the prior covariance matrix. See equation (20) in |Leclercqetal2019|_.

        Returns
        -------
        K : array, double, dimension=(S,S)
            RBF kernel

        """
        import numpy as np
        from sklearn.metrics.pairwise import rbf_kernel
        k_s=self.k_s
        if self.log_kmodes:
            K=rbf_kernel(np.log(k_s).reshape(-1,1), np.log(k_s).reshape(-1,1), gamma=self.gamma_log)
        else:
            K=rbf_kernel(k_s.reshape(-1,1), k_s.reshape(-1,1), gamma=self.gamma)

        return K

    def _get_cosmic_variance(self):
        """Gets the cosmic variance part of the prior covariance matrix, :math:`\\textbf{u}\\textbf{u}^\intercal`. See equations (21)-(22) in |Leclercqetal2019|_.

        Returns
        -------
        V : array, double, dimension=(S,S)
            cosmic variance matrix

        """
        import numpy as np
        k_s=self.k_s
        alpha_cv=self.alpha_cv
        S_k=alpha_cv/np.sqrt(k_s**3)
        V=np.outer(np.ones_like(k_s)+S_k,np.ones_like(k_s)+S_k)

        return V

    def _get_covariance(self):
        """Gets the full prior covariance matrix. See equation (22) in |Leclercqetal2019|_.

        Returns
        -------
        S : array, double, dimension=(S,S)
            covariance matrix of the prior

        """
        import numpy as np
        K=self.rbf
        V=self.cv
        S=self.theta_norm**2 * np.multiply(V, K)

        return S

    def _get_inverse_covariance(self):
        """Gets the inverse covariance matrix.

        Returns
        -------
        inv_S : array, double, dimension=(S,S)
            inverse covariance matrix of the prior

        """
        import numpy as np
        from pyselfi.utils import regular_inv

        inv_S = regular_inv(self.covariance, self.EPS_K, self.EPS_residual)

        return inv_S

    def compute(self):
        """Computes the prior (covariance matrix and its inverse).
        """
        self.rbf = self._get_rbf()
        self.cv = self._get_cosmic_variance()
        self.covariance = self._get_covariance()
        self.inv_covariance = self._get_inverse_covariance()

    def logpdf(self,theta,theta_mean,theta_covariance,theta_icov):
        """Returns the log prior probability at a given point in parameter space. See equation (23) in |Leclercqetal2019|_.

        Parameters
        ----------
        theta : array, double, dimension=S
            evaluation point in parameter space
        theta_mean : array, double, dimension=S
            prior mean
        theta_covariance : array, double, dimension=(S,S)
            prior covariance
        theta_icov : array, double, dimension=(S,S)
            inverse prior covariance

        Returns
        -------
        logpdf : double
            log prior probability

        """
        import numpy as np
        diff = theta-theta_mean
        return -diff.dot(theta_icov).dot(diff)/2. - np.linalg.slogdet(2*np.pi*theta_covariance)[1]/2.

    def save(self,fname):
        """Saves the prior to an output file.

        Parameters
        ----------
        fname : :obj:`str`
            output filename

        """
        import h5py as h5
        import numpy as np
        from ctypes import c_double
        from pyselfi.utils import PrintMessage, save_replace_dataset, save_replace_attr

        PrintMessage(3, "Writing prior in data file '{}'...".format(fname))

        with h5.File(fname, 'r+') as hf:
            save_replace_attr(hf, '/prior/theta_norm', self.theta_norm, dtype=c_double)
            save_replace_attr(hf, '/prior/k_corr', self.k_corr, dtype=c_double)
            save_replace_attr(hf, '/prior/alpha_cv', self.alpha_cv, dtype=c_double)
            save_replace_dataset(hf, '/prior/k_s', data=self.k_s, maxshape=(None,), dtype=c_double)
            save_replace_dataset(hf, '/prior/mean', data=self.mean, maxshape=(None,), dtype=c_double)
            save_replace_dataset(hf, '/prior/rbf', data=self.rbf, maxshape=(None, None), dtype=c_double)
            save_replace_dataset(hf, '/prior/cv', data=self.cv, maxshape=(None, None), dtype=c_double)
            save_replace_dataset(hf, '/prior/covariance', data=self.covariance, maxshape=(None, None), dtype=c_double)
            save_replace_dataset(hf, '/prior/inv_covariance', data=self.inv_covariance, maxshape=(None, None), dtype=c_double)

        PrintMessage(3, "Writing prior in data file '{}' done.".format(fname))

    @classmethod
    def load(cls,fname):
        """Loads the prior from an input file.

        Parameters
        ----------
        fname : :obj:`str`
            input filename

        Returns
        -------
        prior : :obj:`prior`
            loaded prior object

        """
        import h5py as h5
        import numpy as np
        from ctypes import c_double
        from pyselfi.utils import PrintMessage
        PrintMessage(3, "Reading prior in data file '{}'...".format(fname))

        with h5.File(fname,'r') as hf:
            mean=np.array(hf.get('/prior/mean'), dtype=c_double)
            k_s = np.array(hf.get('/prior/k_s'), dtype=c_double)
            theta_norm=hf.attrs['/prior/theta_norm']
            k_corr=hf.attrs['/prior/k_corr']
            alpha_cv=hf.attrs['/prior/alpha_cv']
            prior=cls(k_s,mean,theta_norm,k_corr,alpha_cv)
            prior.rbf = np.array(hf.get('/prior/rbf'), dtype=c_double)
            prior.cv = np.array(hf.get('/prior/cv'), dtype=c_double)
            prior.covariance = np.array(hf.get('/prior/covariance'), dtype=c_double)
            prior.inv_covariance = np.array(hf.get('/prior/inv_covariance'), dtype=c_double)

        PrintMessage(3, "Reading prior in data file '{}' done.".format(fname))
        return prior
#end class(prior)
