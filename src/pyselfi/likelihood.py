#!/usr/bin/env python
#-------------------------------------------------------------------------------------
# pySELFI v2.0 -- src/pyselfi/likelihood.py
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

"""Routines related to the SELFI effective likelihood.

.. _Leclercqetal2019: https://arxiv.org/abs/1902.10149

.. |Leclercqetal2019| replace:: Leclercq *et al*. (2019)
"""

__author__  = "Florent Leclercq"
__version__ = "2.0"
__date__    = "2018-2023"
__license__ = "GPLv3"

class likelihood(object):
    """This class represents the SELFI likelihood. See equations (11)-(15) and (17)-(18) in |Leclercqetal2019|_ for expressions.

    Attributes
    ----------
    blackbox : :obj:`blackbox`
        the blackbox simulator, defined as explained in the `SELFI documentation <../usage/new_blackbox.html>`__
    theta_0 : array, double, dimension=S
        expansion point in parameter space
    Ne : int
        number of simulations at the expansion point, to compute the covariance matrix.
    Ns : int
        number of simulations per expansion direction in parameter space. (Ne may be larger than Ns, in which case only the first Ns simulations at the expansion point are used to compute the gradient)
    Delta_theta : double
        step size for finite differencing to compute the blackbox gradient in parameter space

    """

    # Initialization
    def __init__(self,blackbox,theta_0,Ne,Ns,Delta_theta):
        """Initializes a likelihood object.
        """
        self.blackbox=blackbox
        self.theta_0=theta_0
        self.S=theta_0.shape[0]
        self.P=blackbox.P
        self.Ne=Ne
        self.Ns=Ns
        self.Delta_theta=Delta_theta
        self.EPS_Sigma=1e-7
        self.EPS_residual=1e-3

    def _get_perturb_theta(self,d,h):
        """Computes the perturbed theta = theta_0 + Delta_theta. See equation (18) in |Leclercqetal2019|_.

        Parameters
        ----------
        d : int
            index giving the direction in parameter space, from 0 to S
        h : double
            step size in parameter space

        Returns
        -------
        theta : double
            theta_0 if d=0, or theta_0 + Delta_theta if d is from 1 to S

        """
        import numpy as np
        theta_0=self.theta_0
        if d==0:
            return theta_0
        else:
            Delta_theta = np.zeros_like(theta_0)
            Delta_theta[d-1] = h
            return theta_0 + Delta_theta

    def _compute_Phi_d(self,pool_fname,d,N=None,h=None):
        """Computes a pool of blackbox simulations at one of the required points in parameter space.

        Parameters
        ----------
        pool_fname : :obj:`str`
            filename of the pool
        d : int
            direction in parameter space, from 0 to S
        N : int, optional, default=class value for Ne if d=0, or Ns otherwise
            number of simulations per expansion direction in parameter space
        h : double, optional, default=class value
            step size for finite differencing to compute the blackbox gradient in parameter space

        """
        N=N or (self.Ne if d == 0 else self.Ns)
        h=h or self.Delta_theta
        theta_d = self._get_perturb_theta(d,h)
        Phi_d = self.blackbox.compute_pool(theta_d, d, pool_fname, N)
        setattr(self, 'Phi_'+str(d), Phi_d)

    @classmethod
    def _covariance(cls,Phi):
        """Returns the estimated covariance of a set of summaries, with the prefactor for a virtual signal field. See equations (13)-(14) in |Leclercqetal2019|_.

        Parameters
        ----------
        Phi : array, double, dimension=(N,P)
            a set of N simulations

        Returns
        -------
        Sigma^hat_theta' : array, double, dimension=(P,P)
            (N+1)/N * Sigma^hat_theta where Sigma^hat_theta is the unbiased covariance estimator

        """
        import numpy as np
        from ctypes import c_double
        # Setup the correction factor for the evaluated covariance matrix
        N=Phi.shape[0]
        covariance_prefactor=(N+1)/N

        # Return the evaluated covariance (using the unbiased estimator)
        return covariance_prefactor * np.cov(Phi.T,ddof=1).astype(c_double)

    def compute_C_0(self):
        """Computes the estimated covariance of a set of summaries at the expansion point, C_0, and log|2piC_0|. Assumes that Phi_0 has already been computed.
        """
        import numpy as np
        self.C_0 = self._covariance(self.Phi_0.Phi)
        self.logdet2piC_0 = np.linalg.slogdet(2*np.pi*self.C_0)[1]

    @classmethod
    def _inverse_covariance(cls,Phi,covariance,EPS_Sigma=1e-7,EPS_residual=1e-3):
        """Computes the unbiased inverse of an estimated covariance matrix, using the `Hartlap, Simon & Schneider (2007) <https://arxiv.org/abs/astro-ph/0608064>`__ correction factor. See equations (13) and (15) in |Leclercqetal2019|_.

        Parameters
        ----------
        Phi : array, double, dimension=(N,P)
            a set of N simulations
        covariance : array, double, dimension=(P,P)
            estimated covariance matrix
        EPS_Sigma : double, optional, default=1e-7
            epsilon to be added to the diagonal of the covariance matrix, if necessary
        EPS_residual : double, optional, default=1e-3
            maximum residual in covariance*covariance^{-1} before attempting regularization

        Returns
        -------
        Sigma^hat_theta^{-1}' : array, double, dimension=(P,P)
            alpha*(Sigma^hat_theta')^{-1} where (Sigma^hat_theta')^{-1} is the inverse sample covariance

        """
        from pyselfi.utils import regular_inv
        # Setup the Hartlap, Simon & Schneider (2007) correction factor
        N=Phi.shape[0]
        P=Phi.shape[1]
        alpha = (N-P-2.)/(N-1.) if P<N-2 else 1.

        # Compute the inverse sample covariance
        inv_covariance = regular_inv(covariance, EPS_Sigma, EPS_residual)

        return alpha*inv_covariance

    def compute_inv_C_0(self):
        """Computes the estimated inverse covariance of a set of summaries at the expansion point, inv_C_0. Assumes that Phi_0 has already been computed.
        """
        self.inv_C_0 = self._inverse_covariance(self.Phi_0.Phi, self.C_0, self.EPS_Sigma, self.EPS_residual)

    @classmethod
    def _average(cls,Phi,N):
        """Returns the average of the first N of a set of summaries.

        Parameters
        ----------
        Phi : array, double, dimension=(NN,P)
            a set of NN simulations
        N : int
            take the average of the first N simulations

        """
        import numpy as np
        return np.mean(Phi[:N],axis=0)

    def _compute_f_d(self,d,Ns):
        """Computes the average blackbox at one of the required points in parameter space. Assumes Phi_d has already been computed.

        Parameters
        ----------
        d : int
            direction in parameter space, from 0 to S
        Ns : int
            number of simulations to be used

        """
        f_d = self._average(getattr(self, 'Phi_'+str(d)).Phi, Ns)
        setattr(self, 'f_'+str(d), f_d)

    def _compute_grad_d(self,d,h):
        """Computes the gradient in one direction of parameter space, see equation (18) in |Leclercqetal2019|_. Assumes that f_0 and f_d have already been computed.

        Parameters
        ----------
        d : int
            direction in parameter space, from 1 to S
        h : double
            step size for finite differencing

        """
        f_0 = self.f_0
        f_d = getattr(self, 'f_'+str(d))
        return (f_d-f_0)/h

    def _gradient(self,Ns=None,h=None):
        """Returns the gradient of the linearised blackbox.

        Parameters
        ----------
        Ns : int, optional, default=class value
            number of simulations per expansion direction in parameter space
        h : double, optional, default=class value
            step size for finite differencing to compute the blackbox gradient in parameter space

        Returns
        -------
        grad_f : array, double, dimension=(P,S)
            gradient of the linearised blackbox

        """
        import numpy as np
        from ctypes import c_double

        Ns=Ns or self.Ns
        h=h or self.Delta_theta
        grad_f = np.zeros((self.S, self.P), dtype=c_double)

        # Compute the average blackbox along all directions
        for d in range(self.S+1):
            self._compute_f_d(d,Ns)

        # Evaluate the gradients along all directions
        for d in range(1,self.S+1):
            grad_f[d-1] = self._compute_grad_d(d,h)

        grad_f = grad_f.T
        return grad_f

    def compute_gradient(self,Ns=None,h=None):
        """Computes the gradient of the linearised blackbox.

        Parameters
        ----------
        Ns : int, optional, default=class value
            number of simulations per expansion direction in parameter space
        h : double, optional, default=class value
            step size for finite differencing to compute the blackbox gradient in parameter space

        """
        self.grad_f = self._gradient(Ns,h)

    def run_simulations(self, pool_prefix, pool_suffix, d=None, Ne=None, Ns=None, h=None):
        """Runs some or all the necessary simulations to compute the likelihood. All the necessary simulations amount to Ne + Ns*S blackbox evaluations.

        Parameters
        ----------
        pool_prefix : :obj:`str`
            address and prefix of the filenames for simulation pools
        pool_suffix : :obj:`str`
            suffix of the filenames for simulation pools
        d : int or array of int, optional, default=None
            directions in parameter space, from 0 to S. If set to None, all simulations are run or loaded
        Ne : int, optional, default=class value
            number of simulations at the expansion point, to compute the covariance matrix
        Ns : int, optional, default=class value
            number of simulations per expansion direction in parameter space. (Ne may be larger than Ns, in which case only the first Ns simulations at the expansion point are used to compute the gradient)
        h : double, optional, default=class value
            step size for finite differencing to compute the blackbox gradient in parameter space

        """
        # define or get the required directions to be treated
        import numpy as np
        if d is None:
            d=list(np.arange(0,self.S+1))
        else:
            d=list(np.atleast_1d(d))

        # run simulations at the expansion point
        if 0 in d:
            pool_fname=pool_prefix+"0"+pool_suffix
            self._compute_Phi_d(pool_fname, 0, Ne, h)
            d.remove(0)

        # run simulations along all directions
        for this_d in d:
            pool_fname=pool_prefix+str(this_d)+pool_suffix
            self._compute_Phi_d(pool_fname, this_d, Ns, h)

    def compute(self,Ns=None,h=None):
        """Computes all the necessary quantities to characterize the likelihood (estimated covariance matrix, its inverse, and the gradient).

        Parameters
        ----------
        Ns : int, optional, default=class value
            number of simulations per expansion direction in parameter space
        h : double, optional, default=class value
            step size for finite differencing to compute the blackbox gradient in parameter space

        """
        self.compute_C_0()
        self.compute_inv_C_0()
        self.compute_gradient(Ns,h)

    def data_model(self,theta):
        """Evaluates the linearised data model at a given point in parameter space, see equation (17) in |Leclercqetal2019|_.

        Parameters
        ----------
        theta : array, double, dimension=S
            evaluation point in parameter space

        """
        f_0 = self.f_0
        grad_f = self.grad_f
        theta_0 = self.theta_0
        f_theta = f_0 + grad_f.dot(theta-theta_0)
        return f_theta

    def logpdf(self,theta,phi_obs):
        """Returns the log likelihood probability at a given point in parameter space, see equation (19) in |Leclercqetal2019|_.

        Parameters
        ----------
        theta : array, double, dimension=S
            evaluation point in parameter space
        phi_obs : array, double, dimension=P
            observed summaries vector

        Returns
        -------
        logpdf : double
            log likelihood probability

        """
        phi_pred = self.data_model(theta)
        diff = phi_obs-phi_pred
        inv_C_0 = self.inv_C_0
        logdet2piC_0 = self.logdet2piC_0
        return -diff.dot(inv_C_0).dot(diff)/2. - logdet2piC_0/2.

    def save(self,fname):
        """Saves the likelihood to an output file.

        Parameters
        ----------
        fname : :obj:`str`
            output filename

        """
        import h5py as h5
        from ctypes import c_double, c_double
        from pyselfi.utils import PrintMessage, save_replace_dataset, save_replace_attr
        PrintMessage(3, "Writing likelihood in data file '{}'...".format(fname))

        with h5.File(fname, 'r+') as hf:
            for d in range(self.S+1):
                save_replace_dataset(hf, '/likelihood/f_'+str(d), data=getattr(self,'f_'+str(d)).astype(c_double), maxshape=(None,), dtype=c_double)
            save_replace_dataset(hf, '/likelihood/C_0', data=self.C_0.astype(c_double), maxshape=(None, None), dtype=c_double)
            save_replace_attr(hf, '/likelihood/logdet2piC_0', self.logdet2piC_0, dtype=c_double)
            save_replace_dataset(hf, '/likelihood/inv_C_0', data=self.inv_C_0.astype(c_double), maxshape=(None, None), dtype=c_double)
            save_replace_dataset(hf, '/likelihood/grad_f', data=self.grad_f.astype(c_double), maxshape=(None, None), dtype=c_double)

        PrintMessage(3, "Writing likelihood in data file '{}' done.".format(fname))

    def load(self,fname):
        """Loads the likelihood from an input file.

        Parameters
        ----------
        fname : :obj:`str`
            input filename

        """
        import numpy as np
        import h5py as h5
        from ctypes import c_double
        from pyselfi.utils import PrintMessage
        PrintMessage(3, "Reading likelihood in data file '{}'...".format(fname))

        with h5.File(fname, 'r') as hf:
            for d in range(self.S+1):
                setattr(self, 'f_'+str(d), np.array(hf.get('/likelihood/f_'+str(d)), dtype=c_double).astype(c_double))
            self.C_0 = np.array(hf.get('/likelihood/C_0'), dtype=c_double).astype(c_double)
            self.logdet2piC_0 = hf.attrs['/likelihood/logdet2piC_0']
            self.inv_C_0 = np.array(hf.get('/likelihood/inv_C_0'), dtype=c_double).astype(c_double)
            self.grad_f = np.array(hf.get('/likelihood/grad_f'), dtype=c_double).astype(c_double)

        PrintMessage(3, "Reading likelihood from in file '{}' done.".format(fname))
#end class(likelihood)
