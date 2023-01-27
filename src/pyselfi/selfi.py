#!/usr/bin/env python
#-------------------------------------------------------------------------------------
# pySELFI v2.0 -- src/pyselfi/selfi.py
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

"""Main class of the SELFI code.
"""

__author__  = "Florent Leclercq"
__version__ = "2.0"
__date__    = "2018-2023"
__license__ = "GPLv3"

class selfi(object):
    """Main class of the SELFI code.

    Attributes
    ----------
    fname : :obj:`str`
        name and address of the selfi output file on disk
    pool_prefix : :obj:`str`
        address and prefix of the simulation pools
    pool_suffix : :obj:`str`
        suffix of the simulation pools
    prior : :obj:`prior`
        the prior, defined as explained in the `SELFI documentation <../usage/new_model.html>`__
    blackbox : :obj:`blackbox`
        the blackbox simulator, defined as explained in the `SELFI documentation <../usage/new_blackbox.html>`__
    theta_0 : array, double, dimension=S
        expansion point in parameter space
    Ne : int
        number of simulations at the expansion point, to compute the covariance matrix
    Ns : int
        number of simulations per expansion direction in parameter space. (Ne may be larger than Ns, in which case only the first Ns simulations at the expansion point are used to compute the gradient)
    Delta_theta : double
        step size for finite differencing to compute the blackbox gradient in parameter space
    phi_obs : array, double, dimension=P
        the vector of observed data

    """

    # Initialization
    def __init__(self,fname,pool_prefix,pool_suffix,prior,blackbox,theta_0,Ne,Ns,Delta_theta,phi_obs):
        """Initializes a selfi object.
        """
        import h5py as h5
        from pathlib import Path

        # Output files
        self.fname=fname
        f = Path(fname)
        if not f.is_file():
            h5.File(fname,"w")
        self.pool_prefix=pool_prefix
        self.pool_suffix=pool_suffix

        # Observed data
        self.phi_obs=phi_obs

        # Prior setup
        self.prior=prior

        # Blackbox setup
        from pyselfi.likelihood import likelihood
        self.likelihood=likelihood(blackbox,theta_0,Ne,Ns,Delta_theta)

    ##############################
    # Methods related to the prior
    ##############################
    def compute_prior(self):
        """Computes the prior.
        """
        self.prior.compute()

    def save_prior(self):
        """Saves the prior to the SELFI output file.
        """
        self.prior.save(self.fname)

    def load_prior(self):
        """Loads the prior from the SELFI output file.
        """
        self.prior.load(self.fname)

    ###################################
    # Methods related to the likelihood
    ###################################
    def run_simulations(self, d=None, pool_prefix=None, pool_suffix=None, Ne=None, Ns=None, h=None):
        """Runs the necessary simulations for the likelihood.

        Parameters
        ----------
        d : int or array of int, optional, default=None
            directions in parameter space, from 0 to S. If set to None, all simulations are run or loaded
        pool_prefix : :obj:`str`, optional, default=class value
            address and prefix of the filenames for simulation pools
        pool_suffix : :obj:`str`, optional, default=class value
            suffix of the filenames for simulation pools
        Ne : int, optional, default=class value
            number of simulations at the expansion point, to compute the covariance matrix.
        Ns : int, optional, default=class value
            number of simulations per expansion direction in parameter space
        h : double, optional, default=class value
            step size for finite differencing to compute the blackbox gradient in parameter space

        """
        pool_prefix=pool_prefix or self.pool_prefix
        pool_suffix=pool_suffix or self.pool_suffix
        self.likelihood.run_simulations(pool_prefix,pool_suffix,d,Ne,Ns,h)

    def compute_likelihood(self,Ns=None,h=None):
        """Computes the likelihood, assuming that the necessary simulations are available.

        Parameters
        ----------
        Ns : int, optional, default=class value
            number of simulations per expansion direction in parameter space
        h : double, optional, default=class value
            step size for finite differencing to compute the blackbox gradient in parameter space

        """
        self.likelihood.compute(Ns,h)

    def save_likelihood(self):
        """Saves the likelihood to the SELFI output file.
        """
        self.likelihood.save(self.fname)

    def load_likelihood(self):
        """Loads the likelihood from the SELFI output file.
        """
        self.likelihood.load(self.fname)

    ##################################
    # Methods related to the posterior
    ##################################
    def compute_posterior(self):
        """Computes the SELFI posterior, assuming that the prior and likelihood have been computed.
        """
        from pyselfi.posterior import posterior
        self.posterior=posterior(self.prior.mean,self.prior.covariance,self.prior.inv_covariance,self.likelihood.f_0,self.likelihood.C_0,self.likelihood.inv_C_0,self.likelihood.grad_f,self.phi_obs)
        self.posterior.compute()

    def save_posterior(self):
        """Saves the posterior to the SELFI output file.
        """
        self.posterior.save(self.fname)

    def load_posterior(self):
        """Loads the posterior from the SELFI output file.
        """
        from pyselfi.posterior import posterior
        self.posterior=posterior.load(self.fname)
#end class(selfi)
