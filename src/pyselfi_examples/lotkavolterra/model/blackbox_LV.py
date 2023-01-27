#!/usr/bin/env python
#-------------------------------------------------------------------------------------
# pySELFI v2.0 -- src/pyselfi_examples/lotkavolterra/model/blackbox_LV.py
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

"""A blackbox to simulate prey and predator populations and observations from the Lotka-Volterra model. See section III in |Leclercq2022|_ for a description of the Lotka-Volterra data model.

.. _Leclercq2022: https://arxiv.org/abs/2209.11057

.. |Leclercq2022| replace:: Leclercq (2022)
"""

__author__  = "Florent Leclercq"
__version__ = "2.0"
__date__    = "2018-2023"
__license__ = "GPLv3"


class blackbox(object):
    """This class represents a SELFI blackbox.

    Attributes
    ----------
    P : int
        number of output summary statistics. This is the only mandatory argument
    X0 : double
        initial condition: number of prey
    Y0 : double
        initial condition: number of predators
    seed : int
        (first) random seed to be used for the realizations
    fixnoise : bool
        fix the noise realization? if True the seed used will always be seed
    model : int, optional, default=0
        0= correct data model; 1=misspecified data model
    threshold : double
        maximum number of individuals (prey or predators) that can be observed
    mask : array, double, dimension=n
        a mask containing zeros and ones
    Xefficiency : array, double, dimension=len(t)
        detection efficiency of prey as a function of time
    Yefficiency : array, double, dimension=len(t)
        detection efficiency of prey as a function of time
    obsP : double
        rate of prey misses due to correlation between prey and predators
    obsQ : double
        rate of predator misses due to correlation between prey and predators
    obsR : double
        strength of demographic noise
    obsS : double
        overall strength of observational noise
    obsT : double
        strength of non-diagonal term in observational noise

    """

    # Initialization
    def __init__(self, P, X0, Y0, seed, fixnoise, **metaparams):
        """Initializes the blackbox object.
        """
        self.P=P
        self.X0=X0
        self.Y0=Y0
        self.seed=seed
        self.fixnoise=fixnoise
        self.metaparams = metaparams

    def Xobs(self, S, mask):
        """Evalutes the timesteps at which prey observations are defined.

        Parameters
        ----------
        S : int
            number of input parameters
        mask : array, double, dimension=S
            a mask containing zeros and ones

        Returns
        -------
        Xobs : array, double, dimension=PX
            array containing timesteps, such that PX=np.sum(mask)

        """
        import numpy as np
        return np.arange(S//2)[np.where(mask==1)]

    def Yobs(self, S):
        """Evalutes the timesteps at which predator observations are defined.

        Parameters
        ----------
        S : int
            number of input parameters

        Returns
        -------
        Yobs : array, double, dimension=S/2
            array containing timesteps (np.arange(S//2))

        """
        import numpy as np
        return np.arange(S//2)

    def make_signal(self, theta_true):
        """Evaluates the blackbox to make mock unobserved signal.

        Parameters
        ----------
        theta_true : array, double, dimension=S
            vector of true input function

        Returns
        -------
        theta_signal : array, double, dimension=S
            vector of unobserved signal

        """
        import numpy as np
        from lotkavolterra_simulator import LVobserver

        X0 = self.X0
        Y0 = self.Y0
        S = len(theta_true)
        Xtrue = theta_true[0:S//2]
        Ytrue = theta_true[S//2:S]
        LVo = LVobserver(Xtrue, Ytrue, X0, Y0)
        Xtsignal, Ytsignal = LVo.make_signal(**self.metaparams)
        theta_signal = np.concatenate((Xtsignal, Ytsignal))
        return theta_signal

    def make_demographic_noise_cov(self, theta_true):
        """Evaluates the demographic noise covariance matrix.

        Parameters
        ----------
        theta_true : array, double, dimension=S
            vector of true input function

        Returns
        -------
        D : array, double, dimension=(S,S)
            demographic noise covariance matrix

        """
        from lotkavolterra_simulator import LVobserver

        X0 = self.X0
        Y0 = self.Y0
        S = len(theta_true)
        Xtrue = theta_true[0:S//2]
        Ytrue = theta_true[S//2:S]
        LVo = LVobserver(Xtrue, Ytrue, X0, Y0)
        return LVo.Dnoise_cov(**self.metaparams)

    def make_observational_noise_cov(self, theta_true):
        """Evaluates the observational noise covariance matrix.

        Parameters
        ----------
        theta_true : array, double, dimension=S
            vector of true input function

        Returns
        -------
        O : array, double, dimension=(S,S)
            observational noise covariance matrix

        """
        from lotkavolterra_simulator import LVobserver

        X0 = self.X0
        Y0 = self.Y0
        S = len(theta_true)
        Xtrue = theta_true[0:S//2]
        Ytrue = theta_true[S//2:S]
        LVo = LVobserver(Xtrue, Ytrue, X0, Y0)
        return LVo.Onoise_cov(Xtrue, Ytrue, **self.metaparams)

    def make_data(self, theta_true, seed=None):
        """Evaluates the blackbox to make mock data.

        Parameters
        ----------
        theta_true : array, double, dimension=S
            vector of true input parameters
        seed : int, optional, default=None
            value of the seed to generate the realization, uses current state if set to None

        Returns
        -------
        Phi : array, double, dimension=P
            vector of summary statistics

        """
        return self.evaluate(theta_true, seed=seed)

    def evaluate(self, theta, seed=None, i=0, N=0):
        """Evaluates the blackbox from an input vector of parameters.

        Parameters
        ----------
        theta : array, double, dimension=S
            vector of input parameters
        seed : int, optional, default=None
            value of the seed to generate the realization, uses current state if set to None
        i : int, optional, default=0
            current evaluation index of the blackbox
        N : int, optional, default=0
            total number of evaluations of the blackbox

        Returns
        -------
        Phi : array, double, dimension=P
            vector of summary statistics

        """
        import numpy as np
        from pyselfi.utils import PrintMessage, PrintValue
        from lotkavolterra_simulator import LVobserver

        if seed is not None:
            PrintValue("seed", seed)
            np.random.seed(seed)

        X0 = self.X0
        Y0 = self.Y0
        S = len(theta)
        Xtrue = theta[0:S//2]
        Ytrue = theta[S//2:S]
        LVo = LVobserver(Xtrue, Ytrue, X0, Y0)
        LVo.observe(**self.metaparams)
        Xdata, Ydata = LVo.Xdata, LVo.Ydata
        Phi = np.concatenate((Xdata, Ydata))

        return Phi

    def compute_pool(self, theta, d, pool_fname, N):
        """Computes a pool of realizations of the blackbox.

        Parameters
        ----------
        theta : array, double, dimension=S
            vector of input parameters
        d : int
            direction in parameter space, from 0 to S
        pool_fname : :obj:`str`
            pool file name
        N : int
            number of realizations required

        Returns
        -------
        p : :obj:`pool`
            simulation pool

        """
        def get_current_seed(seed,fixnoise,i):
            if seed is None:
                this_seed=None
            else:
                if fixnoise:
                    this_seed=seed
                else:
                    this_seed=seed+i

            return this_seed

        from pyselfi.pool import pool
        p=pool(pool_fname,N)

        # Run N evaluations of the blackbox at the desired point
        while not p.finished:
            i = p.N_sims+1
            this_seed = get_current_seed(self.seed,self.fixnoise,i)
            Phi = self.evaluate(theta,this_seed,i,N)
            p.add_sim(Phi)
            # It can be convenient to define a variable save_frequency, otherwise one can just call p.save() after every blackbox evaluation
            if hasattr(self,"save_frequency") and i%self.save_frequency==0:
                p.save()
        p.save()

        return p
#end class(blackbox)
