#!/usr/bin/env python
#-------------------------------------------------------------------------------------
# pySELFI v2.0 -- src/pyselfi/lotkavolterra/selfi.py
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

"""SELFI for inference of the Lotka-Volterra model: main class.

.. _Leclercqetal2019: https://arxiv.org/abs/1902.10149

.. |Leclercqetal2019| replace:: Leclercq *et al*. (2019)
"""

__author__  = "Florent Leclercq"
__version__ = "2.0"
__date__    = "2018-2023"
__license__ = "GPLv3"

from pyselfi import selfi

class lotkavolterra_selfi(selfi):
    """Main class for Lotka-Volterra inference with SELFI. Child of the main class of the SELFI code.

    Attributes
    ----------
    fname : :obj:`str`
        name and address of the selfi output file on disk
    pool_prefix : :obj:`str`
        address and prefix of the simulation pools
    pool_suffix : :obj:`str`
        suffix of the simulation pools
    prior : :obj:`prior`
        the prior, defined as explained in the `SELFI documentation <../../usage/new_model.html>`__
    blackbox : :obj:`blackbox`
        the blackbox simulator, defined as explained in the `SELFI documentation <../../usage/new_blackbox.html>`__
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
    LVsimulator : :obj:`LVsimulator`
        a Lotka-Volterra simulator (instance of the class LVsimulator)
    """
    def __init__(self,fname,pool_prefix,pool_suffix,prior,blackbox,theta_0,Ne,Ns,Delta_theta,phi_obs,LVsimulator):
        """Initializes a selfi object.
        """
        super().__init__(fname,pool_prefix,pool_suffix,prior,blackbox,theta_0,Ne,Ns,Delta_theta,phi_obs)
        self.LVsimulator=LVsimulator

    def Mahalanobis_ms_check(self, theta=None):
        """Computes a quantitative check for model misspecification: the
        Mahalanobis distance between the expansion point (prior mean) :math:`\\theta_\mathrm{0}`
        and a point :math:`\\theta` (usually the posterior mean :math:`\\gamma`), using the prior
        inverse covariance matrix :math:`\\textbf{S}`.

        Parameters
        ----------
        theta : array, double, dimension=S, optional, default=posterior.mean
            a point in parameter space

        Returns
        -------
        dist : double
            Mahalanobis distance between prior and theta/posterior

        """
        import numpy as np
        from scipy.spatial.distance import cdist
        theta = theta or self.posterior.mean
        dist = cdist(np.atleast_2d(self.prior.mean), np.atleast_2d(theta), metric='mahalanobis', VI=self.prior.inv_covariance).item(0)
        return dist

    #################
    # Prior optimizer
    #################

    def loglikelihood_hyperparams(self, theta_fid, alpha_norm, t_smooth, t_chaos):
        """Returns the log pdf of the likelihood for hyperparameters :math:`\\alpha_\mathrm{norm}`, :math:`t_\mathrm{smooth}`, :math:`t_\mathrm{chaos}`, defined as the sum of the log posterior pdf at :math:`\\boldsymbol{\\theta}_\mathrm{fid,i}` for :math:`0 \leq i \leq N_\mathrm{fid}-1`:

        :math:`-2 \\log p(\\alpha_\mathrm{norm}, t_\mathrm{smooth}, t_\mathrm{chaos}) \\equiv \\sum_{i=1}^{N_\mathrm{fid}} [ \\log \left| 2\\pi\\boldsymbol{\Gamma} \\right| + (\\boldsymbol{\\theta}_\mathrm{fid,i}-\\boldsymbol{\\gamma})^\intercal \\boldsymbol{\Gamma}^{-1} (\\boldsymbol{\\theta}_\mathrm{fid,i}-\\boldsymbol{\\gamma}) ]`.

        Parameters
        ----------
        theta_fiducial : array, double, dimension=(N_fid,S)
            set of fiducial points in parameter space to optimize the prior
        alpha_norm : double
            alpha_norm value where to evaluate the likelihood
        t_smooth : double
            t_smooth value where to evaluate the likelihood
        t_chaos : double
            t_chaos value where to evaluate the likelihood

        Returns
        -------
        logpdf : double
            log pdf of the likelihood at the evaluation point

        """
        import numpy as np
        self.prior.set_alpha_norm(alpha_norm)
        self.prior.set_t_smooth(t_smooth)
        self.prior.set_t_chaos(t_chaos)
        self.compute_prior()
        self.compute_posterior()

        return np.mean([self.posterior.logpdf(theta_fid[i], self.posterior.mean, self.posterior.covariance, self.posterior.inv_covariance) for i in range(len(theta_fid))])

    def logposterior_hyperparameters(self, theta_fid, alpha_norm, t_smooth, t_chaos):
        """Returns the log pdf of the posterior for hyperparameters :math:`\\alpha_\mathrm{norm}`, :math:`t_\mathrm{smooth}`, :math:`t_\mathrm{chaos}`. In this case, it is the same as :func:`loglikelihood_hyperparams`, but it would be possible
        add a log prior here.

        Parameters
        ----------
        theta_fiducial : array, double, dimension=(N_fid,S)
            set of fiducial points in parameter space to optimize the prior
        alpha_norm : double
            alpha_norm value where to evaluate the likelihood
        t_smooth : double
            t_smooth value where to evaluate the likelihood
        t_chaos : double
            t_chaos value where to evaluate the likelihood

        Returns
        -------
        logpdf : double
            log pdf of the posterior at the evaluation point

        """
        return self.loglikelihood_hyperparams(theta_fid, alpha_norm, t_smooth, t_chaos)

    def optimize_prior(self, theta_fid, x0=None, alpha_norm_min=1e-10, alpha_norm_max=20,
                       t_smooth_min=1e-10, t_smooth_max=20, t_chaos_min=1e-10, t_chaos_max=20,
                       method='L-BFGS-B', options={'maxiter':50, 'ftol':1e-10, 'gtol':1e-10, 'eps':1e-6, 'disp':False}):
        """Optimizes the SELFI prior. See section II.E. in |Leclercqetal2019|_.

        Parameters
        ----------
        theta_fid : array, double, dimension=S
            fiducial point in parameter space to optimize the prior
        x0 : array, double, dimension=3, optional, default=(prior.get_alpha_norm(),prior.get_t_smooth(),prior.get_t_chaos())
            starting point in parameter space
        alpha_norm_min : double, optional, default=1e-10
            boundary in parameter space
        alpha_norm_max : double, optional, default=20
            boundary in parameter space
        t_smooth_min : double, optional, default=1e-10
            boundary in parameter space
        t_smooth_max : double, optional, default=20
            boundary in parameter space
        t_chaos_min : double, optional, default=1e-10
            boundary in parameter space
        t_chaos_max : double, optional, default=20
            boundary in parameter space
        method : :obj:`str`, *optional, default='L-BFGS-B'*
            optimization method. See documentation of scipy.optimize.minimize
        options : dictionary, optional, default={'maxiter':50, 'ftol':1e-10, 'gtol':1e-10, 'eps':1e-6, 'disp':False}
            optimization options. See documentation of scipy.optimize.minimize

        """
        from scipy.optimize import minimize
        def potential(x):
            import numpy as np
            alpha_norm, t_smooth, t_chaos = x
            return -self.logposterior_hyperparameters(theta_fid, alpha_norm, t_smooth, t_chaos)

        x0=x0 or [self.prior.get_alpha_norm(), self.prior.get_t_smooth(), self.prior.get_t_chaos()]
        res = minimize(potential, x0, method=method,
            options=options, bounds=[[max(1e-10,alpha_norm_min),alpha_norm_max],[max(1e-10,t_smooth_min),t_smooth_max],[max(1e-10,t_chaos_min),t_chaos_max]])
        print(res)

        self.prior.set_alpha_norm(res.x[0])
        self.prior.set_t_smooth(res.x[1])
        self.prior.set_t_chaos(res.x[2])
        self.compute_prior()
        self.compute_posterior()

    ###################
    # Score compression
    ###################

    def compute_grad_f_omega(self, t, t_s, tres, Delta_omega):
        """Computes the gradient of the observations with respect to the input parameters
        :math:`\\boldsymbol{\\omega}`, :math:`\\nabla_\\boldsymbol{\\omega} \\boldsymbol{f}`.
        First, compute :math:`\\nabla_\\boldsymbol{\\omega} \\boldsymbol{\\theta}` using :func:`compute_dtheta_domega` of the :obj:`LVsimulator`, then combine with the usual SELFI blackbox gradient :math:`\\nabla \\boldsymbol{f}` so that :math:`\\nabla_\\boldsymbol{\\omega} \\boldsymbol{f} = \\nabla \\boldsymbol{f} \cdot \\nabla_\\boldsymbol{\\omega} \\boldsymbol{\\theta}`.

        Parameters
        ----------
        t : array, double, dimension=n
            array of time values where we approximate X and Y values
            timestep at each iteration is given by t[n+1] - t[n].
        t_s : array, double, dimension=S
            array of the support time values where theta is defined
        tres : double
            time resolution, defined such that:
            t = np.linspace(0,tmax-1,tmax*tres+1)
            t_s = np.arange(tmax)
        Delta_omega : double
            step size for finite differences

        """
        grad_f = self.likelihood.grad_f

        grad_theta = self.LVsimulator.compute_dtheta_domega(t, t_s, tres, Delta_omega)
        self.grad_theta = grad_theta

        grad_f_omega = grad_theta.T.dot(grad_f.T).T
        self.grad_f_omega = grad_f_omega

    def compute_grad_f_omega_direct(self, t, t_s, tres, Delta_omega, pool_prefix, N):
        """Computes the gradient of the observations with respect to the input parameters
        :math:`\\boldsymbol{\\omega}`, :math:`\\nabla_\\boldsymbol{\\omega} \\boldsymbol{f}`,
        using direct simulations (instead of the SELFI framework).

        Parameters
        ----------
        t : array, double, dimension=n
            array of time values where we approximate X and Y values
            timestep at each iteration is given by t[n+1] - t[n].
        t_s : array, double, dimension=S
            array of the support time values where theta is defined
        tres : double
            time resolution, defined such that:
            t = np.linspace(0,tmax-1,tmax*tres+1)
            t_s = np.arange(tmax)
        Delta_omega : double
            step size for finite differences
        pool_prefix : :obj:`str`
            prefix for the filename of simulation pools
        N : int
            number of realizations required

        """
        import numpy as np
        omega0 = self.LVsimulator.omega
        P = self.likelihood.blackbox.P
        grad_f_omega = np.zeros((P,4))
        theta_0 = self.LVsimulator.compute_theta(omega0, t, t_s, tres)
        pool_0 = self.likelihood.blackbox.compute_pool(theta_0, 0, pool_prefix+"0.h5", N)
        f0 = pool_0.Phi.mean(axis=0)

        for n in range(4):
            omega_p1, omega_p2, omega_p3, omega_m1, omega_m2, omega_m3 = np.copy(omega0), np.copy(omega0), np.copy(omega0), np.copy(omega0), np.copy(omega0), np.copy(omega0)
            omega_p1[n] += 1*Delta_omega
            omega_p2[n] += 2*Delta_omega
            omega_p3[n] += 3*Delta_omega
            omega_m1[n] -= 1*Delta_omega
            omega_m2[n] -= 2*Delta_omega
            omega_m3[n] -= 3*Delta_omega
            theta_p1 = self.LVsimulator.compute_theta(omega_p1, t, t_s, tres)
            theta_p2 = self.LVsimulator.compute_theta(omega_p2, t, t_s, tres)
            theta_p3 = self.LVsimulator.compute_theta(omega_p3, t, t_s, tres)
            theta_m1 = self.LVsimulator.compute_theta(omega_m1, t, t_s, tres)
            theta_m2 = self.LVsimulator.compute_theta(omega_m2, t, t_s, tres)
            theta_m3 = self.LVsimulator.compute_theta(omega_m3, t, t_s, tres)
            pool_p1 = self.likelihood.blackbox.compute_pool(theta_p1, +1, pool_prefix+"omega"+str(n+1)+"_p1.h5", N)
            pool_p2 = self.likelihood.blackbox.compute_pool(theta_p2, +2, pool_prefix+"omega"+str(n+1)+"_p2.h5", N)
            pool_p3 = self.likelihood.blackbox.compute_pool(theta_p3, +3, pool_prefix+"omega"+str(n+1)+"_p3.h5", N)
            pool_m1 = self.likelihood.blackbox.compute_pool(theta_m1, -1, pool_prefix+"omega"+str(n+1)+"_m1.h5", N)
            pool_m2 = self.likelihood.blackbox.compute_pool(theta_m2, -2, pool_prefix+"omega"+str(n+1)+"_m2.h5", N)
            pool_m3 = self.likelihood.blackbox.compute_pool(theta_m3, -3, pool_prefix+"omega"+str(n+1)+"_m3.h5", N)
            fp1 = pool_p1.Phi.mean(axis=0)
            fp2 = pool_p2.Phi.mean(axis=0)
            fp3 = pool_p3.Phi.mean(axis=0)
            fm1 = pool_m1.Phi.mean(axis=0)
            fm2 = pool_m2.Phi.mean(axis=0)
            fm3 = pool_m3.Phi.mean(axis=0)

            #grad_f_omega[:,n]=(fp1 - f0)/Delta_omega
            #grad_f_omega[:,n]=(-1/2.*fm1 + 1/2.*fp1)/Delta_omega
            #grad_f_omega[:,n]=(1/12.*fm2 - 2/3.*fm1 + 2/3.*fp1 - 1/12.*fp2)/Delta_omega
            grad_f_omega[:,n]=(-1/60.*fm3 + 3/20.*fm2 - 3/4.*fm1 + 3/4.*fp1 - 3/20.*fp2 + 1/60.*fp3)/Delta_omega

        return grad_f_omega

    def compute_Fisher_matrix(self):
        """Computes the Fisher matrix of the problem, once all necessary quantities are here.
        """
        from pyselfi.utils import regular_inv
        grad_f_omega = self.grad_f_omega
        f0 = self.likelihood.f_0
        self.f0 = f0

        C_0 = self.likelihood.C_0
        self.C_0 = C_0

        F_0 = grad_f_omega.T.dot(C_0).dot(grad_f_omega)
        self.F_0 = F_0

        invF_0 = regular_inv(F_0)
        self.invF_0 = invF_0

        projPhiOmega = invF_0.dot(grad_f_omega.T).dot(C_0)
        self.projPhiOmega = projPhiOmega

    def score_compression(self, phi):
        """Compress a data vector :math:`\\boldsymbol{\\Phi}` using score compression (the optimal quasi-maximum
        likelihood estimator).

        Parameters
        ----------
        phi : array, double, dimension=P
            input data vector to be compressed

        Returns
        -------
        omega_tilde : array, double, dimension=4
            output compressed summaries
        """
        omega0 = self.LVsimulator.omega
        f0 = self.f0
        projPhiOmega = self.projPhiOmega
        omega_tilde = omega0 + projPhiOmega.dot(phi-f0)
        return omega_tilde

    def FisherRao_distance(self, omega_tilde1, omega_tilde2):
        """Computes the Fisher-Rao distance between two compressed summaries
        :math:`\\tilde{\\boldsymbol{\\omega}}_1` and :math:`\\tilde{\\boldsymbol{\\omega}}_2`.

        Parameters
        ----------
        omega_tilde1 : array, double, dimension=4
            first vector of compressed summaries

        omega_tilde2 : array, double, dimension=4
            second vector of compressed summaries

        Returns
        -------
        dist : double
            Fisher-Rao distance between omega_tilde1 and
            omega_tilde2, using the Fisher matrix of the problem
        """
        import numpy as np
        from scipy.spatial.distance import cdist
        #dist = np.sqrt((omega_tilde1-omega_tilde2).T.dot(self.F_0).dot(omega_tilde1-omega_tilde2))
        dist = cdist(np.atleast_2d(omega_tilde1), np.atleast_2d(omega_tilde2), metric='mahalanobis', VI=self.F_0).item(0)
        return dist

    ###########
    # Utilities
    ###########

    def save(self):
        """Saves the data compression information to the SELFI output file.

        """
        import h5py as h5
        from ctypes import c_double
        from pyselfi.utils import PrintMessage, save_replace_dataset
        fname=self.fname
        PrintMessage(3, "Writing data compression information in data file '{}'...".format(fname))

        with h5.File(fname, 'r+') as hf:
            save_replace_dataset(hf, '/compression/grad_f_omega', data=self.grad_f_omega, maxshape=(None, None), dtype=c_double)
            save_replace_dataset(hf, '/compression/f0', data=self.f0, maxshape=(None,), dtype=c_double)
            save_replace_dataset(hf, '/compression/C_0', data=self.C_0, maxshape=(None, None), dtype=c_double)
            save_replace_dataset(hf, '/compression/F_0', data=self.F_0, maxshape=(None, None), dtype=c_double)
            save_replace_dataset(hf, '/compression/invF_0', data=self.invF_0, maxshape=(None, None), dtype=c_double)
            save_replace_dataset(hf, '/compression/projPhiOmega', data=self.projPhiOmega, maxshape=(None, None), dtype=c_double)

        PrintMessage(3, "Writing data compression information in data file '{}' done.".format(fname))

    def load(self,fname):
        """Loads the data compression information from the SELFI input file.

        Parameters
        ----------
        fname : :obj:`str`
            input filename

        """
        import h5py as h5
        import numpy as np
        from ctypes import c_double
        from pyselfi.utils import PrintMessage
        PrintMessage(3, "Reading data compression information in data file '{}'...".format(fname))

        with h5.File(fname,'r') as hf:
            self.grad_f_omega=np.array(hf.get('/compression/grad_f_omega'), dtype=c_double)
            self.f0=np.array(hf.get('/compression/f0'), dtype=c_double)
            self.C_0=np.array(hf.get('/compression/C_0'), dtype=c_double)
            self.F_0=np.array(hf.get('/compression/F_0'), dtype=c_double)
            self.invF_0=np.array(hf.get('/compression/invF_0'), dtype=c_double)
            self.projPhiOmega=np.array(hf.get('/compression/projPhiOmega'), dtype=c_double)

        PrintMessage(3, "Reading data compression information in data file '{}' done.".format(fname))


# end class(child_selfi)
