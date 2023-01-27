#!/usr/bin/env python
#-------------------------------------------------------------------------------------
# pySELFI v2.0 -- src/pyselfi/lotkavolterra/prior.py
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

"""Routines related to the SELFI prior for time-dependent smooth functions.

.. _Leclercq2022: https://arxiv.org/abs/2209.11057

.. |Leclercq2022| replace:: Leclercq (2022)
"""

__author__  = "Florent Leclercq"
__version__ = "2.0"
__date__    = "2018-2023"
__license__ = "GPLv3"

class lotkavolterra_prior(object):
    """This class represents the SELFI prior for time-dependent smooth functions.

    The prior is a Gaussian with mean :math:`\\boldsymbol{\\theta}_0` and covariance matrix :math:`\\textbf{S} = \\alpha_\mathrm{norm}^2 \\textbf{V} \\circ \\textbf{K}` (Hadamard product), such that:
    :math:`\\textbf{K} = \\begin{pmatrix} \\textbf{K}_x & \\textbf{0} \\\\ \\textbf{0} & \\textbf{K}_y \\end{pmatrix}`, :math:`(\\textbf{K}_z)_{ij} = -\\dfrac{1}{2} \\left( \\dfrac{t_i-t_j}{t_\mathrm{smooth}} \\right)^2`, :math:`\\textbf{V} = \\begin{pmatrix} x_0 \\textbf{u}\\textbf{u}^\\intercal & \\textbf{0} \\\\ \\textbf{0} & y_0 \\textbf{u}\\textbf{u}^\\intercal \\end{pmatrix}`, :math:`(\\textbf{u})_i = 1 + \\dfrac{t_i}{t_\mathrm{chaos}}`.


    Attributes
    ----------
    t_s : array, double, dimension=S/2
        array of timesteps
    theta_0 : array, double, dimension=S
        prior mean and blackbox expansion point
    X0 : double
        initial condition (number of prey at t=0)
    Y0 : double
        initial condition (number of predators at t=0)
    alpha_norm : array, double, dimension=2
        overall amplitude of the prior covariance matrix (prey, predators)
    t_smooth : array, double, dimension=2
        typical time parameter for the smoothness prior (prey, predators)
    t_chaos : array, double, dimension=2
        typical time parameter for the time-dependent uncertainty (prey, predators)

    """

    # Initialization
    def __init__(self,t_s,X0,Y0,theta_0,alpha_norm,t_smooth,t_chaos):
        """Initializes a prior object.
        """
        self.t_s=t_s
        self.mean=theta_0
        self.X0=X0
        self.Y0=Y0
        self.set_alpha_norm(alpha_norm)
        self.set_t_smooth(t_smooth)
        self.set_t_chaos(t_chaos)
        self.EPS_K=1e-7
        self.EPS_residual=1e-3

    def set_alpha_norm(self, alpha_norm):
        """Sets alpha_norm as a 2-vector: :math:`\\alpha_\mathrm{norm} [x_0,y_0]`.
        """
        import numpy as np
        self.alpha_norm=alpha_norm*np.array([self.X0,self.Y0])

    def get_alpha_norm(self):
        """Gets alpha_norm (overall amplitude): :math:`\\alpha_\mathrm{norm}[0]/x_0`.
        """
        return self.alpha_norm[0]/self.X0

    def set_t_smooth(self, t_smooth):
        """Sets t_smooth as a 2-vector: :math:`t_\mathrm{smooth}[1,1]`.
        """
        import numpy as np
        self.t_smooth=t_smooth*np.array([1,1])

    def get_t_smooth(self):
        """Gets t_smooth (overall amplitude): :math:`t_\mathrm{smooth}[0]`.
        """
        return self.t_smooth[0]

    def set_t_chaos(self, t_chaos):
        """Sets t_chaos as a 2-vector: :math:`t_\mathrm{chaos} [1,1]`.
        """
        import numpy as np
        self.t_chaos=t_chaos*np.array([1,1])

    def get_t_chaos(self):
        """Gets t_chaos (overall amplitude): :math:`t_\mathrm{chaos}[0]`.
        """
        return self.t_chaos[0]

    @property
    def gammaX(self):
        """double: Defined by :math:`\gamma_X \equiv 1/(2 t_\mathrm{smooth}[0]^2)`.
        """
        return 1/(2*self.t_smooth[0]**2)

    @property
    def gammaY(self):
        """double: Defined by :math:`\gamma_Y \equiv 1/(2 t_\mathrm{smooth}[1]^2)`.
        """
        return 1/(2*self.t_smooth[1]**2)

    def _get_rbf(self):
        """Gets the radial basis function (RBF) part of the prior covariance matrix, :math:`\\textbf{K}`.

        Returns
        -------
        K : array, double, dimension=(S,S)
            RBF kernel

        """
        import numpy as np
        from sklearn.metrics.pairwise import rbf_kernel
        t_s=self.t_s

        KX=rbf_kernel(t_s.reshape(-1,1), t_s.reshape(-1,1), gamma=self.gammaX)
        KY=rbf_kernel(t_s.reshape(-1,1), t_s.reshape(-1,1), gamma=self.gammaY)

        return np.block([[KX, np.zeros_like(KX)], [np.zeros_like(KY), KY]])

    def _get_time_variance(self):
        """Gets the time-dependent part of the prior covariance matrix, :math:`\\textbf{V} = \\textbf{u}\\textbf{u}^\intercal`.

        Returns
        -------
        V : array, double, dimension=(S,S)
            time variance matrix

        """
        import numpy as np
        t_s=self.t_s
        alpha_norm=self.alpha_norm
        t_chaos=self.t_chaos

        SX_t=t_s/t_chaos[0]
        SY_t=t_s/t_chaos[1]
        VX=alpha_norm[0]**2 * np.outer(np.ones_like(t_s)+SX_t,np.ones_like(t_s)+SX_t)
        VY=alpha_norm[1]**2 * np.outer(np.ones_like(t_s)+SY_t,np.ones_like(t_s)+SY_t)

        return np.block([[VX, np.zeros_like(VX)], [np.zeros_like(VY), VY]])

    def _get_covariance(self):
        """Gets the full prior covariance matrix as :math:`\\textbf{S}`.

        Returns
        -------
        S : array, double, dimension=(S,S)
            covariance matrix of the prior

        """
        import numpy as np
        K=self.rbf
        V=self.td
        S=np.multiply(V, K)

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
        """Computes the prior (mean, covariance matrix and its inverse).
        """
        import numpy as np
        self.rbf = self._get_rbf()
        self.td = self._get_time_variance()
        self.covariance = self._get_covariance()
        self.inv_covariance = self._get_inverse_covariance()

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
            save_replace_dataset(hf, '/prior/alpha_norm', data=self.t_s, maxshape=(None,), dtype=c_double)
            save_replace_dataset(hf, '/prior/t_smooth', data=self.t_s, maxshape=(None,), dtype=c_double)
            save_replace_dataset(hf, '/prior/t_chaos', data=self.t_s, maxshape=(None,), dtype=c_double)
            save_replace_dataset(hf, '/prior/t_s', data=self.t_s, maxshape=(None,), dtype=c_double)
            save_replace_dataset(hf, '/prior/mean', data=self.mean, maxshape=(None,), dtype=c_double)
            save_replace_dataset(hf, '/prior/rbf', data=self.rbf, maxshape=(None, None), dtype=c_double)
            save_replace_dataset(hf, '/prior/td', data=self.td, maxshape=(None, None), dtype=c_double)
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
            mean = np.array(hf.get('/prior/mean'), dtype=c_double)
            t_s = np.array(hf.get('/prior/t_s'), dtype=c_double)
            alpha_norm = np.array(hf.get('/prior/alpha_norm'), dtype=c_double)
            t_smooth = np.array(hf.get('/prior/t_smooth'), dtype=c_double)
            t_chaos = np.array(hf.get('/prior/t_chaos'), dtype=c_double)
            prior=cls(t_s,mean,alpha_norm,t_smooth,t_chaos)
            prior.rbf = np.array(hf.get('/prior/rbf'), dtype=c_double)
            prior.td = np.array(hf.get('/prior/td'), dtype=c_double)
            prior.covariance = np.array(hf.get('/prior/covariance'), dtype=c_double)
            prior.inv_covariance = np.array(hf.get('/prior/inv_covariance'), dtype=c_double)

        PrintMessage(3, "Reading prior in data file '{}' done.".format(fname))
        return prior
#end class(prior)
