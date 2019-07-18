#!/usr/bin/env python
#-------------------------------------------------------------------------------------
# pySELFI v1.0 -- examples/grf/model/blackbox_GRF.py
# Copyright (C) 2019-2019 Florent Leclercq.
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

"""A blackbox to generate Gaussian random fields and
evaluate their power spectrum, using the Simbelmynë code
"""

__author__  = "Florent Leclercq"
__version__ = "1.0"
__date__    = "2018-2019"
__license__ = "GPLv3"

class blackbox(object):
    """This class represents a SELFI blackbox
    """
    
    # Initialization
    def __init__(self, P, **kwargs):
        """Initializes the blackbox object

        Parameters
        ----------
        P (int) : number of output summary statistics. This is the only mandatory argument
        **kwargs (optional, dictionary) : other arguments

        """
        self.P=P
        for key, value in kwargs.items():
            setattr(self,key,value)
        
        import numpy as np
        if self.seednoise is None and not self.fixnoise:
            self.seednoise = np.random.randint(1e7)
    
    def get_powerspectrum_from_cosmo(self, cosmo, fname_powerspectrum, force=False):
        """Loads or computes the power spectrum from input cosmological parameters

        Parameters
        ----------
        cosmo (dictionary) : cosmological parameters (and some infrastructure parameters)
        fname_powerspectrum (string) : name of input/output power spectrum file
        force (optional, bool, default=False) : force recomputation?

        Returns
        -------
        P (PowerSpectrum) : power spectrum object

        """
        from power import PowerSpectrum
        from os.path import exists
        if exists(fname_powerspectrum) and not force:
            P=PowerSpectrum.read(fname_powerspectrum)
        else:
            G_sim=self.G_sim
            L0=G_sim.L0; L1=G_sim.L1; L2=G_sim.L2; N0=G_sim.N0; N1=G_sim.N2; N2=G_sim.N2
            P=PowerSpectrum(L0,L1,L2,N0,N1,N2,cosmo)
            P.write(fname_powerspectrum)
        return P
    
    def get_powerspectrum_from_theta(self, theta):
        """Returns a power spectrum from its value at support wavenumbers, by performing a Spline interpolation

        Parameters
        ----------
        theta (array, double, dimension=S) : vector of power spectrum values at the support wavenumbers

        Returns
        ------
        P (PowerSpectrum) : power spectrum object

        """
        from power import PowerSpectrum
        from scipy.interpolate import InterpolatedUnivariateSpline
        from os.path import exists
        G_sim=self.G_sim
        P=self.theta2P(theta)
        Spline=InterpolatedUnivariateSpline(self.k_s, P, k=5)
        powerspectrum=Spline(G_sim.k_modes)
        powerspectrum[0]=0. # fix zero-mode by hand
        P=PowerSpectrum.from_FourierGrid(G_sim,powerspectrum=powerspectrum)
        return P
    
    def aux_blackbox(self, P, seedphases=None, seednoise=None):
        """Auxiliary routine for the GRF blackbox: generates a noisy realization from an input power spectrum object, and returns its normalized estimated power spectrum

        Parameters
        ----------
        P (PowerSpectrum) : input power spectrum object
        seedphases (optional, int, default=None) : value of the seed to generate the phases of the Gaussian random field
        seednoise (optional, int, default=None) : value of the seed to generate the noise realization, uses current state if set to None

        Returns
        -------
        Phi (array, double, size=P) : vector of summary statistics

        """
        import numpy as np
        import scipy.stats as ss
        from pysbmy import c_double, c_float_t
        from field import Field
        from correlations import get_autocorrelation
        
        G_sim=self.G_sim
        P_ss=self.P_ss
        corner0=self.corner0
        corner1=self.corner1
        corner2=self.corner2
        a_min=self.a_min
        L0=G_sim.L0; L1=G_sim.L1; L2=G_sim.L2; N0=G_sim.N0; N1=G_sim.N2; N2=G_sim.N2
        
        # Generate Gaussian random field
        F=Field.GRF(L0,L1,L2,corner0,corner1,corner2,N0,N1,N2,P,a_min,seedphases)

        # Add noise
        if seednoise is not None:
            np.random.seed(seednoise)
        F.data += self.noise_std * ss.norm.rvs(size=(N0,N1,N2)).astype(c_float_t)

        # Get power spectrum
        Pk,Vk = get_autocorrelation(F,P_ss.FourierGrid)
        
        # Normalize the blackbox output
        Phi = Pk/P_ss.powerspectrum
        
        # Convert to c_double
        Phi = Phi.astype(c_double)
        return Phi
    
    def make_data(self, cosmo, fname_powerspectrum, seedphases=None, seednoise=None, force=False):
        """Evaluates the GRF blackbox to make mock data, from input cosmological parameters

        Parameters
        ----------
        cosmo (dictionary) : cosmological parameters (and some infrastructure parameters)
        fname_powerspectrum (string) : name of input/output power spectrum file
        seedphases (optional, int, default=None) : value of the seed to generate the phases of the Gaussian random field
        seednoise (optional, int, default=None) : value of the seed to generate the noise realization, uses current state if set to None
        force (optional, bool, default=False) : force recomputation?

        Returns
        -------
        Phi (array, double, size=P) : vector of summary statistics

        """
        from pyselfi.utils import PrintMessage, INDENT, UNINDENT
        
        PrintMessage(3, "Making mock data...")
        INDENT()
        
        # Generate input initial power spectrum from cosmological parameters
        P_data = self.get_powerspectrum_from_cosmo(cosmo,fname_powerspectrum,force)
        
        # Call auxiliary blackbox method
        Phi = self.aux_blackbox(P_data,seedphases,seednoise)
        
        UNINDENT()
        PrintMessage(3, "Making mock data done.")
        return Phi
    
    def evaluate(self, theta, seedphases=None, seednoise=None, i=0, N=0):
        """Evaluates the GRF blackbox from an input vector of power spectrum coefficients at the support wavenumbers

        Parameters
        ----------
        theta (array, double, dimension=S) : vector of power spectrum values at the support wavenumbers
        seedphases (optional, int, default=None) : value of the seed to generate the phases of the Gaussian random field
        seednoise (optional, int, default=None) : value of the seed to generate the noise realization, uses current state if set to None
        i (optional, int, default=0) : current evaluation index of the blackbox
        N (optional, int, default=0) : total number of evaluations of the blackbox

        Returns
        -------
        Phi (array, double, size=P) : vector of summary statistics

        """
        from pyselfi.utils import PrintMessage, PrintValue, INDENT, UNINDENT

        PrintMessage(3, "Evaluating blackbox ({}/{})...".format(i,N))
        INDENT()
        
        PrintValue("seedphases", seedphases)
        PrintValue("seednoise", seednoise)
        
        # Interpolate P to get P_in
        P_in = self.get_powerspectrum_from_theta(theta)
        
        # Call auxiliary blackbox method
        Phi = self.aux_blackbox(P_in,seedphases,seednoise)
        
        UNINDENT()
        PrintMessage(3, "Evaluating blackbox ({}/{}) done.".format(i,N))
        
        return Phi

    def compute_pool(self, theta, d, pool_fname, N):
        """Computes a pool of realizations of the GRF blackbox
        A method compute_pool with this prototype is the only mandatory method

        Parameters
        ----------
        theta (array, double, dimension=S) : vector of power spectrum values at the support wavenumbers
        d (int) : direction in parameter space, from 0 to S
        pool_fname (string) : pool file name
        N (int) : number of realizations required

        """
        def get_current_seed(seedphases,fixphases,seednoise,fixnoise,i):
            if seedphases is None:
                this_seedphases=None
            else:
                if fixphases:
                    this_seedphases=seedphases
                else:
                    this_seedphases=seedphases+i

            if fixnoise:
                this_seednoise=seednoise
            else:
                this_seednoise=seednoise+i
                
            return this_seedphases,this_seednoise
        
        from pyselfi.pool import pool
        p=pool(pool_fname,N)
        
        # Run N evaluations of the blackbox at the desired point
        while not p.finished:
            i = p.N_sims+1
            this_seedphases, this_seednoise = get_current_seed(self.seedphases,self.fixphases,self.seednoise,self.fixnoise,i)
            phi = self.evaluate(theta,this_seedphases,this_seednoise,i,N)
            p.add_sim(phi)
            if i%self.save_frequency==0:
                p.save()
        p.save()
        
        return p
#end class(blackbox)
