#!/usr/bin/env python
#-------------------------------------------------------------------------------------
# pySELFI v1.0 -- examples/simbelmyne/model/blackbox_SBMY.py
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

"""A blackbox to generate mock galaxy observations and
evaluate their power spectrum, using the Simbelynë code
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
    
    def save_cosmo(self, cosmo, fname_cosmo, force_cosmo=False):
        """Saves cosmological parameters in json format

        Parameters
        ----------
        cosmo (dictionary) : cosmological parameters (and some infrastructure parameters) to be saved
        fname_cosmo (string) : name of the output json file
        force_cosmo (optional, bool, default=False) : overwrite if the file already exists?

        """
        import json
        from os.path import exists
        if not exists(fname_cosmo) or force_cosmo:
            with open(fname_cosmo, 'w') as fp:
                json.dump(cosmo, fp)
    
    def get_powerspectrum_from_cosmo(self, cosmo, fname_powerspectrum, force_powerspectrum=False):
        """Loads or computes the power spectrum from input cosmological parameters

        Parameters
        ----------
        cosmo (dictionary) : cosmological parameters (and some infrastructure parameters)
        fname_powerspectrum (string) : name of input/output power spectrum file
        force_powerspectrum (optional, bool, default=False) : force recomputation?

        """
        from os.path import exists
        from pysbmy.power import PowerSpectrum
        if exists(fname_powerspectrum) and not force_powerspectrum:
            P=PowerSpectrum.read(fname_powerspectrum)
        else:
            G_sim=self.G_sim
            L0=G_sim.L0; L1=G_sim.L1; L2=G_sim.L2; N0=G_sim.N0; N1=G_sim.N2; N2=G_sim.N2
            P=PowerSpectrum(L0,L1,L2,N0,N1,N2,cosmo)
            P.write(fname_powerspectrum)
    
    def get_powerspectrum_from_theta(self, theta, fname_powerspectrum, force_powerspectrum=False):
        """Returns a power spectrum from its value at support wavenumbers,
        by performing a Spline interpolation

        Parameters
        ----------
        theta (array, double, dimension=S) : vector of power spectrum values at the support wavenumbers
        fname_powerspectrum (string) : name of input/output power spectrum file
        force_powerspectrum (optional, bool, default=False) : force recomputation?

        """
        from os.path import exists
        from scipy.interpolate import InterpolatedUnivariateSpline
        from pysbmy.power import PowerSpectrum
        if exists(fname_powerspectrum) and not force_powerspectrum:
            P=PowerSpectrum.read(fname_powerspectrum)
        else:
            G_sim=self.G_sim
            P=self.theta2P(theta)
            Spline=InterpolatedUnivariateSpline(self.k_s, P, k=5)
            powerspectrum=Spline(G_sim.k_modes)
            powerspectrum[0]=0. # fix zero-mode by hand
            P=PowerSpectrum.from_FourierGrid(G_sim,powerspectrum=powerspectrum)
            P.write(fname_powerspectrum)
        
    def setup_parfiles(self, d, fname_simparfile, fname_mockparfile, fname_powerspectrum, fname_RngStateLPT, fname_whitenoise, fname_outputinitialdensity, fname_outputrealspacedensity, fname_outputdensity, fname_RngStateMock, fname_outputmock, fname_outputss, force_parfiles=False):
        """Setup Simbelynë parameter file given the necessary inputs (see Simbelynë documentation for details)

        Parameters
        ----------
        d (int) : index giving the direction in parameter space: -1 for mock data, 0 for the expansion point, or from 1 to S
        fname_simparfile (string) : name of output simulation parameter file
        fname_mockparfile (string) : name of output mock parameter file
        fname_powerspectrum (string) : name of input power spectrum file
        fname_RngStateLPT (string) : name of output random number generator state file for LPT module
        fname_whitenoise (string) : name of output white noise file
        fname_outputinitialdensity (string) : name of output initial density field file
        fname_outputrealspacedensity (string) : name of output real-space density field file
        fname_outputdensity (string) : name of output redshift-space density field file
        fname_RngStateMock (string) : name of output random number generator state file for Mock module
        fname_outputmock (string) : name of output mock field file
        fname_outputss (string) : name of output summary statistics file
        force_parfiles (optional, bool, default=False) : overwrite if files already exists?

        """
        from os.path import exists
        from pysbmy import param_file
        
        Particles=self.Np0
        Mesh=self.G_sim.N0
        BoxSize=self.G_sim.L0
        corner0=self.corner0
        corner1=self.corner1
        corner2=self.corner2
        InputRngStateLPT=self.fsimdir+"/RngStateLPT.rng"
        OutputRngStateLPT=self.fsimdir+"/RngStateLPT.rng"
        WriteInitialConditions=0
        RedshiftLPT=19.
        
        ModulePMCOLA=1
        EvolutionMode=2
        ParticleMesh=self.Npm0
        NumberOfTimeSteps=20
        RedshiftFCs=0
        WriteFinalSnapshot=0
        WriteFinalDensity=0

        InputRngStateMocks=self.fsimdir+"/RngStateMocks.rng"
        OutputRngStateMocks=self.fsimdir+"/RngStateMocks.rng"
        NoiseModel=0
        NumberOfNoiseRealizations=1
        InputSurveyGeometry=self.fsimdir+"/"+self.fname_inputsurveygeometry
        InputSummaryStatskGrid=self.fdir+"G_ss.h5"
        WriteMocks=0
        
        if not exists(fname_simparfile) or force_parfiles:
            if d==-1: # mock data
                S=param_file(ModuleLPT=1,
                            Particles=Particles,
                            Mesh=Mesh,
                            BoxSize=BoxSize,
                            corner0=corner0,
                            corner1=corner1,
                            corner2=corner2,
                            InputRngStateLPT=InputRngStateLPT,
                            OutputRngStateLPT=OutputRngStateLPT,
                            ICsMode=0,
                            WriteICsRngState=1,
                            OutputICsRngState=fname_RngStateLPT,
                            WriteWhiteNoise=1,
                            OutputWhiteNoise=fname_whitenoise,
                            WriteInitialConditions=1,
                            OutputInitialConditions=fname_outputinitialdensity,
                            InputPowerSpectrum=fname_powerspectrum,
                            RedshiftLPT=RedshiftLPT,
                            WriteLPTSnapshot=0,
                            WriteLPTDensity=0,
                            
                            ModulePMCOLA=ModulePMCOLA,
                            EvolutionMode=EvolutionMode,
                            ParticleMesh=ParticleMesh,
                            NumberOfTimeSteps=NumberOfTimeSteps,
                            RedshiftFCs=RedshiftFCs,
                            WriteFinalSnapshot=WriteFinalSnapshot,
                            WriteFinalDensity=1,
                            OutputFinalDensity=fname_outputrealspacedensity,
                            
                            ModuleRSD=1,
                            DoNonLinearMapping=1,
                            WriteRSSnapshot=0,
                            WriteRSDensity=1,
                            OutputRSDensity=fname_outputdensity
                            )
            elif d==0: # expansion point
                S=param_file(ModuleLPT=1,
                            Particles=Particles,
                            Mesh=Mesh,
                            BoxSize=BoxSize,
                            corner0=corner0,
                            corner1=corner1,
                            corner2=corner2,
                            InputRngStateLPT=InputRngStateLPT,
                            OutputRngStateLPT=OutputRngStateLPT,
                            ICsMode=0,
                            WriteICsRngState=1,
                            OutputICsRngState=fname_RngStateLPT,
                            WriteWhiteNoise=1,
                            OutputWhiteNoise=fname_whitenoise,
                            WriteInitialConditions=WriteInitialConditions,
                            OutputInitialConditions=fname_outputinitialdensity,
                            InputPowerSpectrum=fname_powerspectrum,
                            RedshiftLPT=RedshiftLPT,
                            WriteLPTSnapshot=0,
                            WriteLPTDensity=0,
                            
                            ModulePMCOLA=ModulePMCOLA,
                            EvolutionMode=EvolutionMode,
                            ParticleMesh=ParticleMesh,
                            NumberOfTimeSteps=NumberOfTimeSteps,
                            RedshiftFCs=RedshiftFCs,
                            WriteFinalSnapshot=WriteFinalSnapshot,
                            WriteFinalDensity=WriteFinalDensity,
                            OutputFinalDensity=fname_outputrealspacedensity,
                            
                            ModuleRSD=1,
                            DoNonLinearMapping=1,
                            WriteRSSnapshot=0,
                            WriteRSDensity=1,
                            OutputRSDensity=fname_outputdensity
                            )
            else: # points needed for gradients
                S=param_file(ModuleLPT=1,
                            Particles=Particles,
                            Mesh=Mesh,
                            BoxSize=BoxSize,
                            corner0=corner0,
                            corner1=corner1,
                            corner2=corner2,
                            InputRngStateLPT=InputRngStateLPT,
                            OutputRngStateLPT=OutputRngStateLPT,
                            ICsMode=1,
                            WriteICsRngState=0,
                            InputWhiteNoise=fname_whitenoise,
                            WriteInitialConditions=0,
                            InputPowerSpectrum=fname_powerspectrum,
                            RedshiftLPT=RedshiftLPT,
                            WriteLPTSnapshot=0,
                            WriteLPTDensity=0,
                            
                            ModulePMCOLA=ModulePMCOLA,
                            EvolutionMode=EvolutionMode,
                            ParticleMesh=ParticleMesh,
                            NumberOfTimeSteps=NumberOfTimeSteps,
                            RedshiftFCs=RedshiftFCs,
                            WriteFinalSnapshot=WriteFinalSnapshot,
                            WriteFinalDensity=WriteFinalDensity,
                            
                            ModuleRSD=1,
                            DoNonLinearMapping=1,
                            WriteRSSnapshot=0,
                            WriteRSDensity=1,
                            OutputRSDensity=fname_outputdensity
                            )                
            S.write(fname_simparfile)
        
        if not exists(fname_mockparfile) or force_parfiles:
            if d<=0: # mock data or expansion point
                S=param_file(ModuleLPT=0,
                            
                            ModulePMCOLA=0,
                            
                            ModuleRSD=0,
                            
                            ModuleMocks=1,
                            InputRngStateMocks=InputRngStateMocks,
                            OutputRngStateMocks=OutputRngStateMocks,
                            WriteMocksRngState=1,
                            OutputMocksRngState=fname_RngStateMock,
                            NoiseModel=NoiseModel,
                            NumberOfNoiseRealizations=NumberOfNoiseRealizations,
                            InputDensityMocks=fname_outputdensity,
                            InputSurveyGeometry=InputSurveyGeometry,
                            InputSummaryStatskGrid=InputSummaryStatskGrid,
                            WriteMocks=WriteMocks,
                            OutputMockBase=fname_outputmock,
                            WriteSummaryStats=1,
                            OutputSummaryStats=fname_outputss
                            )
            else: # points needed for gradients
                S=param_file(ModuleLPT=0,
                            
                            ModulePMCOLA=0,
                            
                            ModuleRSD=0,
                            
                            ModuleMocks=1,
                            InputRngStateMocks=fname_RngStateMock,
                            WriteMocksRngState=0,
                            NoiseModel=NoiseModel,
                            NumberOfNoiseRealizations=NumberOfNoiseRealizations,
                            InputDensityMocks=fname_outputdensity,
                            InputSurveyGeometry=InputSurveyGeometry,
                            InputSummaryStatskGrid=InputSummaryStatskGrid,
                            WriteMocks=WriteMocks,
                            OutputMockBase=fname_outputmock,
                            WriteSummaryStats=1,
                            OutputSummaryStats=fname_outputss
                            )
            S.write(fname_mockparfile)
    
    def run_sim(self, fname_simparfile, fname_simlogs, fname_outputdensity, force_sim=False):
        """Run simulation with Simbelynë

        Parameters
        ----------
        fname_simparfile (string) : name of the input parameter file
        fname_simlogs (string) : name of the output Simbelynë logs
        fname_outputdensity (string) : name of the output density field to be written
        force_sim (optional, bool, default=False) : force recomputation if output density already exists?

        """
        from os.path import exists
        if not exists(fname_outputdensity) or force_sim:
            from pysbmy import pySbmy
            pySbmy(fname_simparfile, fname_simlogs)
           
    def run_mock(self, fname_mockparfile, fname_mocklogs, fname_outputss, force_mock=False):
        """Run simulation with Simbelynë for mock data

        Parameters
        ----------
        fname_mockparfile (string) : name of the input parameter file
        fname_mocklogs (string) : name of the output Simbelynë logs
        fname_outputss (string) : name of the output summary statistics file to be written
        force_mock (optional, bool, default=False) : force recomputation if output summary statistics file already exists?

        """
        from os.path import exists
        if not exists(fname_outputss) or force_mock:
            from pysbmy import pySbmy
            pySbmy(fname_mockparfile, fname_mocklogs)
    
    def compute_Phi(self, fname_outputss):
        """Compute SELFI summary statistics (Phi) from Simbelynë output files

        Parameters
        ----------
        fname_outputss (string) : name of Simbelynë output summary statistics file

        Returns
        -------
        Phi (array, double, size=P) : vector of summary statistics

        """
        import h5py as h5
        from pysbmy import c_float
        
        # Read the Simbelynë output
        with h5.File(fname_outputss) as f:
            # Normalize the output
            Phi = f['scalars']['Pk'][0][0]/self.P_ss.powerspectrum
            
            # Convert to c_float
            Phi = Phi.astype(c_float)
        
        return Phi
    
    def compute_Phi_DM(self, fname_outputdensity):
        """Compute SELFI summary statistics (Phi) from Simbelynë output files
        Alternative routine to compute_Phi using the dark matter field instead of galaxies, for testing purposes

        Parameters
        ----------
        fname_outputdensity (string) : name of Simbelynë output density field file

        Returns
        -------
        Phi (array, double, size=P) : vector of summary statistics

        """
        from pysbmy.field import read_field
        from pysbmy.correlations import get_autocorrelation
        
        A=read_field(fname_outputdensity)
        Pk,Vk=get_autocorrelation(A, self.G_ss)
        return Pk/self.P_ss.powerspectrum
        
    def clean_sbmy_outputs(self, fname_outputdensity, fname_outputmock, fname_outputss):
        """Remove Simbelynë output files to save disk space

        Parameters
        ----------
        fname_outputdensity (string) : name of output redshift-space density field file
        fname_outputmock (string) : name of output mock field file
        fname_outputss (string) : name of output summary statistics file

        """
        from os.path import exists
        from os import remove
        
        # Remove the Simbelynë outputs
        if exists(fname_outputdensity): remove(fname_outputdensity)
        if exists(fname_outputmock): remove(fname_outputmock)
        if exists(fname_outputss): remove(fname_outputss)
    
    def aux_blackbox(self, d, fname_powerspectrum, fname_simparfile, fname_mockparfile, fname_RngStateLPT, fname_whitenoise, fname_outputinitialdensity, fname_outputrealspacedensity, fname_outputdensity, fname_RngStateMock, fname_outputmock, fname_outputss, fname_simlogs, fname_mocklogs, force_parfiles=False, force_sim=False, force_mock=False):
        """Auxiliary routine for the Simbelynë blackbox: generates a noisy realization from an input power spectrum object, and returns its normalized estimated power spectrum

        Parameters
        ----------
        d (int) : index giving the direction in parameter space: -1 for mock data, 0 for the expansion point, or from 1 to S
        fname_powerspectrum (string) : name of input power spectrum file
        fname_simparfile (string) : name of output simulation parameter file
        fname_mockparfile (string) : name of output mock parameter file
        fname_RngStateLPT (string) : name of output random number generator state file for LPT module
        fname_whitenoise (string) : name of output white noise file
        fname_outputinitialdensity (string) : name of output initial density field file
        fname_outputrealspacedensity (string) : name of output real-space density field file
        fname_outputdensity (string) : name of output redshift-space density field file
        fname_RngStateMock (string) : name of output random number generator state file for Mock module
        fname_outputmock (string) : name of output mock field file
        fname_outputss (string) : name of output summary statistics file
        fname_simlogs (string) : name of the output Simbelynë logs
        fname_mocklogs (string) : name of the output Simbelynë logs
        force_parfiles (optional, bool, default=False) : overwrite if parameter files already exists?
        force_sim (optional, bool, default=False) : force recomputation if output density already exists?
        force_mock (optional, bool, default=False) : force recomputation if output mock summary statistics file already exists?

        """
        # Prepare parameter files
        self.setup_parfiles(d, fname_simparfile, fname_mockparfile, fname_powerspectrum, fname_RngStateLPT, fname_whitenoise, fname_outputinitialdensity, fname_outputrealspacedensity, fname_outputdensity, fname_RngStateMock, fname_outputmock, fname_outputss, force_parfiles)
        
        # Run simulations
        self.run_sim(fname_simparfile, fname_simlogs, fname_outputdensity, force_sim)
        self.run_mock(fname_mockparfile, fname_mocklogs, fname_outputss, force_mock)
        Phi = self.compute_Phi(fname_outputss)
        #Phi = self.compute_Phi_DM(fname_outputdensity) # TEST: Alternatively, the dark matter field can be used instead of galaxies
        
        return Phi
    
    def make_data(self, cosmo, i=0, force_powerspectrum=False, force_parfiles=False, force_sim=False, force_mock=False, force_cosmo=False):
        """Evaluates the Simbelynë blackbox to make mock data, from input cosmological parameters

        Parameters
        ----------
        cosmo (dictionary) : cosmological parameters (and some infrastructure parameters)
        i (optional, int, default=0) : current evaluation index of the blackbox
        force_powerspectrum (optional, bool, default=False) : force recomputation of the power spectrum?
        force_parfiles (optional, bool, default=False) : overwrite if parameter files already exists?
        force_sim (optional, bool, default=False) : force recomputation if output density already exists?
        force_mock (optional, bool, default=False) : force recomputation if output mock summary statistics file already exists?

        Returns
        -------
        Phi (array, double, size=P) : vector of summary statistics

        """
        from pyselfi.utils import PrintMessage, INDENT, UNINDENT
        
        PrintMessage(3, "Making mock data...")
        INDENT()

        # Define filenames
        fname_cosmo = self.fsimdir+"/data/input_cosmo_"+str(i)+".json"
        fname_powerspectrum = self.fsimdir+"/data/input_power_"+str(i)+".h5"
        fname_simparfile = self.fsimdir+"/data/sim_"+str(i)+".sbmy"
        fname_mockparfile = self.fsimdir+"/data/mock_"+str(i)+".sbmy"
        fname_RngStateLPT = self.fsimdir+"/data/RngStateLPT_"+str(i)+".rng"
        fname_whitenoise = self.fsimdir+"/data/initial_density_white_"+str(i)+".h5"
        fname_outputinitialdensity = self.fsimdir+"/data/initial_density_"+str(i)+".h5"
        fname_outputrealspacedensity = self.fsimdir+"/data/output_realspace_density_"+str(i)+".h5"
        fname_outputdensity = self.fsimdir+"/data/output_density_"+str(i)+".h5"
        fname_RngStateMock = self.fsimdir+"/data/RngStateMocks_"+str(i)+".rng"
        fname_outputmock = self.fsimdir+"/data/output_mock_"+str(i)+"_"
        fname_outputss = self.fsimdir+"/data/output_ss_"+str(i)+".h5"
        fname_simlogs = self.fsimdir+"/data/logs_sim_"+str(i)+".txt"
        fname_mocklogs = self.fsimdir+"/data/logs_mock_"+str(i)+".txt"
        
        # Save cosmological parameters
        self.save_cosmo(cosmo, fname_cosmo, force_cosmo)
        
        # Generate input initial power spectrum from cosmological parameters
        self.get_powerspectrum_from_cosmo(cosmo, fname_powerspectrum, force_powerspectrum)
        
        # Call auxiliary blackbox method
        Phi = self.aux_blackbox(-1, fname_powerspectrum, fname_simparfile, fname_mockparfile, fname_RngStateLPT, fname_whitenoise, fname_outputinitialdensity, fname_outputrealspacedensity, fname_outputdensity, fname_RngStateMock, fname_outputmock, fname_outputss, fname_simlogs, fname_mocklogs, force_parfiles, force_sim, force_mock)
        
        UNINDENT()
        PrintMessage(3, "Making mock data done.")
        return Phi
    
    def evaluate(self, theta, d, i=0, N=0, force_powerspectrum=False, force_parfiles=False, force_sim=False, force_mock=False, remove_sbmy=True):
        """Evaluates the Simbelynë blackbox from an input vector of power spectrum coefficients at the support wavenumbers

        Parameters
        ----------
        theta (array, double, dimension=S) : vector of power spectrum values at the support wavenumbers
        d (int) : index giving the direction in parameter space: -1 for mock data, 0 for the expansion point, or from 1 to S
        i (optional, int, default=0) : current evaluation index of the blackbox
        N (optional, int, default=0) : total number of evaluations of the blackbox
        force_powerspectrum (optional, bool, default=False) : force recomputation of the power spectrum?
        force_parfiles (optional, bool, default=False) : overwrite if parameter files already exists?
        force_sim (optional, bool, default=False) : force recomputation if output density already exists?
        force_mock (optional, bool, default=False) : force recomputation if output mock summary statistics file already exists?
        remove_sbmy (optional, bool, default=True) : remove Simbelynë output files from disk?

        Returns
        -------
        Phi (array, double, size=P) : vector of summary statistics

        """
        from pyselfi.utils import PrintMessage, PrintValue, INDENT, UNINDENT

        PrintMessage(3, "Evaluating blackbox ({}/{})...".format(i,N))
        INDENT()
        
        # Define filenames
        fname_powerspectrum = self.fsimdir+"/d"+str(d)+"/input_power_d"+str(d)+".h5"
        fname_simparfile = self.fsimdir+"/d"+str(d)+"/sim_d"+str(d)+"_p"+str(i-1)+".sbmy"
        fname_mockparfile = self.fsimdir+"/d"+str(d)+"/mock_d"+str(d)+"_p"+str(i-1)+".sbmy"
        fname_RngStateLPT = self.fsimdir+"/RngStateLPT_p"+str(i-1)+".rng"
        fname_whitenoise = self.fsimdir+"/initial_density_white_p"+str(i-1)+".h5"
        fname_outputinitialdensity = self.fsimdir+"/initial_density_p"+str(i-1)+".h5"
        fname_outputrealspacedensity = self.fsimdir+"/d"+str(d)+"/output_realspace_density_d"+str(d)+"_p"+str(i-1)+".h5"
        fname_outputdensity = self.fsimdir+"/d"+str(d)+"/output_density_d"+str(d)+"_p"+str(i-1)+".h5"
        fname_RngStateMock = self.fsimdir+"/RngStateMocks_p"+str(i-1)+".rng"
        fname_outputmock = self.fsimdir+"/d"+str(d)+"/output_mock_d"+str(d)+"_p"+str(i-1)+"_"
        fname_outputss = self.fsimdir+"/d"+str(d)+"/output_ss_d"+str(d)+"_p"+str(i-1)+".h5"
        fname_simlogs = self.fsimdir+"/d"+str(d)+"/logs_sim_d"+str(d)+"_p"+str(i-1)+".txt"
        fname_mocklogs = self.fsimdir+"/d"+str(d)+"/logs_mock_d"+str(d)+"_p"+str(i-1)+".txt"
        
        # Interpolate P to get P_in
        self.get_powerspectrum_from_theta(theta, fname_powerspectrum, force_powerspectrum)
        
        # Call auxiliary blackbox method
        Phi = self.aux_blackbox(d, fname_powerspectrum, fname_simparfile, fname_mockparfile, fname_RngStateLPT, fname_whitenoise, fname_outputinitialdensity, fname_outputrealspacedensity, fname_outputdensity, fname_RngStateMock, fname_outputmock, fname_outputss, fname_simlogs, fname_mocklogs, force_parfiles, force_sim, force_mock)
        
        # Clean Simbelynë outputs
        if remove_sbmy:
            self.clean_sbmy_outputs(fname_outputdensity, fname_outputmock, fname_outputss)
        
        UNINDENT()
        PrintMessage(3, "Evaluating blackbox ({}/{}) done.".format(i,N))
        
        return Phi

    def compute_pool(self, theta, d, pool_fname, N, force_powerspectrum=False, force_parfiles=False, force_sim=False, force_mock=False):
        """Computes a pool of realizations of the Simbelynë blackbox
        A method compute_pool with this prototype is the only mandatory method

        Parameters
        ----------
        theta (array, double, dimension=S) : vector of power spectrum values at the support wavenumbers
        d (int) : direction in parameter space, from 0 to S
        pool_fname (string) : pool file name
        N (int) : number of realizations required
        force_powerspectrum (optional, bool, default=False) : force recomputation of the power spectrum?
        force_parfiles (optional, bool, default=False) : overwrite if parameter files already exists?
        force_sim (optional, bool, default=False) : force recomputation if output density already exists?
        force_mock (optional, bool, default=False) : force recomputation if output mock summary statistics file already exists?

        """
        from pyselfi.pool import pool
        p=pool(pool_fname,N)
        
        # Run N evaluations of the blackbox at the desired point
        while not p.finished:
            i = p.N_sims+1
            phi = self.evaluate(theta,d,i,N,force_powerspectrum,force_parfiles,force_sim,force_mock)
            p.add_sim(phi)
            if i%self.save_frequency==0:
                p.save()
        p.save()
        
        return p
    
    def make_survey_geometry(self, N_CAT, cosmo, force=False):
        """Produces a mock survey geometry file

        Parameters
        ----------
        N_CAT (int) : number of subcatalogs
        cosmo (dictionary) : cosmological parameters (and some infrastructure parameters)
        force (optional, bool, default=False) : force recomputation if file already exists?

        """
        from os.path import exists
        fname = self.fsimdir+"/"+self.fname_inputsurveygeometry
        if not exists(fname) or force:    
            import numpy as np
            from pysbmy.survey_geometry import SurveyGeometry
            from pysbmy.cosmology import redshift_from_comoving_distance, luminosity_distance

            # basic setup
            G_sim=self.G_sim
            L0=G_sim.L0; L1=G_sim.L1; L2=G_sim.L2; N0=G_sim.N0; N1=G_sim.N2; N2=G_sim.N2
            corner0=self.corner0; corner1=self.corner1; corner2=self.corner2

            # observational setup
            bright_cut=np.array(((0.),))
            faint_cut=np.array(((0.),))
            rmin=np.array(((0.),))
            rmax=np.array(((L0/2*np.sqrt(3)),))
            zmin=redshift_from_comoving_distance(rmin, **cosmo)
            zmax=redshift_from_comoving_distance(rmax, **cosmo)
            N_BIAS=3
            nbar=2.0e-3
            Nmean=nbar*(L0*L1*L2)/(N0*N1*N2)
            galaxy_bias_mean=np.array(((1.2,0.,0.),))
            galaxy_bias_std=np.array(((0.,0.,0.),))
            galaxy_nmean_mean=np.array(((Nmean),))
            galaxy_nmean_std=np.array(((0.),))
            galaxy_sel_window=np.ones((N_CAT,N0,N1,N2))

            # setup coordinates
            x=np.arange(corner0,corner0+L0,L0/N0)
            x=np.tile(x,(N2,N1,1))+1/2.*L0/N0

            y=np.arange(corner1,corner1+L1,L1/N1)
            y=np.tile(y,(N2,N0,1))
            y=np.swapaxes(y,1,2)+1/2.*L1/N1

            z=np.arange(corner2,corner2+L2,L2/N2)
            z=np.tile(z,(N1,N0,1))
            z=np.swapaxes(z,0,2)+1/2.*L2/N2

            r=np.sqrt(x**2+y**2+z**2)
            b=np.arcsin(z/r)*180/np.pi
            l=np.arctan2(y,x)*180/np.pi
            l[np.where(l<0)]+=360.
                        
            # setup basic mask
            mask=np.ones((N0,N1,N2))
            mask[np.where((b>-self.b_cut)*(b<self.b_cut))]=0.

            # Schechter selection function
            def integrand_luminosity(x, alpha):
                from math import exp
                return x**alpha * exp(-x)
            
            def integral_luminosity(alpha, x_min, x_max):
                from scipy.integrate import quad
                return quad(integrand_luminosity, x_min, x_max, args=(alpha))[0]
            
            def computeSchechterPhi1(faint_absolute_magnitude_cut, bright_absolute_magnitude_cut, faint_apparent_magnitude_cut, bright_apparent_magnitude_cut, Mstar, alpha):
                from math import log10, pow

                xl1 = pow(10.0, 0.4*(Mstar - faint_absolute_magnitude_cut))
                xu1 = pow(10.0, 0.4*(Mstar - bright_absolute_magnitude_cut))

                Phi1 = integral_luminosity(alpha, xl1, xu1)
                
                return Phi1
                
            def computeSchechterCompleteness(d_lum, cosmo, faint_absolute_magnitude_cut, bright_absolute_magnitude_cut, faint_apparent_magnitude_cut, bright_apparent_magnitude_cut, Mstar, alpha, Phi1, corr=0.):
                from math import log10, pow
                
                #corr = zcorrection(redshift)
                corr = 0.
                
                absolute_mu0 = faint_apparent_magnitude_cut  - 5 * log10(d_lum) - 25 - corr;
                absolute_ml0 = bright_apparent_magnitude_cut - 5 * log10(d_lum) - 25 - corr;
                
                abmu = min(absolute_mu0, faint_absolute_magnitude_cut)
                abml = max(absolute_ml0, bright_absolute_magnitude_cut)
                
                abmu = max(abmu,abml)

                xl0 = pow(10.0, 0.4*(Mstar- abmu))
                xu0 = pow(10.0, 0.4*(Mstar- abml))
                                
                Phi0 = integral_luminosity(alpha, xl0, xu0)
                
                return max(0.0, Phi0/Phi1)
        
            selection=np.vectorize(computeSchechterCompleteness, excluded=['**cosmo', 'faint_absolute_magnitude_cut', 'bright_absolute_magnitude_cut', 'faint_apparent_magnitude_cut', 'bright_apparent_magnitude_cut','Mstar','alpha','Phi1','corr'])

            # setup galaxy selection window
            from scipy.interpolate import UnivariateSpline
            N = 1000
            r_arr = np.linspace(rmin+1e-7,rmax,N)
            z_arr = redshift_from_comoving_distance(r_arr, **cosmo)
            d_lum_arr = luminosity_distance(z_arr, **cosmo)
            
            bright_apparent_magnitude_cut = self.bright_apparent_magnitude_cut
            faint_apparent_magnitude_cut = self.faint_apparent_magnitude_cut
            bright_absolute_magnitude_cut = self.bright_absolute_magnitude_cut
            faint_absolute_magnitude_cut = self.faint_absolute_magnitude_cut
            Mstar = self.Mstar
            alpha = self.alpha
            Phi1 = computeSchechterPhi1(faint_absolute_magnitude_cut, bright_absolute_magnitude_cut, faint_apparent_magnitude_cut, bright_apparent_magnitude_cut, Mstar, alpha)

            selection_arr = selection(d_lum_arr, cosmo,
                                            faint_absolute_magnitude_cut,
                                            bright_absolute_magnitude_cut,
                                            faint_apparent_magnitude_cut,
                                            bright_apparent_magnitude_cut,
                                            Mstar, alpha, Phi1)
            spl = UnivariateSpline(r_arr, selection_arr)
            
            galaxy_sel_window[0]*=spl(r)*mask
            galaxy_sel_window[0]/=galaxy_sel_window[0].max()

            S=SurveyGeometry(L0,L1,L2,corner0,corner1,corner2,N0,N1,N2,N_CAT,cosmo,bright_cut,faint_cut,rmin,rmax,zmin,zmax,N_BIAS,galaxy_bias_mean,galaxy_bias_std,galaxy_nmean_mean,galaxy_nmean_std,galaxy_sel_window)
            S.write(fname)
    
    def make_ss_k_grid(self, force=False):
        """Produces the summary statistics Fourier grid

        Parameters
        ----------
        force (optional, bool, default=False) : force recomputation if file already exists?

        """
        from os.path import exists
        fname = self.fsimdir+"/"+self.fname_inputsskgrid
        if not exists(fname) or force:
            self.G_ss.write(fname)
#end class(blackbox)
