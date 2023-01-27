#!/usr/bin/env python
#-------------------------------------------------------------------------------------
# pySELFI v2.0 -- src/pyselfi_examples/simbelmyne/model/setup_SBMY.py
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

"""Setup script for the pySELFI Galaxy survey example
"""

__author__  = "Florent Leclercq"
__version__ = "2.0"
__date__    = "2018-2023"
__license__ = "GPLv3"

import numpy as np
import h5py as h5
from os.path import dirname, realpath, exists
from pysbmy.power import PowerSpectrum, FourierGrid, get_Pk
from pyselfi_examples.simbelmyne.model.blackbox_SBMY import blackbox
from pyselfi.power_spectrum.prior import power_spectrum_prior
from pyselfi.power_spectrum.selfi import power_spectrum_selfi

# Iteration
i_iter=0

# Basic configuration
fdir=realpath(dirname(realpath(__file__))+"/../sims/")+"/"
force=False
L0=L1=L2=1000.
corner0=corner1=corner2=-500.
N0=N1=N2=256
Np0=Np1=Np2=512
Npm0=Npm1=Npm2=1024
if not exists(fdir+"G_sim.h5") or force:
    G_sim=FourierGrid(L0,L1,L2,N0,N1,N2)
    G_sim.write(fdir+"G_sim.h5")
else:
    G_sim=FourierGrid.read(fdir+"G_sim.h5")
k_th=G_sim.k_modes

# Support wavenumbers: take the first N_Gks bins from the simulation, then log-spaced up to k_max
S=100
N_Gks=8
k_max=1.4
k_s=np.concatenate([G_sim.k_modes[1:N_Gks+1],np.logspace(np.log10(G_sim.k_modes[N_Gks+1]),np.log10(k_max),S-N_Gks)])
k_s[0]-=1e-6
k_min=k_s[0]

# Expansion point: a BBKS power spectrum with Planck cosmology for the 0th iteration
planck_mean=np.array((0.6774, 0.04860, 0.3089, 0.9667, 0.8159))
planck_cov=np.diag(np.array(((0.0046)**2, (0.00030)**2, (0.0062)**2, (0.0040)**2, (0.0086)**2)))
omega=planck_mean
cosmo_exp={'h':omega[0], 'Omega_r':0., 'Omega_q':1.-omega[2], 'Omega_b':omega[1], 'Omega_m':omega[2], 'm_ncdm':0., 'Omega_k':0., 'tau_reio':0.066, 'n_s':omega[3], 'sigma8':omega[4], 'w0_fld':-1., 'wa_fld':0., 'k_max':10.0, 'WhichSpectrum':"BBKS"}
if not exists(fdir+"P_exp.h5") or force:
    P_exp=PowerSpectrum(L0,L1,L2,N0,N1,N2,cosmo_exp)
    P_exp.write(fdir+"P_exp.h5")
else:
    P_exp=PowerSpectrum.read(fdir+"P_exp.h5")

if i_iter==0:
    theta_0=np.ones(S)
else:
    with h5.File(fdir+"abc_iter"+str(i_iter-1)+".h5") as hf:
        theta_0=np.array(hf["posterior"]["mean"])

# Scaling of the inference variable
if not exists(fdir+"P_0.npy") or force:
    P_0=get_Pk(k_s,cosmo_exp)
    np.save(fdir+"P_0",P_0)
else:
    P_0=np.load(fdir+"P_0.npy")
def theta2P(theta):
    return theta*P_0

def P2theta(P):
    return P/P_0

# Prior hyperparameters
theta_norm=0.20
k_corr=0.015
alpha_cv=0.0008847681471875314
log_kmodes=False
prior=power_spectrum_prior(k_s,theta_0,theta_norm,k_corr,alpha_cv,log_kmodes)

# Finite differencing setup
Ns=100
Ne=150
Delta_theta=1e-2

# Simulator setup
fname_inputsurveygeometry="input_survey_geometry.h5"
P=50
k_ss_min=2.0e-2
k_ss_max=0.5
k_ss_max_offset=0.008
trim_threshold=100
k_ss=np.logspace(np.log10(k_ss_min),np.log10(k_ss_max),P,dtype=np.float32)
if not exists(fdir+"G_ss.h5") or force:
    G_ss=FourierGrid(L0,L1,L2,N0,N1,N2,k_modes=k_ss,kmax=k_ss_max+k_ss_max_offset,trim_bins=True,trim_threshold=trim_threshold)
    G_ss.write(fdir+"G_ss.h5")
else:
    G_ss=FourierGrid.read(fdir+"G_ss.h5")
P=G_ss.NUM_MODES
k_ss=G_ss.k_modes
if not exists(fdir+"P_ss.h5") or force:
    P_0_ss=get_Pk(k_ss,cosmo_exp)
    P_ss=PowerSpectrum.from_FourierGrid(G_ss,powerspectrum=P_0_ss,cosmo=cosmo_exp)
    P_ss.write(fdir+"P_ss.h5")
else:
    P_ss=PowerSpectrum.read(fdir+"P_ss.h5")

b_cut=10.
bright_apparent_magnitude_cut = 0
faint_apparent_magnitude_cut = 18.5
bright_absolute_magnitude_cut = -25.00
faint_absolute_magnitude_cut = -21.00
Mstar = -20.44
alpha = -1.05

save_frequency=1
blackbox=blackbox(P,theta2P=theta2P,k_s=k_s,G_sim=G_sim,G_ss=G_ss,P_ss=P_ss,corner0=corner0,corner1=corner1,corner2=corner2,Np0=Np0,Npm0=Npm0,fdir=fdir,fsimdir=fdir,fname_inputsurveygeometry=fname_inputsurveygeometry,b_cut=b_cut,bright_apparent_magnitude_cut=bright_apparent_magnitude_cut,faint_apparent_magnitude_cut=faint_apparent_magnitude_cut,bright_absolute_magnitude_cut=bright_absolute_magnitude_cut,faint_absolute_magnitude_cut=faint_absolute_magnitude_cut,Mstar=Mstar,alpha=alpha,save_frequency=save_frequency)

# Survey geometry
blackbox.make_survey_geometry(1,cosmo_exp,force=force)

# Observed data
np.random.seed(181249)
omega=np.random.multivariate_normal(planck_mean, planck_cov)
cosmo_obs={'h':omega[0], 'Omega_r':0., 'Omega_q':1.-omega[2], 'Omega_b':omega[1], 'Omega_m':omega[2], 'm_ncdm':0., 'Omega_k':0., 'tau_reio':0.066, 'n_s':omega[3], 'sigma8':omega[4], 'w0_fld':-1., 'wa_fld':0., 'k_max':10.0, 'WhichSpectrum':"EH"}

if not exists(fdir+"phi_obs.npy") or force:
    phi_obs=blackbox.make_data(cosmo_obs,0,force_powerspectrum=force,force_parfiles=force,force_sim=force,force_mock=force,force_cosmo=force)
    np.save(fdir+"phi_obs",phi_obs)
else:
    phi_obs=np.load(fdir+"phi_obs.npy")

if not exists(fdir+"groundtruth.npz") or force:
    groundtruth=PowerSpectrum(L0,L1,L2,N0,N1,N2,cosmo_obs).powerspectrum
    theta_groundtruth_th=np.insert(groundtruth[1:]/P_exp.powerspectrum[1:],0,1.)
    theta_groundtruth=P2theta(get_Pk(k_s,cosmo_obs))
    np.savez(fdir+"groundtruth",groundtruth=groundtruth,theta_groundtruth_th=theta_groundtruth_th,theta_groundtruth=theta_groundtruth)
else:
    A=np.load(fdir+"groundtruth.npz")
    groundtruth=A["groundtruth"]
    theta_groundtruth_th=A["theta_groundtruth_th"]
    theta_groundtruth=A["theta_groundtruth"]

# Optimization of the prior: k-range, hyperpriors, and fiducial power spectrum
k_opt_min=0.
k_opt_max=k_max
theta_norm_mean=0.2
theta_norm_std=0.3
k_corr_mean=0.020
k_corr_std=0.015
omega=planck_mean
cosmo_fid={'h':omega[0], 'Omega_r':0., 'Omega_q':1.-omega[2], 'Omega_b':omega[1], 'Omega_m':omega[2], 'm_ncdm':0., 'Omega_k':0., 'tau_reio':0.066, 'n_s':omega[3], 'sigma8':omega[4], 'w0_fld':-1., 'wa_fld':0., 'k_max':10.0, 'WhichSpectrum':"EH"}
if not exists(fdir+"fiducial.npz") or force:
    fiducial=PowerSpectrum(L0,L1,L2,N0,N1,N2,cosmo_fid).powerspectrum
    theta_fiducial_th=np.insert(fiducial[1:]/P_exp.powerspectrum[1:],0,1.)
    theta_fiducial=P2theta(get_Pk(k_s,cosmo_fid))
    np.savez(fdir+"fiducial",fiducial=fiducial,theta_fiducial_th=theta_fiducial_th,theta_fiducial=theta_fiducial)
else:
    A=np.load(fdir+"fiducial.npz")
    fiducial=A["fiducial"]
    theta_fiducial_th=A["theta_fiducial_th"]
    theta_fiducial=A["theta_fiducial"]

# Output files
fname=fdir+"selfi_SBMY_iter"+str(i_iter)+".h5"
pool_prefix=fdir+"pools_SBMY_iter"+str(i_iter)+"/pool_"
pool_suffix=".h5"

selfi_SBMY=power_spectrum_selfi(fname,pool_prefix,pool_suffix,prior,blackbox,theta_0,Ne,Ns,Delta_theta,phi_obs)
