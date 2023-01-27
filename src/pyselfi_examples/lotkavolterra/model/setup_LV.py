#!/usr/bin/env python
#-------------------------------------------------------------------------------------
# pySELFI v2.0 -- src/pyselfi_examples/lotkavolterra/model/setup_LV.py
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

"""Setup script for the pySELFI Gaussian random field example
"""

__author__  = "Florent Leclercq"
__version__ = "2.0"
__date__    = "2018-2023"
__license__ = "GPLv3"

import numpy as np
import scipy.stats as ss
import h5py as h5
from os.path import dirname, realpath, exists
from lotkavolterra_simulator import LVsimulator
from pyselfi_examples.lotkavolterra.model.blackbox_LV import blackbox
from pyselfi.lotkavolterra.prior import lotkavolterra_prior
from pyselfi.lotkavolterra.selfi import lotkavolterra_selfi

# Basic configuration
fdir=realpath(dirname(realpath(__file__))+"/../sims/")+"/"
tmax=50
tres=100
t = np.linspace(0,tmax-1,tmax*tres+1)
t_s = np.arange(tmax)
fixnoise=False
seed=1

# ------ THEORY ---------
# Groundtruth values
Delta_t = 2.
alpha_gt, beta_gt, gamma_gt, delta_gt = 1/Delta_t*np.array([1.1, 0.4, 0.4, 0.1])
X0, Y0 = 10, 5
gt = np.array([alpha_gt, beta_gt, gamma_gt, delta_gt])

# Define simulator instance
LVt = LVsimulator(X0, Y0, alpha_gt, beta_gt, gamma_gt, delta_gt)

# Solve the ODEs
Xtheory, Ytheory = LVt.EEuler(t)

# ------- DATA ---------
# Simulate true data
Xtrue, Ytrue = LVt.get_XY(t_s, tres, Xtheory, Ytheory)
theta_true = LVt.get_theta(Xtrue, Ytrue)

# Define a mask
masks=[np.array([4, 5, 6]), np.array([19, 20]), np.array([35, 36, 37])]
mask=np.array([0 if n in np.concatenate(masks) else 1 for n in range(tmax)])

# Define blackbox instances
S=2*tmax
PX=np.sum(mask)
PY=len(Ytrue)
P=PX+PY
mask=mask
# Correct data model. Parameters controling the non-linearity of the model are obsP, obsQ, obsT
threshold=11
obsP=0.05
obsQ=0.01
obsR0=0.15
obsS=0.05
obsT=0.2
def X_efficiency(t):
    return 1/2.*(1+np.cos(t/(1.7*Delta_t)-1/3*np.pi))
def Y_efficiency(t):
    return 1.
Xefficiency0 = X_efficiency(t_s)
Yefficiency0 = Y_efficiency(t_s)
bbA=blackbox(P, X0=X0, Y0=Y0, seed=seed, fixnoise=fixnoise, model=0, Xefficiency=Xefficiency0, Yefficiency=Yefficiency0, threshold=threshold, mask=mask, obsP=obsP, obsQ=obsQ, obsR=obsR0, obsS=obsS, obsT=obsT)

# Misspecified data model
obsR1=70/100.*obsR0
bbB=blackbox(P, X0=X0, Y0=Y0, seed=seed, fixnoise=fixnoise, model=1, mask=mask, obsR=obsR1)

# Simulate observed data
Xobs, Yobs = bbA.Xobs(S,mask), bbA.Yobs(S)
phi_obs = bbA.make_data(theta_true, seed=seed)
Xdata, Ydata = phi_obs[0:PX], phi_obs[PX:P]

# Check the model for the signal
theta_tsignal = bbA.make_signal(theta_true)
Xtsignal, Ytsignal = theta_tsignal[0:S//2], theta_tsignal[S//2:S]

# Check model for the noise covariance matrix
Dtsignal=bbA.make_demographic_noise_cov(theta_true)
Otsignal=bbA.make_observational_noise_cov(theta_true)

# ------ PRIOR ----------
# Expansion values
#alpha_exp, beta_exp, gamma_exp, delta_exp = alpha_gt, beta_gt, gamma_gt, delta_gt
exp_deviation=3/100.
alpha_exp, beta_exp, gamma_exp, delta_exp = ss.multivariate_normal(mean=gt,cov=(gt*exp_deviation)**2*np.identity(4), seed=seed).rvs()
exp = np.array([alpha_exp, beta_exp, gamma_exp, delta_exp])

# Define simulator instance
LVe = LVsimulator(X0, Y0, alpha_exp, beta_exp, gamma_exp, delta_exp)

# Solve the ODEs for the expansion values to find the expansion function and fiducial function
Xthexp, Ythexp = LVe.EEuler(t)

# Find values of the expansion function and fiducial function
Xexp, Yexp = LVe.get_XY(t_s, tres, Xthexp, Ythexp)
theta_exp = LVe.get_theta(Xexp, Yexp)

# Define observer instance and compute signal at the expansion point
theta_esignal = bbA.make_signal(theta_exp)
Xesignal, Yesignal = theta_esignal[0:S//2], theta_esignal[S//2:S]

# Define and compute the prior
alpha_norm = 0.02
t_smooth = 1.6
t_chaos = 8.2
priorA=lotkavolterra_prior(t_s,X0,Y0,theta_exp,alpha_norm,t_smooth,t_chaos)
priorB=lotkavolterra_prior(t_s,X0,Y0,theta_exp,alpha_norm,t_smooth,t_chaos)

# Simulate an ensemble of fiducial functions drawn around the expansion point
Nfid=100
Xthfid, Ythfid = np.zeros((Nfid,len(t))), np.zeros((Nfid,len(t)))
Xfid, Yfid = np.zeros((Nfid,S//2)), np.zeros((Nfid,S//2))
theta_fid = np.zeros((Nfid,S))
for n in range(Nfid):
    alpha_fid, beta_fid, gamma_fid, delta_fid = ss.multivariate_normal(mean=exp,cov=(exp*exp_deviation)**2*np.identity(4)).rvs()
    fid = np.array([alpha_fid, beta_fid, gamma_fid, delta_fid])
    LVf = LVsimulator(X0, Y0, alpha_fid, beta_fid, gamma_fid, delta_fid)
    Xthfid[n,:], Ythfid[n,:] = LVf.EEuler(t)
    Xfid[n,:], Yfid[n,:] = LVf.get_XY(t_s, tres, Xthfid[n,:], Ythfid[n,:])
    theta_fid[n,:] = LVf.get_theta(Xfid[n,:], Yfid[n,:])

# ---- LIKELIHOOD -------
# Finite differencing setup
Ns=100
Ne=150
Delta_theta=1e-2

# Output files
fnameA=fdir+"selfi_LV_modelA.h5"
pool_prefixA=fdir+"pools_LV_modelA/pool_"
fnameB=fdir+"selfi_LV_modelB.h5"
pool_prefixB=fdir+"pools_LV_modelB/pool_"
pool_suffix=".h5"

selfi_LV_modelA=lotkavolterra_selfi(fnameA,pool_prefixA,pool_suffix,priorA,bbA,theta_exp,Ne,Ns,Delta_theta,phi_obs,LVe)
selfi_LV_modelB=lotkavolterra_selfi(fnameB,pool_prefixB,pool_suffix,priorB,bbB,theta_exp,Ne,Ns,Delta_theta,phi_obs,LVe)

# ----- INFERENCE -------
inference_prior_deviation = exp_deviation
alphamin, alphamax = alpha_exp*(1.-3.5*inference_prior_deviation), alpha_exp*(1+3.5*inference_prior_deviation)
betamin, betamax = beta_exp*(1.-3.5*inference_prior_deviation), beta_exp*(1+3.5*inference_prior_deviation)
gammamin, gammamax = gamma_exp*(1.-3.5*inference_prior_deviation), gamma_exp*(1+3.5*inference_prior_deviation)
deltamin, deltamax = delta_exp*(1.-3.5*inference_prior_deviation), delta_exp*(1+3.5*inference_prior_deviation)
