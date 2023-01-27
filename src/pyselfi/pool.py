#!/usr/bin/env python
#-------------------------------------------------------------------------------------
# pySELFI v2.0 -- src/pyselfi/pool.py
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

"""Routines related to simulation pools.
"""

__author__  = "Florent Leclercq"
__version__ = "2.0"
__date__    = "2018-2023"
__license__ = "GPLv3"

class pool(object):
    """This class represents a pool of simulations at a given point in parameter space.

    Attributes
    ----------
    fname : :obj:`str`
        filename of the pool to load or create
    N_target : int, optional
        number of simulations desired in the pool

    """


    # Initialization
    def __init__(self,fname,N_target=None):
        """Initializes a pool object.
        """
        import h5py as h5
        import numpy as np
        from ctypes import c_int, c_float
        from pyselfi.utils import PrintMessage
        from pathlib import Path

        f = Path(fname)
        if f.is_file():
            PrintMessage(3, "Loading pool from data file '{}'...".format(fname))

            with h5.File(fname,'r') as hf:
                if N_target is None:
                    N_target = hf.attrs['/info/scalars/N_target']
                if '/scalars/Phi' in hf:
                    Phi = np.array(hf.get('/scalars/Phi'), dtype=c_float).astype(c_float)
                else:
                    Phi = None

            PrintMessage(3, "Loading pool from data file '{}' done.".format(fname))
        else:
            PrintMessage(3, "Initializing pool in data file '{}'...".format(fname))

            with h5.File(fname, 'w') as hf:
                hf.attrs.create('/info/scalars/N_target', N_target, dtype=c_int)
                hf.create_dataset('/scalars/Phi', data=np.empty((0,0)), maxshape=(None, None), dtype=c_float)
            Phi = None

            PrintMessage(3, "Initializing pool in data file '{}' done.".format(fname))

        self.fname=fname
        self.N_target=N_target
        self.Phi=Phi

    @property
    def N_sims(self):
        """int: Number of simulations currently available in the pool
        """
        return self.Phi.shape[0] if self.Phi is not None else 0

    @property
    def finished(self):
        """bool: Whether the number of available simulations is larger than the target number
        """
        return self.N_sims>=self.N_target if self.N_target is not None else True

    def save(self):
        """Saves the pool to its output file.
        """
        import h5py as h5
        from ctypes import c_float
        from pyselfi.utils import PrintMessage
        PrintMessage(3, "Save pool to data file '{}'...".format(self.fname))

        with h5.File(self.fname, 'r+') as hf:
            data=hf['/scalars/Phi']
            data.resize(self.Phi.shape)
            data[...]=self.Phi
            hf.attrs['/info/scalars/N_target']=self.N_target

        PrintMessage(3, "Save pool to data file '{}' done.".format(self.fname))

    def add_sim(self,phi):
        """Adds a simulation to the pool.

        Parameters
        ----------
        phi : array, double, dimension=P
            summaries of the simulation to be added ot the pool

        """
        import numpy as np
        if self.N_sims==0:
            self.Phi=np.atleast_2d(phi)
        else:
            self.Phi=np.append(self.Phi,np.atleast_2d(phi),axis=0)
#end class(pool)
