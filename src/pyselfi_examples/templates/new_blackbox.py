#!/usr/bin/env python
#-------------------------------------------------------------------------------------
# pySELFI v2.0 -- src/pyselfi_examples/templates/new_blackbox.py
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

"""Template for a new blackbox for pySELFI
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
    **kwargs : dictionary, optional
        any other argument

    """

    # Initialization
    def __init__(self, P, **kwargs):
        """Initializes the blackbox object.
        """
        self.P=P
        for key, value in kwargs.items():
            setattr(self,key,value)

    def _private_method(self, **args):
        """A private method of the blackbox class.
        """

        # DO SOMETHING
        return result

    def make_auxiliary_needed_product(self, **args):
        """Produces any auxiliary product needed, such as a survey geometry file
        """

        #DO SOMETHING...

    def make_data(self, **args):
        """Evaluates the blackbox to make mock data.

        Returns
        -------
        Phi : array, double, dimension=P
            vector of summary statistics

        """

        # DO SOMETHING
        Phi = self._private_method(args)
        # DO SOMETHING

        return Phi

    def evaluate(self, theta, **args):
        """Evaluates the blackbox from an input vector of parameters.

        Parameters
        ----------
        theta : array, double, dimension=S
            vector of input parameters

        Returns
        -------
        Phi : array, double, dimension=P
            vector of summary statistics

        """

        # DO SOMETHING
        Phi = self._private_method(args)
        # DO SOMETHING

        return Phi

    def compute_pool(self, theta, d, pool_fname, N, **args):
        """Computes a pool of realizations of the blackbox. A method compute_pool with this prototype is the only mandatory method.

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
        from pyselfi.pool import pool
        p=pool(pool_fname,N)

        # Run N evaluations of the blackbox at the desired point
        while not p.finished:
            i = p.N_sims+1
            Phi = self.evaluate(theta, args) # args will contain d,i,N, etc.
            p.add_sim(Phi)
            # It can be convenient to define a variable save_frequency, otherwise one can just call p.save() after every blackbox evaluation
            if i%self.save_frequency==0:
                p.save()
        p.save()

        return p
#end class(blackbox)
