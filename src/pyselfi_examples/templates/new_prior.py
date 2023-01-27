#!/usr/bin/env python
#-------------------------------------------------------------------------------------
# pySELFI v2.0 -- src/pyselfi_examples/templates/new_prior.py
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

"""Template for a new prior for pySELFI
"""

__author__  = "Florent Leclercq"
__version__ = "2.0"
__date__    = "2018-2023"
__license__ = "GPLv3"

class new_prior(object):
    """This class represents a SELFI prior.

    Attributes
    ----------
    **kwargs : dictionary, optional
        any prior dependency, such as support wavenumbers/multipoles, or statistical hyperparameters

    """

    # Initialization
    def __init__(self,kwargs):
        """Initializes a prior object.
        """
        for key, value in kwargs.items():
            setattr(self,key,value)

    @property
    def some_property(self):
        """type: A property of the prior
        """

        # DO SOMETHING
        return result

    def _private_method(self, **args):
        """A private method of the prior class.
        """

        # DO SOMETHING
        return result

    def compute(self):
        """Computes the prior (mean, covariance matrix and its inverse).
        """
        # Set the mandatory attributes, if not they have not already been set by the constructor:
        self.mean = self._get_mean()
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
            # Save interesting hyperparameters as follows, if needed
            save_replace_attr(hf, '/prior/some_hyperparameter', self.some_hyperparameter, dtype=c_double)

            # Save here any other interesting prior attribute as follows, if needed
            save_replace_dataset(hf, '/prior/some_array', data=self.some_array, maxshape=(None,), dtype=c_double)

            # Save mandatory attributes: mean, covariance and inv_covariance
            save_replace_dataset(hf, '/prior/mean', data=self.mean, maxshape=(None,), dtype=c_double)
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
            # Load here any mandatory parameter to call the class constructor as follows, if needed
            some_hyperparameter=hf.attrs['/prior/some_hyperparameter']
            some_array=np.array(hf.get('/prior/some_array'), dtype=c_double)

            # Call the class constructor
            prior=cls(kwargs) #kwargs can depend on some_hyperparameter, some_array, etc.

            # Load any other parameter, in particular the mandatory attributes (mean, covariance and inv_covariance), if they are not already defined at this point
            prior.mean = np.array(hf.get('prior/mean'), dtype=c_double)
            prior.covariance = np.array(hf.get('/prior/covariance'), dtype=c_double)
            prior.inv_covariance = np.array(hf.get('/prior/inv_covariance'), dtype=c_double)

        PrintMessage(3, "Reading prior in data file '{}' done.".format(fname))
        return prior
#end class(prior)
