#!/usr/bin/env python
#-------------------------------------------------------------------------------------
# pySELFI v1.1 -- pyselfi_examples/templates/new_selfi.py
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

"""Template for a new "selfi" for pySELFI
"""

__author__  = "Florent Leclercq"
__version__ = "1.1"
__date__    = "2018-2019"
__license__ = "GPLv3"

from pyselfi import selfi

class child_selfi(selfi):
    """This class represents a child "selfi" object.
    """
    
    def additional_method(self, **kwargs):
        """An additional method for "selfi" objects

        Parameters
        ----------
        **kwargs : dictionary, optional
            any argument


        """
        return result
    
# end class(child_selfi)
