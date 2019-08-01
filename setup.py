#!/usr/bin/env python
#-------------------------------------------------------------------------------------
# pySELFI v1.2 -- setup.py
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

"""Setup script for pyselfi
"""

__author__  = "Florent Leclercq"
__version__ = "1.2"
__date__    = "2018-2019"
__license__ = "GPLv3"

import setuptools
import io

# Open README.md
def readme():
    with io.open('README.md', mode='r', encoding='ISO-8859-1') as f:
        return f.read()

# Open requirements.txt
def requirements():
    with io.open('requirements.txt', mode='r') as f:
        return f.read().splitlines()

optionals = {'doc': ['Sphinx']}

# Setup
setuptools.setup(
    name="pyselfi",
    version="1.2",
    author="Florent Leclercq",
    author_email="florent.leclercq@polytechnique.org",
    description="A python implementation of the Simulator Expansion for Likelihood-Free Inference (SELFI) algorithm",
    long_description=readme(),
    long_description_content_type='text/markdown',
    url="http://pyselfi.florent-leclercq.eu",
    packages=setuptools.find_packages(),
    install_requires=requirements(),
    extras_require=optionals,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
)
