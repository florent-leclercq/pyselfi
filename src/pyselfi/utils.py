#!/usr/bin/env python
#-------------------------------------------------------------------------------------
# pySELFI v2.0 -- src/pyselfi/utils.py
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

"""Generic infrastructure routines.
"""

__author__  = "Florent Leclercq"
__version__ = "2.0"
__date__    = "2018-2023"
__license__ = "GPLv3"

# Global variables for fonts
FONT_BLACK		="\033[0;30m"
FONT_RED		="\033[0;31m"
FONT_LIGHTRED		="\033[0;91m"
FONT_BURGUNDY		="\033[38;5;161m"
FONT_GREEN		="\033[0;32m"
FONT_LIGHTGREEN		="\033[38;5;113m"
FONT_LIGHTERGREEN	="\033[38;5;192m"
FONT_ORANGE		="\033[38;5;215m"
FONT_YELLOW		="\033[38;5;227m"
FONT_BLUE		="\033[0;34m"
FONT_LIGHTBLUE		="\033[38;5;38m"
FONT_PURPLE		="\033[0;35m"
FONT_LIGHTPURPLE	="\033[38;5;147m"
FONT_PINK		="\033[38;5;207m"
FONT_PYTHONPINK		="\033[38;5;219m"
FONT_CYAN		="\033[0;36m"
FONT_LIGHTCYAN		="\033[38;5;117m"
FONT_GREY		="\033[38;5;246m"

FONT_BOLDBLACK		="\033[1;30m"
FONT_BOLDRED		="\033[1;31m"
FONT_BOLDGREEN		="\033[1;32m"
FONT_BOLDYELLOW		="\033[1;33m"
FONT_BOLDBLUE		="\033[1;34m"
FONT_BOLDPURPLE		="\033[1;35m"
FONT_BOLDCYAN		="\033[1;36m"
FONT_BOLDGREY		="\033[1;37m"

FONT_SIMBELMYNE		="\033[1;38;5;157m"
FONT_NORMAL		="\033[00m"

# Global variables for verbosity
COMMAND_VERBOSITY	=1
ERROR_VERBOSITY		=1
INFO_VERBOSITY		=1
MODULE_VERBOSITY	=2
WARNING_VERBOSITY	=3
DEBUG_VERBOSITY		=6
MEMORY_VERBOSITY	=7
SCREEN_VERBOSE_LEVEL	=4
G__ind__		=0

def INDENT():
    """Indents the current level of outputs.
    """
    global G__ind__
    G__ind__+=1
    return G__ind__

def UNINDENT():
    """Unindents the current level of outputs.
    """
    global G__ind__
    G__ind__-=1
    return G__ind__

def PrintLeftType(message_type, FONT_COLOR):
    """Prints the type of output to screen.

    Parameters
    ----------
    message_type : :obj:`str`
        type of message
    FONT_COLOR : :obj:`str`
        font color for this type of message

    """
    from time import localtime, strftime
    import sys
    sys.stdout.write("["+strftime("%H:%M:%S", localtime())+"|"+FONT_COLOR+message_type+FONT_NORMAL+"]")
    for ind in range(G__ind__):
        sys.stdout.write("==")
    sys.stdout.write("|")

def PrintMessage(verbosity, message):
    """Prints a message to screen.

    Parameters
    ----------
    verbosity : int
        verbosity of the message
    message : :obj:`str`
        message

    """
    if(SCREEN_VERBOSE_LEVEL>=verbosity):
        PrintLeftType("STATUS    ", FONT_LIGHTGREEN)
        import sys
        sys.stdout.write("{}\n".format(message))
        sys.stdout.flush()

def PrintInfo(message):
    """Prints an information to screen.

    Parameters
    ----------
    message : :obj:`str`
        message

    """
    if(SCREEN_VERBOSE_LEVEL>=INFO_VERBOSITY):
        PrintLeftType("INFO      ", FONT_LIGHTCYAN)
        import sys
        sys.stdout.write("{}\n".format(message))
        sys.stdout.flush()

def PrintDiagnostic(verbosity, message):
    """Prints a diagnostic to screen.

    Parameters
    ----------
    verbosity : int
        verbosity of the message
    message : :obj:`str`
        message

    """
    if(SCREEN_VERBOSE_LEVEL>=verbosity):
        PrintLeftType("DIAGNOSTIC", FONT_GREY)
        import sys
        sys.stdout.write(FONT_GREY+message+FONT_NORMAL+"\n")

def PrintValue(name, value):
    """Prints the value of a variable to screen.

    Parameters
    ----------
    name : :obj:`str`
        name of variable to be printed
    value : :obj:`str`
        value of variable to be printed

    """
    PrintDiagnostic(DEBUG_VERBOSITY, "{}={}".format(name,value))

def PrintCommandLine(message):
    """Prints a bash command line to screen.

    Parameters
    ----------
    message : :obj:`str`
        message

    """
    if(SCREEN_VERBOSE_LEVEL>=COMMAND_VERBOSITY):
        PrintLeftType("COMMAND   ", FONT_YELLOW)
        import sys
        sys.stdout.write(FONT_YELLOW+message+FONT_NORMAL+"\n")

def PrintPythonCommand():
    """Prints a python command to screen.

    Parameters
    ----------
    message : :obj:`str`
        message

    """
    from sys import argv
    message='python '
    for i in range(len(argv)):
        message += argv[i]
        if i<len(argv)-1:
            message += ' '
    PrintCommandLine(message)

def ExecuteBashCommand(commandString):
    """Prints and executes a bash command.

    Parameters
    ----------
    commandString : :obj:`str`
        command line

    """
    from os import system
    PrintCommandLine(commandString)
    system(commandString)

def ExecuteBashCommandMute(commandString):
    """Executes a bash command without printing it.

    Parameters
    ----------
    commandString : :obj:`str`
        command line

    """
    from os import system
    system(commandString)

def ExecuteBashCommandOutput(commandString):
    """Prints and executes a bash command, returns the output.

    Parameters
    ----------
    commandString : :obj:`str`
        command line

    Returns
    -------
    output : :obj:`str`
        the result of subprocess.check_output

    """
    from subprocess import check_output
    PrintCommandLine(commandString)
    return check_output(commandString.split(" "), universal_newlines=True)

def PrintWarning(message):
    """Prints a warning to screen.

    Parameters
    ----------
    message : :obj:`str`
        message

    """
    if(SCREEN_VERBOSE_LEVEL>=WARNING_VERBOSITY):
        PrintLeftType("WARNING   ", FONT_ORANGE)
        import sys
        sys.stdout.write(FONT_ORANGE+message+FONT_NORMAL+"\n")

def PrintError(message):
    """Prints an error to screen.

    Parameters
    ----------
    message : :obj:`str`
        message

    """
    if(SCREEN_VERBOSE_LEVEL>=ERROR_VERBOSITY):
        PrintLeftType("ERROR     ", FONT_LIGHTRED)
        import sys
        sys.stdout.write(FONT_LIGHTRED+message+FONT_NORMAL+"\n")

def FatalError(message):
    """Prints an error to screen and ends the script.

    Parameters
    ----------
    message : :obj:`str`
        message

    """
    import sys
    PrintError(message)
    sys.exit()

def p_load(file, mmap_mode=None, allow_pickle=True, fix_imports=True, encoding='ASCII'):
    """Loads arrays or pickled objects from ``.npy``, ``.npz`` or pickled files. (Replacement of numpy.load, see its documentation).
    """
    import numpy as np
    PrintMessage(3, "Loading array from file '{}'...".format(file))
    arr = np.load(file, mmap_mode=mmap_mode, allow_pickle=allow_pickle, fix_imports=fix_imports, encoding=encoding)
    PrintMessage(3, "Loading array from file '{}' done.".format(file))
    return arr

def p_loadtxt(file, dtype=None, comments='#', delimiter=None, converters=None, skiprows=0, usecols=None, unpack=True, ndmin=0, encoding='bytes'):
    """Loads data from a text file. Each row in the text file must have the same number of values. (Replacement of numpy.loadtxt, see its documentation). WARNING: contrary to numpy.loadtxt, unpack is set to True by default.
    """
    import numpy as np
    dtype=dtype or np.float
    PrintMessage(3, "Loading array from file '{}'...".format(file))
    arr = np.loadtxt(file, dtype=dtype, comments=comments, delimiter=delimiter, converters=converters, skiprows=skiprows, usecols=usecols, unpack=unpack, ndmin=ndmin, encoding=encoding)
    PrintMessage(3, "Loading array from file '{}' done.".format(file))
    return arr

def p_save(file, arr, allow_pickle=True, fix_imports=True):
    """Saves an array to a binary file in NumPy ``.npy`` format. (Replacement of numpy.save, see its documentation).
    """
    import numpy as np
    file = file.replace(".npy", "")
    PrintMessage(3, "Saving array to file '{}.npy'...".format(file))
    np.save(file,arr)
    PrintMessage(3, "Saving array to file '{}.npy' done.".format(file))
    
def p_savez(file,*args,**kwds):
    """Saves several arrays into a single file in uncompressed ``.npz`` format. (Replacement of numpy.savez, see its documentation).
    """
    import numpy as np
    file = file.replace(".npz", "")
    PrintMessage(3, "Saving arrays to file '{}.npz'...".format(file))
    np.savez(file,*args,**kwds)
    PrintMessage(3, "Saving arrays to file '{}.npz' done.".format(file))
    
def p_savetxt(file, X, fmt='%.18e', delimiter=' ', newline='\n', header='', footer='', comments='# ', encoding=None):
    """Saves an array to a text file. (Replacement of numpy.savetxt, see its documentation).
    """
    import numpy as np
    PrintMessage(3, "Saving array to file '{}'...".format(file))
    np.savetxt(file, X, fmt=fmt, delimiter=delimiter, newline=newline, header=header, footer=footer, comments=comments, encoding=encoding)
    PrintMessage(3, "Saving array to file '{}' done.".format(file))
    
def save_replace_dataset(hf, address, data, maxshape, dtype):
    """Saves or replaces a dataset in a hdf5 file. The dataset must be resizable (see the HDF5 manual).

    Parameters
    ----------
    hf : :obj:`h5py.File`
        hdf5 file (already opened with writing permission)
    address : :obj:`str`
        address of the dataset to be written/replaced
    data : array
        data to be saved
    maxshape : array or None
        maximum size of the array
    dtype : array
        data type

    """
    if address in hf:
        chunk=hf[address]
        chunk.resize(data.shape)
        chunk[...]=data
    else:
        hf.create_dataset(address, data=data, maxshape=maxshape, dtype=dtype)

def save_replace_attr(hf, address, attr, dtype):
    """Saves or replaces an attribute in a hdf5 file.

    Parameters
    ----------
    hf : :obj:`h5py.File`
        hdf5 file (already opened with writing permission)
    address : :obj:`str`
        address of the attribute to be written/replaced
    attr : array
        attribute to be saved
    dtype : array
        data type

    """
    if address in hf.attrs:
        hf.attrs[address]=attr
    else:
        hf.attrs.create(address, attr, dtype=dtype)

def get_indices(Id, N0, N1=None, N2=None):
    """Gets 3D indices from 1D index. Assumes row-major ordering.

    Parameters
    ----------
    Id : int
        1D index
    N0 : int
        size of the array x
    N1 : int, optional, default=N0
        size of the array y
    N2 : int, optional, default=N0
        size of the array z

    Returns
    (i,j,k) : array, int, dimension=3
        corresponding 3D indices

    """
    N1=N1 or N0
    N2=N2 or N0
    i = int((Id/(N1*N2))%N0)
    j = int(((Id-N1*N2*i)/N2)%N1)
    k = int((Id-N2*j-N2*N1*i)%N2)
    return (i,j,k)

def get_index(i,j,k,N1,N2=None,N0=None):
    """Gets 1D index from 3D indices. Assumes row-major ordering. WARNING: N0 is passed last when all of N0,N1,N2 are specified!

    Parameters
    ----------
    i : int
        3D index x
    j : int
        3D index y
    k : int
        3D index z
    N1 : int
        size of the array y
    N2 : int, optional, default=N1
        size of the array z
    N0 : int, optional, default=N1
        size of the array x

    Returns
    -------
    Id : int
        corresponding 1D index

    """
    N2=N2 or N1
    return k+N2*(j+N1*i);

def regular_inv(A, EPS_A=1e-7, EPS_residual=1e-3):
    """Computes the inverse of a matrix. Attempts a regularization by adding an epsilon to the diagonal, if the matrix is ill-conditioned.

    Parameters
    ----------
    A : array-like, square
        input matrix
    EPS_A : double, optional, default=1e-7
        epsilon to be added to the diagonal of A, if necessary
    EPS_residual : double, optional, default=1e-3
        maximum residual in A*A^{-1} before attempting regularization

    Returns
    -------
    A^{-1} : array-like, square
        numerical inverse of A

    """
    import numpy as np
    import scipy.linalg as sla
    
    inv_A = sla.inv(A)
    residual=np.sum(np.fabs(inv_A.dot(A)-np.identity(A.shape[0])))
    if residual>EPS_residual:
        A=A+EPS_A*np.identity(A.shape[0])
        inv_A = sla.inv(A)
        residual=np.sum(np.fabs(inv_A.dot(A)-np.identity(A.shape[0])))
        if residual>1e-3:
            from pyselfi.utils import PrintValue
            PrintValue("condition number of A", np.linalg.cond(A))
            PrintValue("residuals in A*A^{-1}", residual)
            raise ValueError("Matrix inversion is imprecise for A, even after regularization attempt.")
        inv_A = inv_A-EPS_A*np.identity(A.shape[0])
    return inv_A
