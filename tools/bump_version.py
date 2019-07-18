#!/usr/bin/env python

"""Script to automatically bump a pySELFI version
"""

__author__  = "Florent Leclercq"
__version__ = "SELFI_VERSION"
__date__    = "2018-SELFI_YEAR"
__license__ = "GPLv3"

# Define the files to be bumped
python_paths = ["*.py", "pyselfi/*.py", "pyselfi/power_spectrum/*.py", "examples/grf/model/*.py", "examples/simbelmyne/model/*.py"]

import argparse
parser = argparse.ArgumentParser(description="Bumps a version of pySELFI for release.")
parser.add_argument("version", help="Version number")
args = parser.parse_args()

def bump_version(version):
    """Bumps pySELFI source files to a given version

    Parameters
    ----------
    version (string) : the version number, in the format
    {major}.{minor}.{release}

    """
    def replace_in_file(filepath,text_to_search,replacement_text):
        import fileinput
        result=False
        with fileinput.FileInput(filepath, inplace=True) as file:
            for line in file:
                if text_to_search in line:
                    result=True
                print(line.replace(text_to_search, replacement_text, 1), end='')
        return result

    def get_current_year():
        import datetime
        now = datetime.datetime.now()
        return now.year

    def bump_version_year(filepath,filename,version,year):
        a = replace_in_file(filepath, "{}-SELFI_YEAR".format(str(year)), str(year))
        b = replace_in_file(filepath, "SELFI_YEAR", str(year))
        if not a and not b:
            from pyselfi.utils import PrintWarning
            PrintWarning("Couldn't bump file '{}': couldn't find 'SELFI_YEAR'".format(filename))
        c = replace_in_file(filepath, "SELFI_VERSION", version)
        if not c:
            from pyselfi.utils import PrintWarning
            PrintWarning("Couldn't bump file '{}': couldn't find 'SELFI_VERSION'".format(filename))
    
    def line_prepender(filepath, line):
        result=False
        with open(filepath, 'r+') as f:
            content = f.read()
            f.seek(0, 0)
            f.write(line.rstrip('\r\n') + '\n' + content)
            result=True
        return result
        
    def add_header(filepath,filename,filetype,version,year):
        if filetype=="python":
            old_header = "#!/usr/bin/env python"
            header = """#!/usr/bin/env python
#-------------------------------------------------------------------------------------
# pySELFI v{} -- {}
# Copyright (C) 2019-{} Florent Leclercq.
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
#-------------------------------------------------------------------------------------""".format(version,filename,str(year))
            if not replace_in_file(filepath, old_header, header):
                from pyselfi.utils import PrintWarning
                PrintWarning("Couldn't bump file '{}': couldn't find '{}' to prepender header".format(filename,old_header))            
    
    def bump_file(f, filetype, version, year):
        from os.path import join, dirname, realpath
        from pyselfi.utils import PrintMessage
        SELFI_ROOT=realpath(dirname(realpath(__file__))+"/../")
        bump_version_year(join(SELFI_ROOT,f), f, version, year)
        add_header(join(SELFI_ROOT,f), f, filetype, version, year)
        PrintMessage(4,"Bumped '{}' to version {}".format(f,version))
    
    def get_filelist(paths):
        import glob
        from os.path import join, relpath, dirname, realpath
        SELFI_ROOT=realpath(dirname(realpath(__file__))+"/../")
        files = []
        for path in paths:
            files.extend(glob.glob(join(SELFI_ROOT,path)))
        files = [relpath(f, SELFI_ROOT) for f in files]
        return files
    
    from pyselfi.utils import PrintInfo

    PrintInfo("This is the pySELFI bump version script.")
    
    year = get_current_year()
    
    # loop on python files
    for f in get_filelist(python_paths):
        bump_file(f, "python", version, year)
    
    PrintInfo("Files modified successfully, version bumped to {}.".format(version))

bump_version(args.version)
