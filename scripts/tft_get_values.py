#!/usr/bin/python
#######################################################################
##  This file is part of ToFeT.
##  
##  ToFeT is free software: you can redistribute it and/or modify
##  it under the terms of the GNU Lesser General Public License as published by
##  the Free Software Foundation, either version 3 of the License, or
##  (at your option) any later version.
##  
##  ToFeT is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU Lesser General Public License for more details.
##  
##  You should have received a copy of the GNU Lesser General Public License
##  along with ToFeT.  If not, see <http://www.gnu.org/licenses/>.
#######################################################################
"""
:mod:`tft_get_values.py`
==========================

.. moduleauthor:: Joe Kwiatkowski

Return data value(s) from ToFeT output file(s).

Command line usage
-------------------

.. code-block:: bash

    tft_get_values.py  TOFET_OUTPUT  SERIES

where *TOFET_OUTPUT* is a :term:`glob pattern` and \
*SERIES* is one of the following options:
    * *fet_current* (for the FET converged current)
    * *Vg* (for the gate voltage)
    * *hoppers* (for the number of hoppers left at the end of the simulation)
    * *temperature* (...)
    * *displacement_mobility* (mobility from the total displacement and total\
        time)
    * *velocity_mobility* (mobility from the collection times)

"""


import sys
from glob import glob

option_strings = {"fet_current": "> CURRENT (A) =", 
                  "Vg": ": Vg ",
                  "hoppers": "> NUMBER OF HOPPERS LEFT =",
                  "temperature": "  : temp",
                  "displacement_mobility": "> MOBILITY FROM TOTAL DISPLACEMENT AND TOTAL TIME (cm^2/V.s)=",
                  "velocity_mobility": "> MOBILITY FROM COLLECTION TIMES (cm^2/V.s)=",
                 }

def getValuesFromMultipleFiles(args, expand=True):
    try: 
        i = args.index("-o")
    except ValueError:
        print "*** ERROR *** : Didn't find '-o' in your options"
        return 1
    if expand == True:
        filenames = glob(args[0])
    else:
        filenames = args[:i]
    options = args[i + 1:]
    retval = []
    for f in filenames:
        retval.append(main([f] + options))
    return retval

def main(args):
    inFile = open(args[0])
    strings = []
    for option in args[1:]:
        if option not in option_strings.keys():
            print "*** ERROR *** : Don't understand " + option
            exit()
        else:
            strings.append(option_strings[option]) 

    data = {}
    for line in inFile:
        for string in strings:
            if string in line:
                data[string] = line.split()[-1]

    if (len(data) != len(args[1:])):
        print "*** WARNING ***: Couldn't find all of " + repr(args[1:]) + " in " + args[0]
        return 1

    retval = []
    for string in strings:
        retval.append(data[string])
    return retval

if __name__ == '__main__':
    data = getValuesFromMultipleFiles(sys.argv[1:], False)
    if isinstance(data, int):
        exit()
    if isinstance(data[0], list):
        for row in data:
            for item in row:
                print item,
            print ""

