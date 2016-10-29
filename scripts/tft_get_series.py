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
:mod:`tft_get_series.py`
=========================
.. moduleauthor:: Joe Kwiatkowski

Return a series of data from a ToFeT output file.

Command line usage
-------------------

.. code-block:: bash

tft_get_series.py  TOFET_OUTPUT  SERIES

where *SERIES* is one of the following options:
    * occupation_probabilities (for occupation probabilities)
    * transient (for ToF photocurrent transients)
    * occupied_molecules (for molecules occupied at the end of the simulation)
    * energies (for the static energy + average Coulomb energy of each molecule)
"""

import sys
from math import *
from optparse import OptionParser

option_string = {"occupation_probabilities": "> TOTAL OCCUPATION TIMES AND TIMES VISITED",
                 "transient": "> PHOTOCURRENT TRANSIENT",
                 "occupied_molecules": "> OCCUPIED MOLECULES AT END OF SIMULATION",
                 "energies": "> ENERGY (static + Coulomb)",
                 
                }

def main(filename, option):
    parser = OptionParser()
    parser.add_option("--occ", help="Extract occupation data.")
    parser.add_option("--trans", help="Extract photocurrent transient data.")
    (options, args) = parser.parse_args(sys.argv[1:])

    exit()

    try: 
        string = option_string[option]
    except KeyError:
        print "*** ERROR *** : Don't understand option " + option
        return 1
    inFile = open(filename)
    data=[]
    for line in inFile:
        if string in line: 
            break
    else:
        print "*** WARNING ***: Couldn't find " + option + " in " + filename
        return 1
    inFile.next()   #Read in column labels
    for line in inFile:
        if ">" in line: break
        else: data.append(line.strip())
    return data

if __name__ == '__main__':
    data = main(sys.argv[1], sys.argv[2])
    if isinstance(data, list):
        for i in range(len(data)): 
            print data[i] 
