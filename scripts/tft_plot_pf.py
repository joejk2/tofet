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
:mod:`tft_plot_pf.py`
=====================

.. moduleauthor:: Joe Kwiatkowski

From ToFeT output files, extract the data necessary to plot Poole-Frenkel
plots::

    sqrt[F (V/cm)]      mobility (cm^2/V.s)
        
Notice the units of the field are in (V/cm): this makes the experimentalists
happy.

Command-line usage
--------------------

.. code-block:: bash
        
        tft_plot_pf.py  FILES

where *FILES* are specified by a :term:`glob pattern`. 

"""

import sys, os, fnmatch, math


def plot_pf(filename):
    filenames = sys.argv[1:]
    mobilities = []
    sqrt_fields  = []
    for filename in filenames:
        sqrt_field, mobility = -1, -1
        inFile = open(filename)
        for line in inFile:
            if ": fieldZ" in line:
                sqrt_field = math.sqrt( - float (line.split()[2]) )
            if "> MOBILITY FROM COLLECTION TIMES (cm^2/V.s)=" in line:
                mobility = line.split()[6]
            if (sqrt_field!=-1) and (mobility!=-1):  
                sqrt_fields.append(sqrt_field)
                mobilities.append(mobility)
                break;
        if (sqrt_field == -1) or (mobility == -1):  
            print "*** ERROR *** : Didn't find the field or mobility in ", filename
            exit()

    for sqrt_field, mobility in sorted ( zip(sqrt_fields, mobilities) ):
        print "%e\t\t%s" % (sqrt_field*1e4, mobility) #V/cm!


if __name__ == '__main__':
    plot_pf(sys.argv[1])
