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
.. moduleauthor:: Joe Kwiatkowski

Take a file containing *(X, Y, Z)* coordinates, and define collection and
generation regions.

Command-line usage
--------------------

.. code-block:: bash

    tft_define_types.py  XYZ  COL_GEN_FRACTION  [OPTIONS]


where :ref:`XYZ <sec_xyz_file>` is the file containing *(X, Y, Z)* data in three columns, 
and :term:`COL_GEN` is the fraction of the film along the *Z*-axis which should
generate and collect [#f1]_.

Options
------------

.. program:: tft_define_types.py

.. cmdoption:: -e

    This option specifies that site energies are given in the 
    :ref:`xyz file <sec_xyz_file>`. 


.. [#f1] 
    So the total fraction of the film that is a 
    generator *or* a collector is 2 * :term:`COL_GEN`. 

"""

import sys
import numpy
from optparse import OptionParser
from io import print_lists_to_columns, read_file_to_lines


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    usage = "usage: %prog [options] foo.xyz collection_generation_cg_fraction"
    option_parser = OptionParser(usage=usage)
    option_parser.add_option("-e", default=False, action="store_true",\
                         help="If specified, energies expected in the "\
                         "xyz file.")
    (options, args) = option_parser.parse_args(argv)
    if (len(args) != 2):
        option_parser.print_help()
        exit()
    xyz_file = open(args[0])
    cg_fraction = float(args[1])
    if (cg_fraction>0.5) or (cg_fraction<0.0):
        print "Oops: (0.0 < cg_fraction < 0.5) Please!"
        exit()
    #Read in coords
    xyz_lines = read_file_to_lines(args[0])
    x, y, z, E = [], [], [], []
    for line in xyz_lines:
        words = line.split()
        x.append(float(words[0]))
        y.append(float(words[1]))
        z.append(float(words[2]))
        # If necessary, read in energies.
        if options.e:
            E.append(float(words[3]))
    if len(x) == len(y) == len(z):
        pass
    else:
        print where_am_i()[0] + ": badly formatted input in", args[0]
        exit()
    #Define and output collectors and generators
    z_array = numpy.array(z)
    zMax = max(z_array)
    type = []
    for i in range(len(x)):
        if (z[i] <= cg_fraction * zMax):   
            type.append('g')
        elif (z[i] >= (1 - cg_fraction) * zMax): 
            type.append('c')
        else:
            type.append('-')
    if options.e:
        print_lists_to_columns(x, y, z, type, E)
    else:
        print_lists_to_columns(x, y, z, type)
    
    
if __name__ == '__main__':
    main()
