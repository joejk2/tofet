#!/usr/bin/python
#######################################################################
##  This file is part of ToFeT.
##  
##  ToFeT is free software: you can redistribute it and/or modify
##  it under the terms of the GNU Lesser General Public License as 
##  published by the Free Software Foundation, either version 3 of 
##  the License, or (at your option) any later version.
##  
##  ToFeT is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU Lesser General Public License for more details.
##  
##  You should have received a copy of the GNU Lesser General Public 
##  License  along with ToFeT.  If not, 
##  see <http://www.gnu.org/licenses/>.
#######################################################################
"""
:mod:`tft_extract.py`
======================

.. moduleauthor:: Joe Kwiatkowski

Extract data from a ToFeT output file(s).  Can return either single values, or
entire series of data.  There are options to :mod:`tft_extract.py` to extract
just about any data set from a ToFeT output file.

For full options and usage type::

    tft_extract.py --help

Example use
----------------------

Extract temperature and mobility from two files::

    tft_extract.py  --temp  --mob_vel  OUTPUT_FILE_1 OUTPUT_FILE_2

Extract the photocurrent transient from all files which match the 
glob :file:`ideal_transient*.out`::

    tft_extract.py --transient  ideal_transient*.out


Assumptions made
--------------------
When searching for series of data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
In ToFeT's output file, series of data must have the following form::

    ... UNIQUE_HEADER ...
    ... some description of the columns ...
    row_0
    row_1
    row_2
    .
    .
    .
    > 

Note the final ``>``, which is necessary to tell ToFeT when the series has
finished.   
The *UNIQUE_HEADER* (or enough of it to make it unique)
 must be specified within this module source.

When searching for single values
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Values must have the following form::

    ... UNIQUE_LINE_BEGINNING ... value

Note that the value must be the last element of the line.  
The *UNIQUE_LINE_BEGINNING* (or enough of it to make it unique)
 must be specified within this module source.

"""

import sys
from optparse import OptionParser

# Option - string match for a data *series*.
s_option_string = {}
s_option_string['--occ_prob'] = "> TOTAL OCCUPATION TIMES AND TIMES VISITED"
s_option_string['--transient'] = "> PHOTOCURRENT TRANSIENT"
s_option_string['--occ_mol'] = "> OCCUPIED MOLECULES AT END OF SIMULATION"
s_option_string['--energies'] = "> ENERGY (static + Coulomb)"
# Option - string match for a data *value*.
v_option_string = {}
v_option_string["--i_fet"] = "> CURRENT (A) ="
v_option_string["--vg"] = ": Vg "
v_option_string["--reorg"] = ": reorg "
v_option_string["--fieldZ"] = ": fieldZ "
v_option_string["--hop_start"] = ": hoppers "
v_option_string["--maxTime"] = ": maxTime "
v_option_string["--deltaTime"] = ": deltaTime "
v_option_string["--alpha"] = ": alpha "
v_option_string["--tol"] = ": tol "
v_option_string["--hop_end"] = "> NUMBER OF HOPPERS LEFT ="
v_option_string["--temp"] = " : temp"
v_option_string["--mob_dis"] = "> MOBILITY FROM TOTAL DISPLACEMENT AND TOTAL TIME (cm^2/V.s)="
v_option_string["--mob_vel"] = "> MOBILITY FROM COLLECTION TIMES (cm^2/V.s)="


def extract_value_from_file(filename, string_match):
    input = open(filename) 
    for line in input:
        if string_match in line:
            data = line.split()[-1]  
            break
    else:
        print "*** WARNING ***: " + filename + " doesn't contain:\n " +\
              "                " + string_match 
        return 1
    return data
    

def extract_series_from_file(filename, string_match):
    input = open(filename) 
    for line in input:
        if string_match in line:
            break
    else:
        print "*** WARNING ***: " + filename + " doesn't contain:\n " +\
              "                " + string_match 
        return 1
    input.next()   #Read in column labels
    data = "" 
    for line in input:
        if ">" in line: 
            break
        else: 
            data += line
    return data


def wrap_extract(filename, option):
    if option in s_option_string:
        string_match = s_option_string[option]
        extract_from_file = extract_series_from_file
    elif option in v_option_string:
        string_match = v_option_string[option]
        extract_from_file = extract_value_from_file
    else:
        print "*** ERROR ***: Option " + option + " not understood."
        exit()
    return extract_from_file(filename, string_match)



def optparse_wrapper(option, opt_str, value, parser):
    options.append(opt_str)


if __name__ == '__main__': 
    options = []
    parser = OptionParser("%prog: OPTIONS FILE_1 [FILE_2, FILE_3, ...]")
    for option in s_option_string.keys():
        parser.add_option(option, help=s_option_string[option],
                          action='callback', callback=optparse_wrapper)
    for option in v_option_string.keys():
        parser.add_option(option, help=v_option_string[option],
                          action='callback', callback=optparse_wrapper)

    (parser_options, args) = parser.parse_args(sys.argv[1:])
    if len(args) < 1:
        parser.error("Expected at least on file name")
    for filename in args:
        for option in options:
            print wrap_extract(filename, option),
        print ""


