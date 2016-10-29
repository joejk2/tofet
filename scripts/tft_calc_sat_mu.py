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

Calculate field-effect mobilities.

Command-line usage
--------------------

.. code-block:: bash
    
    tft_calc_sat_mu.py  FILES  -l LENGTH  -w WIDTH  [OPTIONS]

where *FILES* are ToFeT output files specified with a :term:`glob pattern`. 
*LENGTH* and *WIDTH* are the length and width of the simulated FET
channel, measured in Angstrom. 

Options
--------

.. program:: tft_calc_mu_sat.py

.. cmdoption:: debug

    Output information useful for debugging.
"""


import sys
from tft_get_values import getValuesFromMultipleFiles
from math import sqrt
from scipy import stats
import numpy

def calc_mu_sat(filenames, W, L, debug=False):
    Vg_current = getValuesFromMultipleFiles(filenames + ["-o", "Vg", "fet_current"], False)
    Vg_current = numpy.array(Vg_current, float)
    Vg = Vg_current[:,0]
    sqrt_current = numpy.sqrt(Vg_current[:,1])
    #sqrt_current *= sqrt(1.60217646e-19)  # hoppers -> Coulombs
    Vg_sqrt_current_gradient = stats.linregress(Vg, sqrt_current)[0]

    hoppers =  getValuesFromMultipleFiles(filenames + ["-o", "hoppers"], False)
    hoppers = numpy.array(hoppers, float)[:,0]
    hoppers *= 1.60217646e-19  # hoppers -> Coulombs
    Vg_hoppers_gradient = stats.linregress(Vg, hoppers)[0]

    C = Vg_hoppers_gradient / (L * W) 

    if debug:
        print "Vg_current:"
        print Vg_current
        print "hoppers:"
        print hoppers
        print "L = ", L, "Ang"
        print "W = ", W, "Ang"
        print "Vg_sqrt_current_gradient = ", Vg_sqrt_current_gradient, "A^0.5/V"
        print "Vg_hoppers_gradient = ", Vg_hoppers_gradient, "C/V"
        print "C = ", C, "C/(V.Ang^2)"

    return  2*L/(W*C) * Vg_sqrt_current_gradient**2 * 1e-16

if __name__ == '__main__':
    try:
        W_index = sys.argv.index("-w")
        L_index = sys.argv.index("-l")
    except ValueError:
        print "*** ERROR *** : Need both -l and -w options in command line!"
        exit()
    try: 
        sys.argv.index("debug")
    except ValueError:
        debug = False
    else:
        debug = True
        sys.argv.remove("debug")
    W = sys.argv[W_index + 1]
    L = sys.argv[L_index + 1]
    sys.argv.remove("-w")
    sys.argv.remove("-l")
    sys.argv.remove(W)
    sys.argv.remove(L)
    mu = calc_mu_sat(sys.argv[1:], float(W), float(L), debug)
    print mu, " cm^2 / V.s"
