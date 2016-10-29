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

Calculate the Electrochemical potential for each molecule.

Command-line usage
--------------------

.. code-block:: bash 

    tft_calc_ecp.py  OUTPUT_FILE 

where *OUTPUT_FILE* contains the energies and occupation probabilities of each
molecule.

"""


import sys
from math import log
import tft_extract

k = 8.617343e-5  # eV / K

def calcFE(filename):
    temp = float(tft_extract.wrap_extract(filename, "--temp"))
    energies = tft_extract.wrap_extract(filename, "--energies").split('\n')
    energies.remove('')
    occ_probs = tft_extract.wrap_extract(filename, "--occ_prob").split('\n')
    occ_probs.remove('')
    if (len(energies) != len(occ_probs)):
        print "*** ERROR ***: Incorrect read of energies " + \
              "and / or occupation probabilities"
        return 1
    x, y, z, E, prob, fermi_energy = [], [], [], [], [], []
    for row in energies:
        words = row.split()
        x.append(words[0])
        y.append(words[1])
        z.append(words[2])
        E.append(float(words[3]))
    for row in occ_probs:
        words = row.split()
        prob.append(float(words[3]))
    for i in range(len(occ_probs)):
        if (prob[i] == 0.0) or (prob[i] == 1.0):  # no singularities thank you!
            fermi_energy.append("nan")
        else:
            fermi_energy.append(E[i] - k * temp * log(1.0 / prob[i] - 1.0))
    return (x, y, z, fermi_energy)
        

if __name__ == '__main__':
    data = calcFE(sys.argv[1])
    
    if not isinstance(data, int):
        x, y, z, fermi_energy = data
        for i in range(len(x)):
            print x[i], y[i], z[i], fermi_energy[i]

