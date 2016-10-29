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
:mod:`tft_run_batch.py`
==========================

.. moduleauthor:: Joe Kwiatkowski

Run ToFeT in batch.

Command-line usage
--------------------

.. code-block:: bash
    
    tft_run_batch.py  XYZ  EDGE  SIM

where, as is normal for a simple :mod:`tft` simulation, *XYZ* is the 
:term:`XYZ file`, *EDGE* is the :term:`Edge file`, and *SIM* is the 
:term:`Sim file`.

After each variable in *SIM*, list all the values you would like the variable to
take, separated by blank space::

    temp 200 250 300 350  

:mod:`tft_run_batch.py` then scans through foo.sim and runs tofet for each value of the
field, putting the output in the files 0.out, 1.out, ....

If you so wish, multiple variables can have multiple values. For example::

    temp 200 250 300 350  
    fieldZ -4e-4 -2e-4  
    reorg 0.132 0.142  

In this more complicated case you might have to look into each output
file to work out which parameters were used in the simulation. 

:mod:`tft_run_batch.py` is less forgiving to comments than the unwrapped tofet. The only
robust way to include comments when using :mod:`tft_run_batch.py` is comment out entire lines
in the :term:`sim file` file by including a # anywhere. Even if # isn`t at the beginning
of the line it comments out the entire line.

If you want to call :mod:`tft_occ` rather than :mod:`tft`, add occ to the argument
list of :mod:`tft_run_batch.py`:

.. code-block:: bash 

    tft_run_batch.py  SIM  XYZ  EDGE  occ

"""
        

import sys
import os
from string import *

binary, sim, xyz, edge, occ = "tft", "", "", "", ""

def setup():
    global binary, sim, xyz, edge, occ
    #Open master sim file
    for arg in sys.argv[1:]:
        if   ".sim" in arg: sim = arg
        elif ".xyz" in arg: xyz = arg
        elif ".edge" in arg: edge = arg
        elif ".occ" in arg: occ = arg

        elif (arg=="occ"):   binary = "tftOccupation"

    #Check that you've got all the input files:
    try: 
        assert ( os.path.isfile(sim) and os.path.isfile(xyz) and os.path.isfile(edge) ) == True
    except:
        print "*** ERROR ***: Can't find one of the input files ***.sim, ***.xyz, ***.edge"
        print sim, xyz, edge
        sys.exit()

    #Read reorg, temp, field etc...
    simIn=open(sim,"r")
    input=[[]]
    for line in simIn:
        if "#" not in line:
            input.append(line.split())
    input = [x for x in input if x]
    return input

def cross(args):
    counter=-1
    ans = [[]]
    for arg in args:
        ans = [x+[y] for x in ans for y in arg[1:]]
    for perm in ans:
        counter+=1
        simOut=open("tmp.sim","w")
        for var in range(len(perm)):
            print >> simOut, args[var][0], perm[var]
        simOut.close()
        # Run ToFeT:
        out = "%d.out" % counter 
        cmd = "%s tmp.sim %s %s %s > %s" % (binary, xyz, edge, occ, out)
        print "ToFeT.py: running simulation %d of %d with %s" % (counter+1,len(ans),binary)
        os.system(cmd)
        # If there was a '.occ' file, then update it
        if occ:
            cmd = "tft_extract.py --occ_mol " + out + " > " + occ  
            os.system(cmd)
    return counter+1

def main():
    input = setup()
    counter = cross(input)
    #Clean up
    os.system("rm -f tmp.sim")

if __name__ == '__main__':
    main()
