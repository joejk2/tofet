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
:mod:`tft_make_cubic_lattice.py`
=================================
.. moduleauthor:: Joe Kwiatkowski

Make a cubic lattice, outputting :term:`XYZ <XYZ file>` and :term:`edge <Edge file>`
files.  

Command-line usage
------------------

.. code-block:: bash
    
    tft_make_cubic_lattice.py  XMAX  YMAX  ZMAX  COL_GEN  J  DE

where *(XMAX, YMAX, ZMAX)* are the maximum dimensions of the sample,
:term:`COL_GEN`
is the fraction of the film which is used for collection and generation,
:term:`J` is
the :term:`transfer integral` between :term:`vertices <vertex>`, and :term:`DE` is the difference in energies
of :term:`vertices <vertex>`.

Currently, only supports the format in which energies are placed in the
:term:`edge file <Edge file>`
"""


#DO WORK
########
def print_xyz(xyz_filename, dimensions, coll_gen_thickness):
    out_file = open(xyz_filename, "w")
    xLength, yLength, zLength = dimensions
    for z in range(zLength):
        #What type of vertex?
        if (z < coll_gen_thickness): #Need > here since z starts at 0
            type='g' 
        elif (z >= zLength - coll_gen_thickness): #Need >= here since z ends at zLength-1
            type='c' 
        else: type='-'
        for y in range(yLength):
            for x in range(xLength):
                out_file.write( "%d\t%d\t%d\t%s\n" % ( x, y, z, type ) )

def print_edges(edge_filename, dimensions, J, DE):
    """Print to out_filename
    - each vertex has neighbours at +x, +y, +z, unless at edges
    - x neighbour at vertexID + 1
    - no neighbour whenever vertexID % xLength = xLength - 1
    - y neighbour at vertexID + xLength
    - no neighbour whenever vertexID % (xLength*yLength) > xLength*yLength - xLength
    - z neighbour at vertexID + xLength*yLength
    - no neighbour whenever vertexID > xLength*yLength*(zLength - 1)"""
    out_file = open(edge_filename, "w")
    xLength, yLength, zLength = dimensions
    edgeCount = 0   #Count number of edges
    for i in range(xLength * yLength * zLength):
        if ( i%xLength != xLength-1 ):
            out_file.write( "%d\t%d\t%f\t%f\n" % (i, i+1, J, DE) )
            edgeCount += 1
        if ( i%(xLength*yLength) < xLength*yLength-xLength ):
            out_file.write( "%d\t%d\t%f\t%f\n" % (i, i+xLength, J, DE) )
            edgeCount += 1
        if ( i < xLength*yLength*(zLength-1) ): 
            out_file.write( "%d\t%d\t%f\t%f\n" % (i, i+xLength*yLength, J, DE) )
            edgeCount += 1

    #Test we've got the expected number of edges
    expectedEdges = 3 * xLength * yLength * zLength \
                  - xLength * yLength \
                  - xLength * zLength \
                  - yLength * zLength 
    try:
        assert edgeCount == expectedEdges
    except AssertionError:
        print "*** ERROR ***: You've got an incorrect number of edges..."
        print "               Number of edges = %d" % edgeCount
        print "               Expected number = %d" % expectedEdges 


if __name__ == "__main__":
    import sys
    xLength = int(sys.argv[1])
    yLength = int(sys.argv[2])
    zLength = int(sys.argv[3])
    col_gen = float(sys.argv[4])
    J = float(sys.argv[5])
    DE = float(sys.argv[6])
    print_xyz("output.xyz", (xLength, yLength, zLength), col_gen)
    print_edges("output.edge", (xLength, yLength, zLength), J, DE)
