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
:mod:`tft_rotate_coords.py`
=============================

.. moduleauthor:: Joe Kwiatkowski

Rotates :term:`xyz file` to align Z-axis perpendicular to a unit cell face. 

Command-line usage
---------------------

.. code-block:: bash

    tft_rotate_coords.py   XYZ  CELL  CROSS_PRODUCT

where *XYZ* is the :term:`XYZ file` (that must end with the 
suffix *.xyz*), *CELL* is a 
file containing the dimensions of the morphology (and must end with the 
suffix *.cell*).

The format of the *CELL* file is::

    a       b       c  
    alpha   beta    gamma

where the first line relates to the length of your primitive basis vectors
(arbitrary units), and the second line the angles between them (in degrees).
For a rectangular lattice this would simply be::

    a       b       c  
    90      90      90

The input cross_product is the vector cross product that defines the
direction you want the field to point along. 

.. note:: 
    It is only physical to apply the electric field perpendicular to 
    the faces of your cell. 

The allowed values of cross_product 
    are *aXb*, *aXc*, *bXa*, *bXc*, *cXa*, *cXb*.
"""

import sys
from numpy import *
from math import radians, atan2, acos, atan


def rotate_coords(cell_filename, xyz_filename, cross_product):
    inCell=open(cell_filename, "r")
    lines=inCell.readlines()
    lengths=lines[0].split()
    angles =lines[1].split()
    a_mag,b_mag,c_mag = double(lengths[0]),double(lengths[1]),double(lengths[2])
    alpha,beta,gamma  = radians(double(angles[0])),radians(double(angles[1])),radians(double(angles[2]))

    #Define ScaleInv (frac -> cart)
    volume=a_mag*b_mag*c_mag*sin(gamma)*sqrt(1-cos(beta)*cos(beta)-cos(alpha)*cos(alpha));
    scaleInv=array([ [a_mag, b_mag*cos(gamma), c_mag*cos(beta)  ],
                   [0,     b_mag*sin(gamma), c_mag*(cos(alpha)-cos(beta)*cos(gamma))/sin(gamma)],
                   [0,     0,                volume/(a_mag*b_mag*sin(gamma))] ])

    #Determine E
    a=dot(scaleInv,[1,0,0])
    b=dot(scaleInv,[0,1,0])
    c=dot(scaleInv,[0,0,1])
    if   "aXb" in cross_product:
        print "E already along z"
        sys.exit()
    elif "aXc" in sys.argv[1:]:
        E=cross(a,c)
    elif "bXa" in sys.argv[1:]:
        E=cross(b,a)
    elif "bXc" in sys.argv[1:]:
        E=cross(b,c)
    elif "cXa" in sys.argv[1:]:
        E=cross(c,a)
    elif "cXb" in sys.argv[1:]:
        E=cross(c,b)
    else:
        print "What direction?"

    # Rotate z onto E
    theta=atan2(E[1],E[0])
    phi  =acos(E[2]/sqrt(dot(E,E)))
    R=array( [ [ cos(theta)*cos(phi), sin(theta)*cos(phi), -sin(phi) ],
               [-sin(theta)         , cos(theta)         ,  0        ],
           [ cos(theta)*sin(phi), sin(theta)*sin(phi),  cos(phi) ] ] )

    # Read XYZ file, and rotate.
    inXYZ=open(xyz_filename,"r")
    coords=array([])
    for line in inXYZ:
        words=line.split()
        coords=append(coords,[(float(words[0]),float(words[1]),float(words[2]))])
        coords=coords.reshape(size(coords)/3,3)

    for i in range (0,len(coords)):
        coords[i]=dot(R,coords[i])

    for i in range(0,len(coords)):
        print "%f\t%f\t%f" % (coords[i,0],coords[i,1],coords[i,2])



if __name__ == '__main__':
    if (len(sys.argv)<4):
        print "Give:"
        print "  1) ***.xyz file"
        print "  2) ***.cell file"
        print "  3) Two primitive vectors whose cross product defines new z"
        sys.exit(-1)

    # Read definition of unit cell
    for arg in sys.argv[1:]:
        if "cell" in arg:
            cell_filename = arg
            break
    else:
        print "Please give file ***.cell"

    # Read in existing coords
    for arg in sys.argv[1:]:
        if "xyz" in arg:
            xyz_filename = arg
            break
    else: 
        print "Need xyz file!"
        sys.exit()

    # Read cross product
    cross_products = ["aXb", "aXc", "bXa", "bXc", "cXa", "cXb"]
    for arg in sys.argv[1:]:
        if arg in cross_products:
            cross_product = arg
            break
    else:
        print "Need a cross product in", cross_products 

    rotate_coords(cell_filename, xyz_filename, cross_product)
