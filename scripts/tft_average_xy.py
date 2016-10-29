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

Take some four-dimensional dataset *(X, Y, Z, VALUE)*.  
Divide the dataset into bins along the *Z* axis and calculate the average
*<VALUE>* for each bin.
Return *(Z, <VALUE>)*, where *Z* is the lower bound of the bin.

Command-line usage
---------------------

.. code-block:: bash

    tft_average_xy.py  DATASET  BINS

where *BINS* is the number of bins into which the data is partitioned along the
*Z*-axis.
"""

def average_xy(filename, bins):
    inFile = open(filename)
    x, y, z, value = [], [], [], []
    avgValue = [.0]*(bins+1)
    norm    = [0]*(bins+1)
    zMax=-1e50
    zMin=+1e50
    for line in inFile:
        words=line.split()
        if (words[3] != "nan"): 
            x.append(float(words[0]))
            y.append(float(words[1]))
            z.append(float(words[2]))
            value.append(float(words[3]))
            if (z[-1] > zMax):
                zMax = z[-1]
            if (z[-1] < zMin):
                zMin = z[-1]

    deltaZ = (zMax-zMin)/(bins-1)
    for i in range(len(z)):
        bin = int(z[i]/deltaZ)
        avgValue[int(z[i]/deltaZ)] += value[i]
        norm[int(z[i]/deltaZ)] += 1

    for i in range(bins):
        if avgValue[i] == "nan":
            print i*deltaZ, "nan"
        elif norm[i] == 0:
            print i*deltaZ, "0" 
        else:
            print i*deltaZ, avgValue[i]/norm[i]
    inFile.close()


if __name__ == '__main__':
    import sys
    average_xy(sys.argv[1], int(sys.argv[2]))

