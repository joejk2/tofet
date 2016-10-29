#!/usr/bin/python 
import sys
import pylab
from scipy.stats import linregress

traps = [20100]

inFile = open(sys.argv[1])
temp = float(sys.argv[2])
dielectric = float(sys.argv[3])
kT = temp * 8.617343e-5

x, y, z, prob_read, prob_calc = [], [], [], [], []
for line in inFile:
    words = line.split()
    x.append(float(words[0]))
    y.append(float(words[1]))
    z.append(float(words[2]))
    prob_read.append(float(words[3]))

for i in range(len(x)):
    E = .0
    for trap in traps:
        if trap != i:
            r_sqd = (x[i] - x[trap])**2\
                  + (y[i] - y[trap])**2\
                  + (z[i] - z[trap])**2
            E += 14.3996442 / (pylab.sqrt(r_sqd) * dielectric)
            prob_calc.append(pylab.exp(-E / kT))
        else: 
            prob_calc.append(0.)
            
prob_calc /=  pylab.sqrt(pylab.dot(prob_calc, prob_calc))
prob_read /= pylab.sqrt(pylab.dot(prob_read, prob_read))

print linregress(prob_read, prob_calc)
pylab.plot(prob_read, prob_calc)
pylab.plot([0, 0.008], [0, 0.008])
pylab.show()
