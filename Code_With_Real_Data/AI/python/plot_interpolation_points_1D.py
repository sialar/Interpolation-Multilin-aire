import matplotlib.colors as colors
import matplotlib.pyplot as plt
import matplotlib.cm as cmx
import numpy as np
import string
import math
import sys

fig = plt.figure()

variables = ["burnup", "fuel_temperature", "moderator_density", "boron_concentration", "xenon_level"]

if len(sys.argv)!=3:
    print "Choose core type and one direction from [burnup, fuel_temperature, moderator_density, boron_concentration, xenon_level]"
    sys.exit(1)

core = sys.argv[1]
direction = int(sys.argv[2])

figureTitle = "Interpolation points projected in direction " + str(direction)
fileName = "../data/interpolation_points.dat"
figureName = "../images/" + core + "/Interpolation_Points/interpolation_points_" + str(direction)

points_file = open(fileName,"r")
lines = points_file.readlines()
points_file.close()

npts = int(lines[0])

x, zero = [], []
for i in range(npts):
    x.append(float(lines[i+1].split(" ")[direction]))
    zero.append(0.0)

plt.scatter(x,zero,color='r',s=10,marker='o',alpha=1)

axes = plt.gca()
xmin = 1.2*min(x)
if xmin == 0:
    xmin = -0.2*max(x)
axes.set_xlim([xmin,1.2*max(x)])

#print set(x)

plt.xlabel(variables[direction])
plt.title(figureTitle)
plt.savefig(figureName)
