import matplotlib.colors as colors
import matplotlib.pyplot as plt
import matplotlib.cm as cmx
import numpy as np
import string
import math
import sys

fig = plt.figure()


if len(sys.argv)!=2:
    print("Invalid number of arguments, choose one direction")
    sys.exit(1)

direction = int(sys.argv[1])

points_file = open("../data/interpolation_points.dat","r")
lines = points_file.readlines()
points_file.close()

figureTitle = "Interpolation points in direction " + str(direction)

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

plt.xlabel(str(direction))
plt.title(figureTitle)
plt.show()
