import matplotlib.colors as colors
import matplotlib.pyplot as plt
import matplotlib.cm as cmx
import numpy as np
import string
import math
import sys

taille = (10,10)
fig = plt.figure(figsize=taille)

variables = ["burnup", "fuel_temperature", "moderator_density", "boron_concentration", "xenon_level"]

if len(sys.argv)!=4:
    print "Choose core type and 2 directions from [burnup, fuel_temperature, moderator_density, boron_concentration, xenon_level]"
    sys.exit(1)

core = sys.argv[1]
direction_x = int(sys.argv[2])
direction_y = int(sys.argv[3])

figureTitle = "Interpolation points projected in directions " + str(direction_x) + " and " + str(direction_y)
fileName = "../data/interpolation_points.dat"
figureName = "../images/" + core + "/Interpolation_Points/interpolation_points_" + str(direction_x) + "_" + str(direction_y)

points_file = open(fileName,"r")
lines = points_file.readlines()
points_file.close()

npts = int(lines[0])

x, y = [], []
for i in range(npts):
    x.append(float(lines[i+1].split(" ")[direction_x]))
    y.append(float(lines[i+1].split(" ")[direction_y]))

plt.scatter(x,y,color='k',s=80,marker='*',alpha=1)

plt.xlabel(variables[direction_x])
plt.ylabel(variables[direction_y])
plt.title(figureTitle)
plt.savefig(figureName)
