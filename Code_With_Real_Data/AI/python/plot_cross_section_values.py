import matplotlib.colors as colors
import matplotlib.pyplot as plt
import matplotlib.cm as cmx
import numpy as np
import string
import math
import sys

def get_cmap(N):
    '''Returns a function that maps each index in 0, 1, ... N-1 to a distinct
    RGB color.'''
    color_norm  = colors.Normalize(vmin=0, vmax=N-1)
    scalar_map = cmx.ScalarMappable(norm=color_norm, cmap='hsv')
    def map_index_to_rgb_color(index):
        return scalar_map.to_rgba(index)
    return map_index_to_rgb_color

taille = (20,20)
fig = plt.figure(figsize=taille)

core = sys.argv[1]
csName = sys.argv[2]

n = 0
if core == "MOX":
    n = 14000
elif core == "UOX":
    n = 14700
else:
    n = 13300

cmap = get_cmap(5)

fileName = "../data/" + core + "/" + csName
results_file = open(fileName,"r")
lines = results_file.readlines()
results_file.close()

parameters = [ "burnup", "fuel temperature", "moderator density", "boron concentration", "xenon level"]
figureTitle = "Cross-section " + " as function of "
figuresDir = "../images/" + core

all_points, all_res = [], []

for axis in range(5):
    values = []
    points = []
    meanvalues = []
    meanpoints = []
    for i in range(n):
        points.append(float(lines[i].split(" ")[axis+1]))
        values.append(float(lines[i].split(" ")[6]))

    mean = 0.0
    for i in range(1,len(points)):
        if points[i]==points[i-1]:
            mean = mean + values[i-1]
        else:
            meanpoints.append(points[i-1])
            meanvalues.append(mean/i)

    #print points, len(points)
    #print meanpoints, len(meanpoints)

    #print values, len(values)
    #print meanvalues, len(meanvalues)

    all_points.append(meanpoints)
    all_res.append(meanvalues)

nbh = int(math.sqrt(5))+1
nbv = int(5/nbh) + 1

#for i in range(5):
#    print nbh,nbv,i
#    ax = fig.add_subplot(nbh,nbv,i+1)
#    ax.set_title(figureTitle + parameters[i])
#    ax.set_xlabel(parameters[i])
#    ax.set_ylabel(csName)


plt.plot(all_points[0], all_res[0], c='k')
plt.show()
