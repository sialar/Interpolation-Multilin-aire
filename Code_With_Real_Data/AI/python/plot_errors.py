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
csName, iters = [], []

figureTitle = "Convergence speed of cross-sections " + core
fileName = "../data/" + core + "/RelativeErrors"
figureName = "../images/" + core + "/RelativeErrors"

results_file = open(fileName,"r")
lines = results_file.readlines()
results_file.close()

niters = int(lines[0])
for i in range(niters):
    iters.append(int(lines[1].split(" ")[i]))
nbcs = int(lines[2])
for i in range(nbcs):
    csName.append(lines[3].split(" ")[i])

cmap = get_cmap(nbcs)
cs_res = {}

nbh = int(math.sqrt(nbcs))
nbv = nbcs/nbh


for i in range(nbcs):
    cs_res[csName[i]] = []
    for j in range(4,len(lines)):
        cs_res[csName[i]].append(abs(float(lines[j].split(" ")[i])))

for i in range(nbcs):
    ax = fig.add_subplot(nbh,nbv,i+1)
    ax.set_title(csName[i])
    plt.plot(iters, cs_res[csName[i]], c='k')


plt.show()
