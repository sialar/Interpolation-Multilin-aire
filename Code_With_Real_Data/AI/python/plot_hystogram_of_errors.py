
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
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

def get_border(errors):
    if len(errors)==0:
        return 0
    return max(max(errors),abs(min(errors)))


fig = plt.figure()

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

figureTitle = "Hystogram of relative errors for the cross-section " + core + "_" + csName
figuresDir = "../images/" + core

ai_errors, tucker_errors, cocagne_errors = [], [], []

for i in range(n):
    cocagne_errors.append(float(lines[i].split(" ")[10]))
    tucker_errors.append(float(lines[i].split(" ")[11]))
    ai_errors.append(float(lines[i].split(" ")[12]))

errors_border = max(get_border(cocagne_errors),get_border(tucker_errors),get_border(ai_errors))
bins = int(errors_border/10)
plt.hist(cocagne_errors, facecolor='b', histtype = 'step', alpha = 1, normed = True)
plt.hist(tucker_errors, facecolor='g', histtype = 'step', alpha = 1, normed = True)
plt.hist(ai_errors, facecolor='r', histtype = 'step', alpha = 1, normed = True)

plt.show()
