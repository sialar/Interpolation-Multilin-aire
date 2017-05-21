import numpy as np
import matplotlib.pyplot as plt
import sys
from random import randint
import matplotlib.cm as cmx
import matplotlib.colors as colors
import warnings

def separator(N):
    for i in range(N):
        print("*",end='')
    print("\n")

def get_cmap(N):
    '''Returns a function that maps each index in 0, 1, ... N-1 to a distinct
    RGB color.'''
    color_norm  = colors.Normalize(vmin=0, vmax=N-1)
    scalar_map = cmx.ScalarMappable(norm=color_norm, cmap='hsv')
    def map_index_to_rgb_color(index):
        return scalar_map.to_rgba(index)
    return map_index_to_rgb_color

def plot_basis_function_progressively(alpha,interp_points,x,y,dt):
    for i in range(int(len(y)/len(x))):
        plt.hold(True)
        plt.plot(x, y[i*len(x):(i+1)*len(x)], c=cmap(i))
        plt.pause(dt)
        separator(206)
        print("x(", i, ") = ", interp_points[i])
        print("alpha(", i,") = ", alpha[i], end="\n\n")


input_file = open( "plot_function.txt", "r")
lines = input_file.readlines()
input_file.close()

nb_functions = int(lines[0].split(" ")[0])
nb_points = int(lines[0].split(" ")[1])
x, y, res, real_res = [], [], [], []
interp_points, alpha, zeros = [], [], []

cmap = get_cmap(nb_functions)

for k in range(1,nb_points+1):
    x.append(float(lines[k].split(" ")[0]))
    zeros.append(0.0)
for i in range(1,nb_functions+1):
    for j in range(1,nb_points+1):
        y.append(float(lines[j].split(" ")[i]))

plt.subplot(211)
for j in range(1,nb_points+1):
    res.append(float(lines[j].split(" ")[nb_functions+1]))
    real_res.append(float(lines[j].split(" ")[nb_functions+2]))
plt.plot(x, zeros, 'k', c='k')
plt.plot(x, res, 'k', c='g')
plt.plot(x, real_res, 'k', c='r')

plt.subplot(212)
plt.axis([-1.1, 1.1, -1.1, 1.1])
plt.plot(x, zeros, 'k', c='k')

for k in range(0,nb_functions):
    alpha.append(float(lines[nb_points+1].split(" ")[k]))

for k in range(0,nb_functions):
    interp_points.append(float(lines[nb_points+2].split(" ")[k]))

plot_basis_function_progressively(alpha,interp_points,x,y,1)
plt.show()
separator(206)
