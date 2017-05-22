import matplotlib.colors as colors
import matplotlib.pyplot as plt
import matplotlib.cm as cmx
import numpy as np
import warnings
import random
import sys

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

def plot_basis_function_progressively(alpha,interp_points,x,y,dt,nb_colors):
    for i in range(int(len(y)/len(x))):
        plt.hold(True)
        print("\t\t\t\t\t\t\t\t\t\t\tx(", i, ") = ", interp_points[i])
        print("\t\t\t\t\t\t\t\t\t\t\talpha(", i,") = ", alpha[i], end="\n\n")
        for j in range(i):
            plt.plot(x, y[j*len(x):(j+1)*len(x)], c='g')

        plt.plot(x, y[i*len(x):(i+1)*len(x)], c='r')
        plt.pause(dt)
        separator(206)

input_file = open( "plot_function.txt", "r")
lines = input_file.readlines()
input_file.close()

nb_functions = int(lines[0].split(" ")[0])
nb_points = int(lines[0].split(" ")[1])
x, y, res, real_res = [], [], [], []
interp_points, alpha, zeros = [], [], []

taille = (25,10)
plt.figure(figsize=taille)

nb_colors = 10
cmap = get_cmap(nb_colors)

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
y_min = - (0.1 + (1 - int(sys.argv[1])))
plt.axis([-1.1, 1.1, y_min, 1.1])
plt.plot(x, zeros, 'k', c='k')

for k in range(0,nb_functions):
    alpha.append(float(lines[nb_points+1].split(" ")[k]))

for k in range(0,nb_functions):
    interp_points.append(float(lines[nb_points+2].split(" ")[k]))

plot_basis_function_progressively(alpha,interp_points,x,y,1.5,nb_colors)
plt.show()
separator(206)
