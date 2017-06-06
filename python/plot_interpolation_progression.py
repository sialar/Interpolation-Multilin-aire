import matplotlib.colors as colors
import matplotlib.pyplot as plt
import matplotlib.cm as cmx
import numpy as np
import warnings
import random
import sys

def get_cmap(N):
    '''Returns a function that maps each index in 0, 1, ... N-1 to a distinct
    RGB color.'''
    color_norm  = colors.Normalize(vmin=0, vmax=N-1)
    scalar_map = cmx.ScalarMappable(norm=color_norm, cmap='hsv')
    def map_index_to_rgb_color(index):
        return scalar_map.to_rgba(index)
    return map_index_to_rgb_color

def plot_interpolation_progressively(x,g,zero,y,alpha,interp_points,z,dt):
    for i in range(int(len(y)/len(x))):
        plt.hold(True)
        plt.clf()

        plt.subplot(212)
        plt.axis([-1.1, 1.1, y_min, 1.1])
        for j in range(i):
            plt.plot(x, y[j*len(x):(j+1)*len(x)], c='y', alpha=0.7*j/i)
            plt.plot(x, y[i*len(x):(i+1)*len(x)], c='r')

        plt.subplot(211)
        plt.axis([-1.1, 1.1, min(z)-0.1, max(z)+0.1])
        plt.plot(x, zero, 'k', c='k')
        plt.plot(x, g, 'k', c='g')
        plt.plot(x, z[i*len(x):(i+1)*len(x)], c='r')
        msg = ("alpha("+ str(interp_points[i]) + ") = " + str(alpha[i]))
        plt.title(msg)

        plt.pause(dt)

def prepare_interpolation_data(lines):
    x, g, zero, y, alpha, interp_points = [], [], [], [], [], []
    for k in range(1,nb_points+1):
        x.append(float(lines[k].split(" ")[0]))
        g.append(float(lines[k].split(" ")[nb_functions+1]))
        zero.append(0.0)
    for i in range(1,nb_functions+1):
        for j in range(1,nb_points+1):
            y.append(float(lines[j].split(" ")[i]))
    for k in range(0,nb_functions):
        alpha.append(float(lines[nb_points+1].split(" ")[k]))
    for k in range(0,nb_functions):
        interp_points.append(float(lines[nb_points+2].split(" ")[k]))
    return x, g, zero, y, alpha, interp_points

def prepare_interpolation_progression(lines):
    z = []
    for i in range(nb_functions):
        for j in range(nb_points):
            z.append(float(lines[j].split(" ")[i]))
    return z

def get_y_min():
    if method==0:
        return -1.1
    elif method==1:
        return -0.1
    else:
        return -0.3

input_file = open( "../data/basis_functions.txt", "r")
lines_basis_functions = input_file.readlines()
input_file.close()

input_file = open( "../data/interpolation_progression.txt", "r")
line_interpolation_progression = input_file.readlines()
input_file.close()

nb_functions = int(lines_basis_functions[0].split(" ")[0])
nb_points = int(lines_basis_functions[0].split(" ")[1])
method = int(lines_basis_functions[0].split(" ")[2])

dt = 0.5
nb_colors = 10
cmap = get_cmap(nb_colors)
taille = (25,10)
plt.figure(figsize=taille)

y_min = get_y_min()
x, g, zero, y, alpha, interp_points = prepare_interpolation_data(lines_basis_functions)
z = prepare_interpolation_progression(line_interpolation_progression)
plot_interpolation_progressively(x,g,zero,y,alpha,interp_points,z,dt)
plt.show()
