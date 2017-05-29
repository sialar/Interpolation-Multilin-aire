from mpl_toolkits.mplot3d import Axes3D
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

def plot_interpolation_progressively(dt):
    for i in range(nb_functions):
        plt.hold(True)
        print("\t\t\t\t\t\t\t\t\t\t\tx(", i, ") = (", nuX[i], ",", nuY[i], ")")
        print("\t\t\t\t\t\t\t\t\t\t\talpha(", i,") = ", alpha[i], end="\n\n")
        plt.clf()

        ax = fig.gca(projection='3d')
        v = []
        for j in range(nb_points):
            v.append(float(line_interpolation_progression[j].split(" ")[i]))
        ax.plot_trisurf(x1, x2, g, linewidth=0.02, antialiased=True, color='g')
        ax.plot_trisurf(x1, x2, v, linewidth=0.02, antialiased=True, color='r')
        ax.set_xlim(-1.1, 1.1)
        ax.set_ylim(-1.1, 1.1)
        ax.set_zlim(-0.1, 1.1)
        plt.pause(dt)
        separator(206)

def prepare_interpolation_data(lines):
    x1, x2, y = [], [], []
    g, zero = [], []
    alpha, nuX, nuY  = [], [], []
    for k in range(1,nb_points+1):
        x1.append(float(lines[k].split(" ")[0]))
        x2.append(float(lines[k].split(" ")[1]))
        g.append(float(lines[k].split(" ")[2*(nb_functions+1)]))
        zero.append(0.0)
    for i in range(2,nb_functions+2):
        for j in range(1,nb_points+1):
            y.append(float(lines[j].split(" ")[i]))
    for k in range(0,nb_functions):
        alpha.append(float(lines[nb_points+1].split(" ")[k]))
        nuX.append(float(lines[nb_points+2].split(" ")[k]))
        nuY.append(float(lines[nb_points+3].split(" ")[k]))
    return x1, x2, y, g, zero, alpha, nuX, nuY

def prepare_interpolation_progression(lines):
    z = []
    for i in range(nb_functions):
        for j in range(nb_points):
            z.append(float(lines[j].split(" ")[i]))
    return z

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
fig = plt.figure(figsize=taille)


x1, x2, y, g, zero, alpha, nuX, nuY = prepare_interpolation_data(lines_basis_functions)
z = prepare_interpolation_progression(line_interpolation_progression)
plot_interpolation_progressively(dt)
plt.show()
