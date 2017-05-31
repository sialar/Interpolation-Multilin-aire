from mpl_toolkits.mplot3d import Axes3D # librairie 3D
import matplotlib.pyplot as plt
import numpy as np
import sys

def prepare_initial_grid(methods,n,m):
    leja_sequence_file = open( "../data/leja_sequence.txt", "r")
    leja_lines = leja_sequence_file.readlines()
    leja_sequence_file.close()

    dichotomy_file = open( "../data/dichotomy_sequence.txt", "r")
    dichotomy_lines = dichotomy_file.readlines()
    dichotomy_file.close()

    leja_sequence_x, leja_sequence_y = [], []
    dichotomy_x, dichotomy_y = [], []
    leja_values, dichotomy_values = [], []
    x, y = [], []
    nb_pts_leja = min(n,m)
    nb_pts_dichotomy = int(dichotomy_lines[0])

    for i in range(nb_pts_leja):
        leja_values.append(float(leja_lines[i]))
    for i in range(nb_pts_dichotomy):
        dichotomy_values.append(float(dichotomy_lines[i+1]))

    if methods[0]==0 and methods[1]==0:
        print(1)
        for v1 in leja_values:
            for v2 in leja_values:
                    leja_sequence_x.append(v1)
                    leja_sequence_y.append(v2)
        x = leja_sequence_x
        y = leja_sequence_y

    elif methods[0]==0 and methods[1]!=0:
        print(2)
        for v1 in leja_values:
            for v2 in dichotomy_values:
                    leja_sequence_x.append(v1)
                    dichotomy_y.append(v2)
        x = leja_sequence_x
        y = dichotomy_y

    elif methods[0]!=0 and methods[1]==0:
        print(3)
        for v1 in dichotomy_values:
            for v2 in leja_values:
                    dichotomy_x.append(v1)
                    leja_sequence_y.append(v2)
        x = dichotomy_x
        y = leja_sequence_y

    elif methods[0]!=0 and methods[1]!=0:
        print(4)
        for v1 in dichotomy_values:
            for v2 in dichotomy_values:
                    dichotomy_x.append(v1)
                    dichotomy_y.append(v2)
        x = dichotomy_x
        y = dichotomy_y

    plt.axis([-1.1, 1.1, -1.1, 1.1])
    plt.scatter(x, y, s=s, c='w')

def plot_picked_points_progressively_2d(indices,points,s,dt):
    plt.hold(True)
    plt.axis([-1.1, 1.1, -1.1, 1.1])
    for i in range(len(points)):
        p = points[i]
        plt.scatter(p[0], p[1], s=s, c='r')
        plt.pause(dt)

def plot_picked_points_progressively_3d(indices,points,s,dt):
    ax = fig.gca(projection='3d')
    ax.set_xlim(-1.1, 1.1)
    ax.set_ylim(-1.1, 1.1)
    ax.set_zlim(-1.1, 1.1)
    plt.hold(True)
    for i in range(len(points)):
        p = points[i]
        ax.scatter(p[0], p[1], p[2], c='r', s=s, marker='o')
        plt.pause(dt)


all_points, picked_points, picked_indices = [], [], []

if len(sys.argv)!=2:
    print("Invalid number of arguments choose the space dimension.")
    sys.exit(0)
else:
    dim = int(sys.argv[1])

taille = (16,8)
fig = plt.figure(figsize=taille)

ai_output_file = open( "../data/path.txt", "r")
lines = ai_output_file.readlines()
ai_output_file.close()

if int(lines[0]):
    l = 1
    n = int(lines[1].split(" ")[0])
    m = int(lines[1].split(" ")[1])
    if dim==3:
        l = int(lines[1].split(" ")[2])

    offset = 4
    s = 80
    dt = 0.01

    methods = []
    for i in range(dim):
        methods.append(int(lines[2].split(" ")[i]))

    nb_points = int(lines[3])

    for i in range(nb_points):
        x = float(lines[i+offset].split(" ")[0])
        y = float(lines[i+offset].split(" ")[1])
        z = 0
        if dim==3:
            z = float(lines[i+offset].split(" ")[2])
        picked_points.append((x,y,z))


    if dim==2:
        prepare_initial_grid(methods,n,m)
        plot_picked_points_progressively_2d(picked_indices, picked_points,s,dt)
    else:
        plot_picked_points_progressively_3d(picked_indices, picked_points,s,dt)
    plt.show()
