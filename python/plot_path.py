from mpl_toolkits.mplot3d import Axes3D # librairie 3D
import matplotlib.pyplot as plt
import numpy as np
import sys

def prepare_initial_grid_for_indices(points,n,m):
    x, y = [], []
    for i in points:
        x.append(i[0])
        y.append(i[1])

    plt.subplot(121)
    plt.axis([-1, n, -1, m])
    print(min(n,m))
    if (max(n,m)<100):
        plt.scatter(x, y, s=s/2, c='w')

def prepare_initial_grid_for_leja_points(n,m):
    leja_sequence_file = open( "../data/leja_sequence.txt", "r")
    line = leja_sequence_file.readlines()
    leja_sequence_file.close()
    leja_sequence_x, leja_sequence_y = [], []
    values = []
    nb_pts = 20
    for i in range(nb_pts):
        values.append(float(line[i]))
    for v1 in values:
        for v2 in values:
            leja_sequence_x.append(v1)
            leja_sequence_y.append(v2)
    plt.subplot(122)
    plt.axis([-1.1, 1.1, -1.1, 1.1])
    plt.scatter(leja_sequence_x, leja_sequence_y, s=s/2, c='w')

def prepare_initial_grid_for_middle_points():
    dichotomy_file = open( "../data/dichotomy_sequence.txt", "r")
    line = dichotomy_file.readlines()
    dichotomy_file.close()
    dichotomy_x, dichotomy_y = [], []
    values = []
    nb_pts = int(line[0])
    for i in range(1,nb_pts+1):
        values.append(float(line[i]))
    for v1 in values:
        for v2 in values:
            dichotomy_x.append(v1)
            dichotomy_y.append(v2)
    plt.axis([-1.1, 1.1, -1.1, 1.1])
    plt.scatter(dichotomy_x, dichotomy_y, s=s, c='w')

def plot_picked_points_progressively_2d(indices,points,s,dt):
    plt.hold(True)
    for i in range(len(points)):
        if method==0:
            nu = indices[i]
            plt.subplot(121)
            plt.axis([-1, n, -1, m])
            plt.scatter(nu[0], nu[1], s=s/2, c='r')
            plt.subplot(122)

        p = points[i]
        plt.axis([-1.1, 1.1, -1.1, 1.1])
        if method==0:
            plt.scatter(p[0], p[1], s=s/2, c='r')
        else:
            plt.scatter(p[0], p[1], s=s, c='r')
        plt.pause(dt)

def plot_picked_points_progressively_3d(indices,points,s,dt):
    if method!=0:
        ax = fig.gca(projection='3d')
        ax.set_xlim(-1.1, 1.1)
        ax.set_ylim(-1.1, 1.1)
        ax.set_zlim(-1.1, 1.1)
    else:
        ax1 = fig.add_subplot(1, 2, 1, projection='3d')
        ax1.set_xlim(-1, n)
        ax1.set_ylim(-1, m)
        ax1.set_zlim(-1, l)
        ax2 = fig.add_subplot(1, 2, 2, projection='3d')
        ax2.set_ylim(-1.1, 1.1)
        ax2.set_xlim(-1.1, 1.1)
        ax2.set_zlim(-1.1, 1.1)
    plt.hold(True)
    for i in range(len(points)):
        p = points[i]
        if method==0:
            nu = indices[i]
            ax1.scatter(nu[0], nu[1], nu[2], s=s/2, c='r', marker='o')
            ax2.scatter(p[0], p[1], p[2], c='r', s=s/2, marker='o')
        else:
            p = points[i]
            ax.scatter(p[0], p[1], p[2], c='r', s=s, marker='o')
        plt.pause(dt)


all_points, picked_points, picked_indices = [], [], []

if len(sys.argv)!=3:
    print("Invalid number of arguments: choose the space dimension and the interpolation method.")
    sys.exit(0)
else:
    dim = int(sys.argv[1])
    method = int(sys.argv[2])

taille = (16,8)
fig = plt.figure(figsize=taille)

ai_output_file = open("../data/path.txt", "r")
lines = ai_output_file.readlines()
ai_output_file.close()

n = int(lines[0].split(" ")[0])
m = int(lines[0].split(" ")[1])

if dim==3:
    l = int(lines[0].split(" ")[2])
else:
    l = 1

nb_points = int(lines[1])
offset = 2
s = 80
dt = 0.05

for i in range(n):
    for j in range(m):
        for k in range(l):
            all_points.append((i,j,k))

for i in range(nb_points):
    x = float(lines[i+offset].split(" ")[0])
    y = float(lines[i+offset].split(" ")[1])
    z = 0
    if dim==3:
        z = float(lines[i+offset].split(" ")[2])
    picked_points.append((x,y,z))
    if method==0:
        a = int(lines[i+offset].split(" ")[dim])
        b = int(lines[i+offset].split(" ")[dim+1])
        c = 0
        if dim==3:
            c = int(lines[i+offset].split(" ")[dim+2])
        picked_indices.append((a,b,c))

if dim==2:
    if method==0:
        prepare_initial_grid_for_indices(all_points,n,m)
        prepare_initial_grid_for_leja_points(n,m)
    else:
        prepare_initial_grid_for_middle_points()

if dim==2:
    plot_picked_points_progressively_2d(picked_indices, picked_points,s,dt)
else:
    plot_picked_points_progressively_3d(picked_indices, picked_points,s,dt)
plt.show()
