from mpl_toolkits.mplot3d import Axes3D # librairie 3D
import matplotlib.pyplot as plt
import numpy as np
import sys

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
    plt.subplot(121)
    plt.axis([-1.1, 1.1, -1.1, 1.1])
    plt.scatter(dichotomy_x, dichotomy_y, s=s, c='w')
    plt.subplot(122)
    plt.axis([-1.1, 1.1, -1.1, 1.1])
    plt.scatter(dichotomy_x, dichotomy_y, s=s, c='w')


def plot_picked_points_progressively_2d(alpha1, alpha2, points1, points2, s, dt):
    plt.hold(True)
    for i in range(nb_points):
        for j in range(i):
            q1 = points1[j]
            plt.subplot(121)
            plt.axis([-1.1, 1.1, -1.1, 1.1])
            plt.scatter(q1[0], q1[1], s=s, c='y')
            q2 = points2[j]
            plt.subplot(122)
            plt.axis([-1.1, 1.1, -1.1, 1.1])
            plt.scatter(q2[0], q2[1], s=s, c='y')

        p1 = points1[i]
        plt.subplot(121)
        plt.axis([-1.1, 1.1, -1.1, 1.1])
        plt.scatter(p1[0], p1[1], s=s, c='k')
        plt.title(alpha1[i])

        p2 = points2[i]
        plt.subplot(122)
        plt.axis([-1.1, 1.1, -1.1, 1.1])
        plt.scatter(p2[0], p2[1], s=s, c='k')
        plt.title(alpha2[i])

        plt.pause(dt)


def plot_picked_points_progressively_3d(alpha1, alpha2, points1, points2, s, dt):
    ax1 = fig.add_subplot(1, 2, 1, projection='3d')
    ax1.set_xlim(-1.1, 1.1)
    ax1.set_ylim(-1.1, 1.1)
    ax1.set_zlim(-1.1, 1.1)
    ax2 = fig.add_subplot(1, 2, 2, projection='3d')
    ax2.set_ylim(-1.1, 1.1)
    ax2.set_xlim(-1.1, 1.1)
    ax2.set_zlim(-1.1, 1.1)
    plt.hold(True)
    for i in range(nb_points):
        for j in range(i):
            q1 = points1[j]
            ax1.scatter(q1[0], q1[1], q1[2], c='y', s=s, marker='o')
            q2 = points2[j]
            ax2.scatter(q2[0], q2[1], q2[2], c='y', s=s, marker='o')
        p1 = points1[i]
        ax1.scatter(p1[0], p1[1], p1[2], c='k', s=s, marker='o')
        #ax1.set_title(alpha1[i])
        p2 = points2[i]
        ax2.scatter(p2[0], p2[1], p2[2], c='k', s=s, marker='o')
        #ax2.set_title(alpha2[i])
        plt.pause(dt)


all_points, picked_points1, picked_points2 = [], [], []

if len(sys.argv)!=3:
    print("Invalid number of arguments: choose the space dimension and the interpolation method.")
    sys.exit(0)
else:
    dim = int(sys.argv[1])
    method = int(sys.argv[2])

taille = (16,8)
fig = plt.figure(figsize=taille)

path_file1 = open( "../data/path.txt", "r")
lines1 = path_file1.readlines()
path_file1.close()

path_file2 = open( "../data/other_path.txt", "r")
lines2 = path_file2.readlines()
path_file2.close()

alpha1, alpha2 = [], []
nb_points = int(lines1[1])
offset = 2
s = 50
dt = 0.1


for i in range(nb_points):
    x1 = float(lines1[i+offset].split(" ")[0])
    y1 = float(lines1[i+offset].split(" ")[1])
    z1 = 0
    x2 = float(lines2[i+offset].split(" ")[0])
    y2 = float(lines2[i+offset].split(" ")[1])
    z2 = 0
    if dim==3:
        z1 = float(lines1[i+offset].split(" ")[2])
        z2 = float(lines2[i+offset].split(" ")[2])
    picked_points1.append((x1,y1,z1))
    picked_points2.append((x2,y2,z2))
    alpha1.append(float(lines1[nb_points+offset].split(" ")[i]))
    alpha2.append(float(lines2[nb_points+offset].split(" ")[i]))

if dim==2:
    prepare_initial_grid_for_middle_points()
    plot_picked_points_progressively_2d(alpha1, alpha2, picked_points1, picked_points2, s, dt)
else:
    plot_picked_points_progressively_3d(alpha1, alpha2, picked_points1, picked_points2, s, dt)

plt.show()
