from mpl_toolkits.mplot3d import Axes3D # librairie 3D
import matplotlib.pyplot as plt
import numpy as np
import sys

def plot_picked_points_progressively_2d(alpha, indices,points,s,dt):
    #plt.hold(True)
    plt.axis([-1.1, 1.1, -1.1, 1.1])
    for i in range(len(points)):
        p = points[i]
        plt.scatter(p[0], p[1], s=s/2, c='r')
        plt.title(alpha[i])
        plt.pause(dt)

def plot_picked_points_progressively_3d(alpha, indices,points,s,dt):
    ax = fig.gca(projection='3d')
    ax.set_xlim(-1.1, 1.1)
    ax.set_ylim(-1.1, 1.1)
    ax.set_zlim(-1.1, 1.1)
    #plt.hold(True)
    for i in range(len(points)):
        p = points[i]
        ax.scatter(p[0], p[1], p[2], c='r', s=s, marker='o')
        ax.set_title(alpha[i])
        plt.pause(dt)


all_points, picked_points, picked_indices = [], [], []

if len(sys.argv)!=2:
    print("Invalid number of arguments, choose the space dimension.")
    sys.exit(0)
else:
    dim = int(sys.argv[1])

taille = (16,8)
fig = plt.figure(figsize=taille)

ai_output_file = open( "../data/path.dat", "r")
lines = ai_output_file.readlines()
ai_output_file.close()

if len(lines)==0:
    sys.exit(1)

l = 1
n = int(lines[0].split(" ")[0])
m = int(lines[0].split(" ")[1])
if dim==3:
    l = int(lines[0].split(" ")[2])

offset = 2
s = 80
dt = 0.05

alpha = []

nb_points = int(lines[1])
for i in range(nb_points):
    x = float(lines[i+offset].split(" ")[0])
    y = float(lines[i+offset].split(" ")[1])
    z = 0
    if dim==3:
        z = float(lines[i+offset].split(" ")[2])
    picked_points.append((x,y,z))
    alpha.append(lines[nb_points+offset].split(";")[i])

if dim==2:
    plot_picked_points_progressively_2d(alpha, picked_indices, picked_points, s, dt)
else:
    plot_picked_points_progressively_3d(alpha, picked_indices, picked_points, s, dt)
plt.show()
