from mpl_toolkits.mplot3d import Axes3D # librairie 3D
import matplotlib.pyplot as plt
import numpy as np

def prepare_initial_grid(points,n,m,l,s,c):
    x, y, z = [], [], []
    for i in points:
        x.append(i[0])
        y.append(i[1])
        z.append(i[2])
    i = 0
    while (i<len(z) and z[i]==0):
        i = i+1
    if i==len(z):
        plt.axis([-1, n, -1, m])
        plt.scatter(x, y, s=s, c=c)
        return 2
    else:
        ax = fig.gca(projection='3d')
        ax.scatter(x, y, z, s=s, c=c)
        return 3


def plot_picked_points_progressively(dim,points,s,c,dt):
    plt.ion()
    for p in points:
        if dim==2:
            plt.scatter(p[0], p[1], s=s, c=c)
        else:
            ax = fig.gca(projection='3d')
            ax.scatter(p[0], p[1], p[2], s=s, c=c)
        plt.pause(dt)
    while True:
        plt.pause(dt)

all_points, picked_points = [], []

fig = plt.figure()

ai_output_file = open( "path.txt", "r")
lines = ai_output_file.readlines()
ai_output_file.close()

dim = int(lines[0])
n = int(lines[1].split(" ")[0])
m = int(lines[1].split(" ")[1])
l = int(lines[1].split(" ")[2])
nb_points = int(lines[2])
offset = 3
s = 100

for i in range(n):
    for j in range(m):
        for k in range(l):
            all_points.append((i,j,k))

for i in range(nb_points):
    x = int(lines[i+offset].split(" ")[0])
    y = int(lines[i+offset].split(" ")[1])
    z = int(lines[i+offset].split(" ")[2])
    picked_points.append((x,y,z))

dim = prepare_initial_grid(all_points,n,m,l,s,'y')
plot_picked_points_progressively(dim, picked_points,s,'k',0.1)
plt.show()
