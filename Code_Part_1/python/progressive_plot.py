import numpy as np
import matplotlib.pyplot as plt

def prepare_initial_grid(points,n,m,s,c):
    plt.axis([-1, n, -1, m])
    x, y = [], []
    for i in points:
        x.append(i[0])
        y.append(i[1])
    plt.scatter(x, y, s=s, c=c)

def plot_picked_points_progressively(points,s,c,dt):
    plt.ion()
    for p in points:
        plt.scatter(p[0], p[1], s=s, c=c)
        plt.title('Progression du chemin d\'indices')
        plt.pause(dt)
    while True:
        plt.pause(dt)

all_points, picked_points = [], []

ai_output_file = open( "path.txt", "r")
lines = ai_output_file.readlines()
ai_output_file.close()

n = int(lines[0].split(" ")[0])
m = int(lines[0].split(" ")[1])
nb_points = int(lines[1])
offset = 2
s = 100

for i in range(n):
    for j in range(m):
        all_points.append((i,j))

for i in range(nb_points):
    x = int(lines[i+offset].split(" ")[0])
    y = int(lines[i+offset].split(" ")[1])
    picked_points.append((x,y))

prepare_initial_grid(all_points,n,m,s,'g')
plot_picked_points_progressively(picked_points,s,'r',0.05)
