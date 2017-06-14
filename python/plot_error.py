from mpl_toolkits.mplot3d import Axes3D # librairie 3D
import matplotlib.pyplot as plt
import numpy as np
import sys


fig = plt.figure()

error_file = open( "../data/interpolation_error.txt", "r")
lines = error_file.readlines()
error_file.close()

n = int(lines[0])
x, y0, y1, y2 = [], [], [], []

for i in range(n):
    x.append(float(lines[i+1].split(" ")[0]))
    y0.append(float(lines[i+1].split(" ")[1]))
    if len(lines[i+1].split(" ")) > 2:
        y1.append(float(lines[i+1].split(" ")[2]))
    if len(lines[i+1].split(" ")) > 3:
        y2.append(float(lines[i+1].split(" ")[3]))



plt.plot(x,y0,c='g')
if len(y1) > 0:
    plt.plot(x,y1,c='r')
if len(y2) > 0:
    plt.plot(x,y2,c='b')
plt.title('Interpolation error')
plt.savefig('../images/interpolation_error.png')
plt.show()
