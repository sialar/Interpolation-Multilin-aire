from __future__ import division
from scipy import *
from pylab import *

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D # librairie 3D

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

data_file = open("interpolation_result.txt", "r")
lines = data_file.readlines()
data_file.close()

xSize = int(lines[0].split(" ")[0])
ySize = int(lines[0].split(" ")[1])

x_fin, y_fin, z_fin, real_z_fin = 0, 0, 0, 0
x,y,z,real_z = [], [], [], []

for i in range(1,xSize*ySize+1):
	x.append(float(lines[i].split(" ")[0]))
	y.append(float(lines[i].split(" ")[1]))
	z.append(float(lines[i].split(" ")[2]))
	real_z.append(float(lines[i].split(" ")[3]))

ax.plot_trisurf(x, y, z, linewidth=0.01, antialiased=True, color='r')
ax.plot_trisurf(x, y, real_z, linewidth=0.01, antialiased=True, color='b')
plt.show()
