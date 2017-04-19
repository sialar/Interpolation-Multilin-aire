from mpl_toolkits.mplot3d import Axes3D # librairie 3D
import matplotlib.pyplot as plt
import numpy as np

fig = plt.figure()
ax = fig.gca(projection='3d')

data_file = open("interpolation_result_2D.txt", "r")
lines = data_file.readlines()
data_file.close()

xSize = int(lines[0].split(" ")[0])
ySize = int(lines[0].split(" ")[1])

data = []
x_tab, y_tab, z_tab, real_z_tab = [], [], [], []

for i in range(1,xSize*ySize+1):
	x = float(lines[i].split(" ")[0])
	y = float(lines[i].split(" ")[1])
	z = float(lines[i].split(" ")[2])
	real_z = float(lines[i].split(" ")[3])
	data.append((x,y,z,real_z))

data = sorted(data, key=lambda col: col[0])

for i in range(len(data)):
	x_tab.append(data[i][0])
	y_tab.append(data[i][1])
	z_tab.append(data[i][2])
	real_z_tab.append(data[i][3])

ax.plot_trisurf(x_tab, y_tab, z_tab, linewidth=0.01, antialiased=True, color='b')
ax.plot_trisurf(x_tab, y_tab, real_z_tab, linewidth=0.01, antialiased=True, color='r')
plt.show()
