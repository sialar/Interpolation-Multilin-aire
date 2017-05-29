from mpl_toolkits.mplot3d import Axes3D # librairie 3D
import matplotlib.pyplot as plt
import numpy as np
import sys


fig = plt.figure()

error_file = open( "../data/interpolation_error.txt", "r")
lines = error_file.readlines()
error_file.close()

n = int(lines[0])
x, y = [], []

for i in range(n):
    x.append(float(lines[i+1].split(" ")[0]))
    y.append(float(lines[i+1].split(" ")[1]))

plt.plot(x,y,c='g')
plt.title('Interpolation error')
plt.savefig('../images/interpolation_error.png')
plt.show()
