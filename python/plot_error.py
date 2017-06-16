from mpl_toolkits.mplot3d import Axes3D # librairie 3D
import matplotlib.pyplot as plt
import numpy as np
import sys


fig = plt.figure()

error_file = open( "../data/interpolation_error.txt", "r")
lines = error_file.readlines()
error_file.close()

n = int(lines[0])
x, y0, y1, y2, y3 = [], [], [], [], []

for i in range(n):
    x.append(float(lines[i+1].split(" ")[0]))
    y0.append(float(lines[i+1].split(" ")[1]))
    y1.append(float(lines[i+1].split(" ")[2]))
    y2.append(float(lines[i+1].split(" ")[3]))
    y3.append(float(lines[i+1].split(" ")[4]))

for i in range(int(sys.argv[1])):
    x.pop(0)
    y0.pop(0)
    y1.pop(1)
    y2.pop(2)
    y3.pop(2)

plt.plot(x,y0,c='g')
plt.plot(x,y1,c='r')
plt.plot(x,y2,c='b')
plt.plot(x,y3,c='k')
plt.title('Interpolation error')
plt.savefig('../images/interpolation_error.png')
plt.show()
