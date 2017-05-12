import numpy as np
import matplotlib.pyplot as plt
import sys

input_file = open( "plot_function.txt", "r")
lines = input_file.readlines()
input_file.close()

n = int(lines[0])
x, y, z, t = [], [], [], []

for i in range(1,n):
    x.append(float(lines[i].split(" ")[0]))
    y.append(float(lines[i].split(" ")[1]))
    z.append(float(lines[i].split(" ")[2]))
    t.append(0)


plt.plot(x, y, 'k', c='r')
plt.plot(x, z, 'k', c='b')
plt.plot(x, t, 'k', c='g')

plt.savefig('lagrange_polynomials.png')
plt.axis([-1, 1, -2, 2])
plt.show()
