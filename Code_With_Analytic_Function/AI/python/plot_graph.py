from mpl_toolkits.mplot3d import Axes3D # librairie 3D
import matplotlib.pyplot as plt
import numpy as np
import sys


fig = plt.figure()

iters = [ 500, 1000, 2000, 3500, 5000, 7500, 10000 ]

ai_3 = [ 1677.8, 330.07, 25.73, 10.6, 16.8, 0.501, 0.246]
ai_3_time = [ 0.16, 0.51, 2.04, 7.2, 13.19, 32.69, 62.03 ]

tu_3 = [ 1512.49, 387.62, 73.9, 20.05, 3.66, ]
tu_3_time = [ 0.03, 0.16, 9.26, 50.65, 154.05, ]


ai_5 = [ 21519, 16827, 8587.88, 4297.72, 1366.71, 567.824]
ai_5_time = [ 0.31, 1.02, 4.6, 15.7, 35.42, 74.56 ]

tu_5 = [ ]
tu_5_time = [ ]


print len(iters[6:]), len(tu_time)
plt.plot(iters, ai_time, c='r', alpha = 0.4)
plt.plot(iters_plus, ai_plus_time, c='r')
plt.plot(iters[6:], tu_time, c='g')
plt.xlabel('Number of evaluation point')
plt.ylabel('Run time (s)')

plt.show()
