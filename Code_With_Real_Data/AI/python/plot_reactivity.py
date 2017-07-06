import matplotlib.pyplot as plt
import numpy as np
import string
import math
import sys

taille = (10,10)
fig = plt.figure(figsize=taille)

core = sys.argv[1]
figureTitle = "Comparaison of approximation errors (pcm) for the reactivity"
fileName = "../data/" + core + "/ReactivityError"
figureName = "../images/" + core + "/ReactivityError"
results_file = open(fileName,"r")
lines = results_file.readlines()
results_file.close()

point = []
co_err, tu_err, ai_err = [], [], []

for i in range(len(lines)):
    point.append(int(lines[i].split(" ")[0]))
    co_err.append(float(lines[i].split(" ")[9]))
    tu_err.append(float(lines[i].split(" ")[10]))
    ai_err.append(float(lines[i].split(" ")[12]))

s = 10
cocagne = plt.scatter(point,co_err,color='b',s=s,marker='o',alpha=0.1)
tucker = plt.scatter(point,tu_err,color='g',s=s,marker='o',alpha=0.1)
ai = plt.scatter(point,ai_err,color='r',s=s,marker='o',alpha=0.1)

plt.legend((cocagne, tucker, ai),
           ('Cocagne', 'Tucker', 'AI'),
           scatterpoints=1,
           loc='lower left',
           ncol=3,
           fontsize=10)

plt.title(figureTitle)
plt.savefig(figureName)

#plt.show()
