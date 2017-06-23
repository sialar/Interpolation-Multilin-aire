# -*- coding: utf-8 -*-

from scipy.optimize import fsolve
#import pylab
import numpy as np
import sys
import warnings
import scipy
import math
from scipy.integrate import quad
from math import *
import matplotlib.pyplot as plt

if (len(sys.argv)>2):
    f = int(sys.argv[2])
if f!=1 and f!=2 and f!=3 and f!=4 and f!=5:
    print "Function number must be 1, 2, 3, 4 or 5"
    sys.exit(1)

perturbations, perturbationsWidth = [], []
nbPerturbations = 0

def polynomialFunction(x,n):
    coefs_file = open( "../AI/data/polynomial_coefs.txt", "r")
    lines = coefs_file.readlines()
    coefs_file.close()
    deg = int(lines[0].split(" ")[1])
    d = int(lines[0].split(" ")[2])
    coefs = []
    res = 0

    for j in range(deg):
        for k in range(d):
            coefs.append(float(lines[j+1].split(" ")[k]))
    for j in range(deg):
        for k in range(d):
                res = res + coefs[j*d+k] * (x[k])**(j+1);
    return res

def readPerturbationFromFile():
    perturb_file = open( "../AI/data/perturbations.txt", "r")
    lines = perturb_file.readlines()
    perturb_file.close()
    nbPerturbations = int(lines[0])
    for i in range(nbPerturbations):
        perturbations.append(lines[1].split(" ")[i])
        perturbationsWidth.append(lines[2].split(" ")[i])

def hat(a,b,fa,fb,t):
    c = (a+b)/2
    if t <= c and t >= a:
        return ((t-a)/(c-a)+fa)
    elif t >= c and t <= b:
        return ((t-b)/(c-b)+fb)
    else:
        return 0

def g(x,n):
    res = 1
    for i in range(1,n):
        res *= exp(-x[i])
    return sqrt(1-x[0]**2) * res

def cosinus(x,n):
    f = 10
    res = 1
    for i in range(1,n):
        res *= exp(-x[i])
    return cos(2*math.pi*f*x[0]) * res

def h(x,n):
    res = 0
    for i in range(n):
        res += x[i]
    if n>1:
        return sin(res)
    else:
        return 2*exp(-x[0]**2 + sin(2*math.pi*x[0]))

def inOneHat(x):
    for i in range(nbPerturbations):
        xc = perturbations[i];
        xe = perturbationsWidth[i];
        if x>=xc-xe/2 and x<=xc+xe/2:
            return i
    return -1


def phi(x,n):
    res = 1
    f = 1
    hatPos = inOneHat(x[0])
    for i in range(1,n):
        res *= exp(-x[i])
    if (hatPos>-1):
        xc = perturbations[hatPos];
        xe = perturbationsWidth[hatPos];
        a = xc - xe/2;
        b = xc + xe/2;
        return hat(a,b,sin(2*M_PI*f*a),sin(2*M_PI*f*b),x(0)) * res;
    else:
        return sin(2*math.pi*f*x[0]) * res;

def sinOfNorm(x,n):
    res = 0
    for i in range(n):
        res += x[i]**2
    return sin(sqrt(res))

def getInterpolatedFunction(x,n):
    if f==1:
        return polynomialFunction(x,n)
    elif f==2:
        return g(x,n)
    elif f==3:
        return sinOfNorm(x,n)
    elif f==4:
        return h(x,n)
    else:
        return phi(x,n)

#z = []
#p = []
#y = np.linspace(-1.0, 1.0, num=1000)
#for i in range(len(y)):
#    p.append(y[i])
#    z.append(phi(p,1))
#    del p[:]
#plt.plot(y,z)
#plt.ylabel('some numbers')
#plt.show()


class function:
    def __init__(self, argsF, functionName, dimension):
        self.args = argsF
        self.name = functionName
        self.nbrArgs = len(argsF)
        self.dimension = dimension

        readPerturbationFromFile()

    def evaluate(self, x):
        a=0.25
        b=0.02
        c=0.05
        d= 0.1
        f = 70.0
        return getInterpolatedFunction(x, self.dimension)
        #return h(x)
        #return sinOfNorm(x)
        #return polynomialFunction(x)
	    #return (c/((x[0]-a)*(x[0]-a)+b + np.sin(x[0]*np.cos(x[1]+x[2]))/f + 0.05)+d/(np.sqrt(x[0]+x[1]+x[2])+1))*(x[1]*x[1]  + 0.05)*(x[2] + 0.01) + (self.args[1]*x[1] + self.args[2]*x[2])*np.exp(x[0]*x[0]+x[1]+x[2])/10 + (np.cos(x[0]))/10
