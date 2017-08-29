# -*- coding: utf-8 -*-

from scipy.optimize import fsolve
#import pylab
import numpy as np
import warnings
import scipy
from scipy.integrate import quad
class function:
    def __init__(self, argsF, functionName, dimension):
        self.args = argsF
        self.name = functionName
        self.nbrArgs = len(argsF)
        self.dimension = dimension

    
    def evaluate(self, x):
        #return self.args[0]*x[0]*x[0]*x[1] + self.args[1]*x[1] + self.args[2]*x[2] + 1
        a=0.25
        b=0.02
        c=0.05
        d= 0.1 #0.3#0.1
        f = 70.0

	return (c/((x[0]-a)*(x[0]-a)+b + np.sin(x[0]*np.cos(x[1]+x[2]))/f + 0.05)+d/(np.sqrt(x[0]+x[1]+x[2])+1))*(x[1]*x[1]  + 0.05)*(x[2] + 0.01) + (self.args[1]*x[1] + self.args[2]*x[2])*np.exp(x[0]*x[0]+x[1]+x[2])/10 + (np.cos(x[0]))/10
 

  