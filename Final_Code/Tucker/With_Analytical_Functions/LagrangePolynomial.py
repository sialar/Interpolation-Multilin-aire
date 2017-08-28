# -*- coding: utf-8 -*-
import math
import numpy as np
import warnings
import copy
#import roots
##=====================================
### This class is used to define a Lagrange polynomial via given values of parameter x and corresponding function values f(x)
##=====================================
### INPUTS:
### + Set of {x} and {f(x)}
##=====================================
### OUTPUTS:
### + Lagrange polynomial
### + Different interpolation methods for : x (a point), x (an array)
##=====================================
class  LagrangePolynomial:
    def __init__(self, x, fx):
        if len(x) != len(fx) :
            warnings.warn("error", "x and fx must be the same size")
        else :
            self.degree = len(x)
            self.axisValues = x
            self.polynomialValues = fx
##===================================================             
    ### ======  x is a point ========            
          
    def getInterpolation(self, x ):
   #given points (self.axisValues, self.polynomialValues), fit a lagrange polynomial and return
        #the value at point x """

        f = 0.0
    
    # sum over points
        m = 0
        while (m < len(self.axisValues)):

        # create the Lagrange basis polynomial for point m        
            l = None

            n = 0
            while (n < len(self.axisValues)):
                if n == m:
                    n += 1
                    continue

                if l == None:
                    l = (x - self.axisValues[n])/(self.axisValues[m] - self.axisValues[n])
                else:
                    l *= (x - self.axisValues[n])/(self.axisValues[m] - self.axisValues[n])

                n += 1

        
            f += self.polynomialValues[m]*l

            m += 1

        return f
##===================================================         
    def getInterpolation_bis(self, x):
        n = len(self.axisValues)
        valueOfPolynome = 0.0
        
        for i in range(n):
            li = 1.0
            indexesOfPoints = range(n)
            indexesOfPoints.remove(indexesOfPoints[i])
            for j in indexesOfPoints:                                  
                li = li*(x - self.axisValues[j])/(self.axisValues[i]-self.axisValues[j])
                            
            valueOfPolynome = valueOfPolynome + self.polynomialValues[i]*li
        return valueOfPolynome            
##=================================================== 
##=================================================== 
    ### ======  x is an array, return an array ======  
    def getInterpolationArr_bis(self, x):
        fArr = []
        for xi in x :
            fArr.append(self.getInterpolation_bis(xi))
        return fArr
##=================================================== 

    def getInterpolationArr(self, x):
        fArr = []
        for xi in x :
            fArr.append(self.getInterpolation(xi))
        return fArr
        

##===================================================  
    def getDerivative(self, x):
        n = len(self.axisValues)
        polynome_derivative = 0.0
        
        for i in range(n):
           
            indexesOfPoints = range(n)
            indexesOfPoints.remove(indexesOfPoints[i])

            
            denominator_i = 1.0
            for j in indexesOfPoints:  
                denominator_i = denominator_i * (self.axisValues[i]-self.axisValues[j])
                
                
            numerator_i = 0.0
            
            for j in indexesOfPoints:  
                numerator_i_derivative_at_xj = 1.0
                indexesOfPoints_ij = indexesOfPoints[:]
                indexesOfPoints_ij.remove(j)

                for k in indexesOfPoints_ij:
                    numerator_i_derivative_at_xj = numerator_i_derivative_at_xj*(x - self.axisValues[k])
                    
                numerator_i  = numerator_i  + numerator_i_derivative_at_xj
                  
            polynome_derivative = polynome_derivative + self.polynomialValues[i]*numerator_i/denominator_i  
            
        return  polynome_derivative
 
