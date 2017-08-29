# -*- coding: utf-8 -*-
import numpy as np
import warnings
from numpy import linalg as LA
import scipy
import multiDimensIntegralMethods as M 
import copy
###=====================================
###=====================================
### This class is used to orthonormalize a set of vecteurs with norm L2 
###=====================================
### INPUTS:
### Set of vectors U
### Integral weight for scalar product
###=====================================
### OUTPUTS:
### Set of orthonormalized vectors 
##=====================================
##=====================================

class orthoNormalL2 :
    def  __init__(self, W_xk) :
        self.W_xk = W_xk        
##=====================================
    def dotL2(self, u1, u2):

        if len(u1) != len(u2) or  len(u1) != len( self.W_xk)  :
            warnings.warn("error",'vector u1, u2 and  weights must be same size')
        else :
            return  np.sum(u1*u2* self.W_xk) 
##=====================================
    def normL2(self, u):
        return np.sqrt(self.dotL2(u,u))
##=====================================
    def normalizeL2(self, u ):
       
        return  u/ self.normL2(u)
##=====================================  
    def projU1_U2(self, u1, u2):
        dotU1U2 = self.dotL2(u1, u2)
        dotU1U1 = self.dotL2(u1, u1) 
        projectU1U2 =  (dotU1U2/dotU1U1)*u1        
        return  projectU1U2
##=====================================
# Source http://en.wikipedia.org/wiki/Gram%E2%80%93Schmidt_process
##=====================================
    def orthonormalVecteurs (self, U):
    ###U : containing vectors to orthonormalize.
        k = len(U)
       
        for i in range(0,k):           
            U[i] = self.normalizeL2(U[i])           
            for j in range(i+1, k):               
                U[j] = U[j] -  self.projU1_U2(U[i], U[j])         
        return U
        
##=====================================  
