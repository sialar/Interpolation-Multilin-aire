# -*- coding: utf-8 -*-
import numpy as np
import warnings
from numpy import linalg as LA
import scipy
import multiDimensIntegralMethods as M 

###=====================================
###=====================================
### This class is used to orthonormalize a set of vecteurs with norm L2 
###=====================================
### INPUTS:
### Set of vectors U
### Integral weights for scalar product
###=====================================
### OUTPUTS:
### Set of orthonormalized vectors 
##=====================================
##=====================================

class orthoNormalL2 :
  ###=====================================
  ### This class is used to orthonormalize a set of vecteurs with norm L2 
  ###=====================================
    def  __init__(self, W_xk) :
        """
       Argumentes: 
        "W_xk" : poid d'intégration pour le produit scalaire dans L2
        """
        self.W_xk = W_xk        
##=====================================
    def dotL2(self, u1, u2):
        """
       Argumentes: 
        "u1" [array] : représente un vecteur
        "u2" [array] : représente un vecteur
       Retourn :
        np.sum(u1*u2* self.W_xk) [un flottant] : produit scalaire dans L2 de u1 et u2
        """

        if len(u1) != len(u2) or  len(u1) != len( self.W_xk)  :
            warnings.warn("error",'vector u1, u2 and  weights must be same size')
        else :
            return  np.sum(u1*u2* self.W_xk) 
##=====================================
    def normL2(self, u):
        """
       Argumentes: 
        "u" [array] : représente un vecteur
      
       Retourn :
        np.sqrt(self.dotL2(u,u)) [un flottant] : la norme L2 de u
        """
        return np.sqrt(self.dotL2(u,u))
##=====================================
    def normalizeL2(self, u ):
        """
       Argumentes: 
        "u" [array] : représente un vecteur
      
       Retourn :
       u/ self.normL2(u) [un array] : la normalisation de vecteur u
        """
       
        return  u/ self.normL2(u)
##=====================================  
    def projU1_U2(self, u1, u2):
        """
       Argumentes: 
        "u1" [array] : représente un vecteur
        "u2" [array] : représente un vecteur
       Retourn :
        projectU1U2 [array] : la projection L2 de u1 sur u2
        """
        dotU1U2 = self.dotL2(u1, u2)
        dotU1U1 = self.dotL2(u1, u1) 
        projectU1U2 =  (dotU1U2/dotU1U1)*u1        
        return  projectU1U2

    def orthonormalVecteurs (self, U):
##=====================================
# Source http://en.wikipedia.org/wiki/Gram%E2%80%93Schmidt_process
##=====================================
        """
     Argumentes: 
       U [une liste de arrays] : liste des vecteurs à orthonormer
     Retourn :
      U [une liste de arrays] : liste des vecteurs orthonormés
        """

        k = len(U)
       
        for i in range(0,k):           
            U[i] = self.normalizeL2(U[i])           
            for j in range(i+1, k):               
                U[j] = U[j] -  self.projU1_U2(U[i], U[j])                 
        
        return U
        
##=====================================  