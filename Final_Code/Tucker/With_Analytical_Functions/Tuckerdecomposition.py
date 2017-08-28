# -*- coding: utf-8 -*-
import numpy as np
import KLdecomposition as KL
import LagrangeInterpolation as LI
import random
import warnings
import copy
##=====================================
### This class is used to define the Tucker decomposition for a given cross-section (via crossSectionName, e.g: macro_totale0)
##=====================================
### INPUTS:
### + List of all basis functions for "d" directions (defined by KLdecomposition class)
### + Reference points (determined by Greedy algorithm) to solve system of coefficients
##=====================================
### OUTPUTS:
### + Coefficients of the Tucker decomposition of the given "crossSectionName"
### + Tucker decomposition which can evaluate any given point
##=====================================
###=====================================
### TO CHANGE (following the test case): 
### listOfKeptPercent for a priori sparsity
###=====================================
class Tucker:
    def __init__(self, listOfBasisFcts, listOfNbOfBasisFcts, refCrossSectionToDetermineCoeffs, listOfInterpolationPoints, crossSectionName, listOfKeptEigenValues):
        
        ### Data for the given cross-section via its "crossSectionName"
        self.listOfBasisFcts = listOfBasisFcts
        self.listOfNbOfBasisFcts  = listOfNbOfBasisFcts 
        self.refCrossSectionToDetermineCoeffs = refCrossSectionToDetermineCoeffs      
        self.listOfInterpolationPoints = listOfInterpolationPoints
        self.crossSectionName = crossSectionName
       
        self.dimension = len(self.listOfBasisFcts)
        
        self.listOfKeptEigenValues = listOfKeptEigenValues
        #print "listOfKeptEigenValues", self.listOfKeptEigenValues #[str(0)] #[0]
        
        #raw_input()
      
        ###=====================================    
        CardinalOfBasisFcts = 1
        for Axis_k in range(self.dimension):
            CardinalOfBasisFcts = CardinalOfBasisFcts*listOfNbOfBasisFcts[Axis_k]
            
        ###===================================== 
        self.CardinalOfBasisFcts = CardinalOfBasisFcts
 ##===================================== 
    def produitArr(self, u):    
        n = len(u)   
        s = 1.0
        for i in range(0,n):
            s = s*u[i]
        return s
##=====================================    
    def transfIndexIntToArray(self, I, d, listOfNbPoints):  
        I = int(I)           
        i = []
        ii = [0]
        for k in range(0,d):
            i.append(ii)
                
        index = range(0,d)
        indexInverse = index[::-1] 
            
        for j in indexInverse:
            if j > 0 :            
                i[j] = int(I/self.produitArr(listOfNbPoints[:j] ))
                I = I - i[j]*self.produitArr(listOfNbPoints[:j] )
                I = int(I)
            
            else :
                i[j] = int(I) 
            
        return i

 ##=====================================
    def lastSquare(self, A, b):           
        sol_np = np.linalg.lstsq(A, b)[0]
        return sol_np
 ##=====================================
    def getCoefficients(self):

        self.matrixA_full = np.zeros(self.CardinalOfBasisFcts*self.CardinalOfBasisFcts).reshape(self.CardinalOfBasisFcts,self.CardinalOfBasisFcts)
        self.vectorOfExactValues = np.zeros(self.CardinalOfBasisFcts) 
        ###==============================================
        ### self.listOfRefPoint used to determine cross-section values in b (on the right hand size of the system A*coeffs = b)
        self.listOfRefPoint = []
        ###==============================================
        for i in range(self.CardinalOfBasisFcts):
            ### ==================== Determine exact value :f(x1,x2,...xd)====================
            iArr = self.transfIndexIntToArray(i, self.dimension, self.listOfNbOfBasisFcts)           
            refPoint = []
           
            for Axis_k in range(self.dimension):
               
                refValue_Axis_k = self.listOfInterpolationPoints[str(Axis_k)][iArr[Axis_k]]
                refPoint.append(refValue_Axis_k)
            
            self.vectorOfExactValues[i]  = self.refCrossSectionToDetermineCoeffs.evaluate(refPoint)
            
            self.listOfRefPoint.append(refPoint)
           
            ### ========================================
            ### Determine the elements of matrix "A" via the evaluation of eigenvectors on reference points
            ### ========================================
           
            for j in range(self.CardinalOfBasisFcts):
                p = 1.0
                jArr = self.transfIndexIntToArray(j, self.dimension, self.listOfNbOfBasisFcts)   
                for Axis_k in range(self.dimension):         
     
                    basisFct =  self.listOfBasisFcts[str(Axis_k)][jArr[Axis_k]]
                    valueOfEigenVector_Axis_k_jk = LI.getInterpolation(basisFct, refPoint[Axis_k])
                   
                    p = p*valueOfEigenVector_Axis_k_jk
                   
                self.matrixA_full[i][j] = p
        ###===================================================
        ### self.FinalTuckerCoeffs is an array, each index [I] associated with an array from indexes of basis functions: I = [i1,i2,...,id]
                                                                                                                                
        self.FinalTuckerCoeffs = np.linalg.solve(self.matrixA_full,self.vectorOfExactValues)
        ### Transform list in I to list in [i1,i2,...,id]
        self.listOfFinalCoefIndexes_arr = []
        for j in range(self.CardinalOfBasisFcts):
            jArr = self.transfIndexIntToArray(j, self.dimension, self.listOfNbOfBasisFcts)         
            self.listOfFinalCoefIndexes_arr.append(jArr)
        
        print "self.FinalTuckerCoeffs", self.FinalTuckerCoeffs
        
        
        File_indexes = open("indexCoef_initial.dat",'w')      
        
        
        for i in range(len(self.listOfFinalCoefIndexes_arr)):
            i0 = self.listOfFinalCoefIndexes_arr[i][0]
            i1 = self.listOfFinalCoefIndexes_arr[i][1]
            i2 = self.listOfFinalCoefIndexes_arr[i][2]
            
            File_indexes.write("%s\t%s\t%s\n"%(i0,i1,i2))
        File_indexes.close()
        #raw_input()
###===================================================
   
##===================================================
    def evaluate(self, point) :
        approxValue = 0.0
   
        for I in range(len(self.listOfFinalCoefIndexes_arr)) :
            d = 1.0
           
            for Axis_k in range(self.dimension) :                
                fct = self.listOfBasisFcts[str(Axis_k)][self.listOfFinalCoefIndexes_arr[I][Axis_k]]
                d = d*LI.getInterpolation(fct, point[Axis_k])
               
            approxValue = approxValue +  self.FinalTuckerCoeffs[I]*d       
     
        return approxValue
 
