# -*- coding: utf-8 -*-

import LagrangePolynomial 
import LagrangeInterpolation as LI
import numpy as np
from numpy import linalg as LA
import scipy
import warnings
  
###===================================================
### Greedy algorithm via a discret version
### ---> Determine a set of interpolation points from a set of given points for a set of vectors (represent for a set of functions) 
###===================================================

def getEmpiricalPoints(listOfFcts, listOfInitialInterpolationPoints, orderInterpolation):  

   
    listOfFinalIndexOfEmpiricalInterpolationPoints = []
    listOfFinalEmpericalInterpolationPoints = []
    
    for m in range(orderInterpolation):
        
        
        listOfResiduaForInterpolation_m = []
        
        listOfInterpolFcts = []
        
        if m == 0:               
            indexOfArgMaxPoint, indexOfFctMax = getMaxOfAListOfFcts(listOfFcts)
            argMaxPoint = listOfInitialInterpolationPoints[indexOfArgMaxPoint]
            
            listOfFinalIndexOfEmpiricalInterpolationPoints.append(indexOfArgMaxPoint)
           
            listOfFinalEmpericalInterpolationPoints.append(argMaxPoint)
        
            for indexFct in range(len(listOfFcts)): 
                basisFctAtArgMaxPoint_m_indexFct = listOfFcts[indexFct][indexOfArgMaxPoint]
                interpolFct_m_indexFct = basisFctAtArgMaxPoint_m_indexFct * (listOfFcts[indexOfFctMax]/listOfFcts[indexOfFctMax][indexOfArgMaxPoint])              
                residua_m_indexFct = listOfFcts[indexFct] - interpolFct_m_indexFct
                listOfResiduaForInterpolation_m.append(residua_m_indexFct)   
            
        if m > 0:               
            indexOfArgMaxPoint, indexOfFctMax = getMaxOfAListOfFcts(listOfResiduaForInterpolation_m_previous)
            argMaxPoint = listOfInitialInterpolationPoints[indexOfArgMaxPoint]  
            
            listOfFinalIndexOfEmpiricalInterpolationPoints.append(indexOfArgMaxPoint)            
            listOfFinalEmpericalInterpolationPoints.append(argMaxPoint)
            
            for indexFct in range(len(listOfFcts)):
                listOfFinalIndexOfEmpiricalInterpolationPoints_array = np.asarray(listOfFinalIndexOfEmpiricalInterpolationPoints)
                listOfInitialInterpolationPoints_array = np.asarray(listOfFinalEmpericalInterpolationPoints)
                
                
                
                valuesOfBasisFctOnInterpolationPoints = listOfFcts[indexFct][listOfFinalIndexOfEmpiricalInterpolationPoints_array]
                
                Lagrangefct = [LagrangePolynomial.LagrangePolynomial(listOfInitialInterpolationPoints_array, valuesOfBasisFctOnInterpolationPoints)]
                
                interpolFct_m_indexFct = LI.getInterpolationArr(Lagrangefct, listOfInitialInterpolationPoints)  
                residua_m_indexFct = listOfFcts[indexFct] - interpolFct_m_indexFct
                listOfResiduaForInterpolation_m.append(residua_m_indexFct)
                
        listOfResiduaForInterpolation_m_previous = listOfResiduaForInterpolation_m  ### Exist after the first loop m=0
    

    return listOfFinalEmpericalInterpolationPoints, listOfFinalIndexOfEmpiricalInterpolationPoints
    

###===================================================
### Find a max and argmax of a set of vectors (represents for a set of functions in sense discret)
###===================================================
def getMaxOfAListOfFcts(listOfFcts):
    n = len(listOfFcts)
    listOfMaxValueFcts = []
    listOfArgmaxIndex = []
    for i in range(n):
        absFct = abs(listOfFcts[i])
       
        maxValue = np.max(absFct)
        argMaxIndex = np.argmax(absFct)
        
        if absFct[argMaxIndex] != maxValue:
            warnings.warn("warning",'not right the argmax')
       
        listOfMaxValueFcts.append(maxValue)
        listOfArgmaxIndex.append(argMaxIndex)
        
        IndexesOfArgMaxPoint = [ii for ii in range(len(absFct)) if absFct[ii]== maxValue]
        
            
    finalMaxValue = max(listOfMaxValueFcts)
    indexFct = np.argmax(listOfMaxValueFcts)
    indexOfArgMaxPoint = listOfArgmaxIndex[indexFct]
    IndexesOfArgMaxPoint = [ii for ii in range(len(listOfMaxValueFcts)) if listOfMaxValueFcts[ii]== finalMaxValue]

    return indexOfArgMaxPoint, indexFct


