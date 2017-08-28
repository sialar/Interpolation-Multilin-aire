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
    """
        focntion qui fournit la liste des points communs pour toutes les sections efficaces qui sont utilisés dans la résolution des coefficients de Tucker 
        
        Argument:
           listOfFcts [une liste d'arrays]: chaque array réprésente un vecteur propre gardé pour une direction donnée qui est utilisé dans la décomposition de Tucker  
           listOfInitialInterpolationPoints [une liste ou array des flottants] : coordonées des points sur lesquels on sélectionne les points par greedy
           orderInterpolation [un entier] : un ordre de l'interpolation de Lagrange et aussi le nombre de points que l'on veut sélectionner 
       retourne:
           listOfFinalEmpericalInterpolationPoints [une liste ou array des flottants] : coordonées des points sélectionnés par greedy
           listOfFinalIndexOfEmpiricalInterpolationPoints [une liste ou array des entiers] : indices des points sélectionnés par greedy
    """
   
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
def getMaxOfAListOfFcts(listOfFcts):
    """
        focntion qui fournit le max et argmax d'un ensemble des vecteurs de même taille (représente un ensemble des fonctions dans le sens discret)
               
        
        Argument:
           listOfFcts [une liste d'arrays]: chaque array réprésente un vecteur propre gardé pour une direction donnée qui est utilisé dans la décomposition de Tucker  
          
       retourne:
           indexOfArgMaxPoint [un entier] : indice du point sur lequel on a max sur toutes les vecteurs (en valeurs absolutes)
           indexFct [un entier] : indice d'une fonction dans l'ensemble sur lequel on atteint le max
    """
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


