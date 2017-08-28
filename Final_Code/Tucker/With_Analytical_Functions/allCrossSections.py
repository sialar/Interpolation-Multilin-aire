# -*- coding: utf-8 -*-
import KLdecomposition as KL
import numpy as np
import greedyAlgo as greedy
import LagrangeInterpolation as LI
#import crossSections

##=====================================
### This method is used to gather all basis functions of all cross-sections, the main goal is to
### + Determine all common points per axis for all basis functions by the greedy algorithm 
##=====================================

def getAllBasisFctsForAllSections(listOfCrossSectionNamesWithArgs, listOfDomainBorders, listOfNbPointsOnEachIntervalForFinalDiscretization,\
            listOfTuckerGridNodesFromDKlibs, listOfMethodsIntegrals, listOfCrossSectionsForTuckerGridsFromDklibs):
            
    listOfBasisFcts = {}
    listOfNbOfBasisFcts = {}
    
    listOfKeptEigenValues = {}
    
    dimension = len(listOfDomainBorders)
    for csName in listOfCrossSectionNamesWithArgs:
        listOfBasisFcts_csName = {}
        listOfNbOfBasisFcts_csName = []
        listOfKeptEigenValues_csName = {}
        
        for Axis_k in range(dimension):
            KLDecomposition_csName_Axis_k = KL.KL(Axis_k, listOfDomainBorders, listOfNbPointsOnEachIntervalForFinalDiscretization[str(Axis_k)],\
            listOfTuckerGridNodesFromDKlibs[str(Axis_k)], listOfMethodsIntegrals[str(Axis_k)], listOfCrossSectionsForTuckerGridsFromDklibs[csName][str(Axis_k)], csName)
            ### Pour appeler la fct qui fournit toutes les fcts de base (orthonormée) dans la direction k
            KLDecomposition_csName_Axis_k.finalOrthoNormalEigens_Axis_k()
            KLDecomposition_csName_Axis_k.getListOfInterpolationFcts()
            
            listOfBasisFcts_csName[str(Axis_k)] = KLDecomposition_csName_Axis_k.listOfBasicFctsUsingLagrangeInterpolation_Axis_k
            
            listOfKeptEigenValues_csName[str(Axis_k)] = KLDecomposition_csName_Axis_k.finalEigVals_Axis_k
            
            listOfNbOfBasisFcts_csName.append(KLDecomposition_csName_Axis_k.nbOfFcts_Axis_k)
            
        listOfBasisFcts[csName] =  listOfBasisFcts_csName
        listOfNbOfBasisFcts[csName] = listOfNbOfBasisFcts_csName
        listOfKeptEigenValues[csName] = listOfKeptEigenValues_csName
        
        CardinalOfBasisFcts = 1
        for Axis_k in range(dimension):
            CardinalOfBasisFcts = CardinalOfBasisFcts*listOfNbOfBasisFcts[csName][Axis_k]  
   
    return listOfBasisFcts,  listOfNbOfBasisFcts, listOfKeptEigenValues
    
###===================================================    
def getListOfBasisFctsByAxis(dimension, listOfCrossSectionNamesWithArgs, listOfBasisFcts,  listOfNbOfBasisFcts):
   
    listOfMaxNbFcts = {}
    listOfAllBasisFcts = {}
    
    for Axis_k in range(dimension):  
        maxOfNbFct_Axis_k = 0
        listOfAllBasisFcts_Axis_k = []
        
        for csName in listOfCrossSectionNamesWithArgs:
            if maxOfNbFct_Axis_k < listOfNbOfBasisFcts[csName][Axis_k]:
                maxOfNbFct_Axis_k = listOfNbOfBasisFcts[csName][Axis_k]
            ## gather all functions for a direction concerned
            listOfAllBasisFcts_Axis_k =  listOfAllBasisFcts_Axis_k + listOfBasisFcts[csName][str(Axis_k)] 

        
        listOfMaxNbFcts[str(Axis_k)] = maxOfNbFct_Axis_k
        listOfAllBasisFcts[str(Axis_k)] = listOfAllBasisFcts_Axis_k
  
    return listOfMaxNbFcts, listOfAllBasisFcts
    
###===================================================   
def getListOfPointsForGreedyAlgo(listOfDomainBorders, listOfNbPointsForGreedyAlgo):
    listOfPointsForGreedyAlgo = {}    
   
    dimension = len(listOfDomainBorders)
    for Axis_k in range(dimension):
        nbExtremes = len(listOfDomainBorders[str(Axis_k)])

        listOfPointsForGreedyAlgo_Axis_k = np.asarray([])
        
        for indexOfBorder in range(nbExtremes-1):
            leftBorder = listOfDomainBorders[str(Axis_k)][indexOfBorder]
            rightBorder = listOfDomainBorders[str(Axis_k)][indexOfBorder + 1]
            listOfPointsForGreedyAlgo_Axis_k_OnSubInterval =  np.linspace(leftBorder,rightBorder,num=listOfNbPointsForGreedyAlgo[Axis_k])
            if indexOfBorder > 0 :
                listOfPointsForGreedyAlgo_Axis_k_OnSubInterval = np.delete(listOfPointsForGreedyAlgo_Axis_k_OnSubInterval, 0)
            listOfPointsForGreedyAlgo_Axis_k = np.concatenate((listOfPointsForGreedyAlgo_Axis_k, listOfPointsForGreedyAlgo_Axis_k_OnSubInterval), axis=0)
        
        listOfPointsForGreedyAlgo[str(Axis_k)] =  listOfPointsForGreedyAlgo_Axis_k
        
    print "========================="      
    print "listOfPointsForGreedyAlgo", listOfPointsForGreedyAlgo   
    print "========================="  
    #raw_input()
    return listOfPointsForGreedyAlgo

###===================================================  
### Construct a list of points/Axis_k, used for all basis functions in this Axis_k 
###=================================================== 
def getFinalEmpiricalInterpolationPoints(dimension, listOfMaxNbFcts, listOfAllBasisFcts, listOfPointsForGreedyAlgo):
    listOfBasisFctOnPointsForGreedyAlgo = {}
    #listOfPointsForGreedyAlgo = listOfTuckerGridNodes
    
    for Axis_k in range(dimension):
        listOfBasisFctOnPointsForGreedyAlgo_k = []
        for fct in listOfAllBasisFcts[str(Axis_k)]:
           
            listOfBasisFctOnPointsForGreedyAlgo_k.append(np.asarray(LI.getInterpolationArr(fct, listOfPointsForGreedyAlgo[str(Axis_k)])))
    
        listOfBasisFctOnPointsForGreedyAlgo[str(Axis_k)] = listOfBasisFctOnPointsForGreedyAlgo_k   
   
   
    ###========== Finding the empirical point for each direction ==== 
    listOfFinalEmpiricalInterpolationPoints = {}    
    
    for Axis_k in range(dimension):       
        listOfFinalEmpiricalInterpolationPoints_Axis_k = greedy.getEmpiricalPoints(listOfBasisFctOnPointsForGreedyAlgo[str(Axis_k)], \
                  listOfPointsForGreedyAlgo[str(Axis_k)], listOfMaxNbFcts[str(Axis_k)])[0]  ### [1] is for indexes of empirical points
                  
        listOfFinalEmpiricalInterpolationPoints[str(Axis_k)] = listOfFinalEmpiricalInterpolationPoints_Axis_k
   
    print "========================="  
    print "listOfFinalEmpiricalInterpolationPoints", listOfFinalEmpiricalInterpolationPoints
    print "========================="
    #raw_input()
    
    return listOfFinalEmpiricalInterpolationPoints
    
###===================================================  
def verifyInterpolPointsInTuckerNodesAndReferenceGrid(dimension, interpolationPoints, TuckerGridNodes, referencePoints):    
    
    for Axis_k in range(dimension):
        interpolPoints = interpolationPoints[str(Axis_k)]
        TuckerNodes = TuckerGridNodes[str(Axis_k)]
        refPoints = referencePoints[str(Axis_k)]
        
        for point in interpolPoints:
            if point not in TuckerNodes:
                print "point =", point
                warnings.warn("error",'empirical point is not in Tucker nodes')
                
            ### Find the point in refPoints which is the nearest point to point 
            j = min(range(len(refPoints)), key=lambda i: abs(refPoints[i]- point))
            refPoint = refPoints[j]
            if point == 0:
                eps = np.power(10.0,-10)
                if abs(refPoint - point) > eps:
                    print "point =", point
                    warnings.warn("error",'empirical point is not in reference points')
            else:
                eps = np.power(10.0,-4)
                if abs((refPoint - point)/point) > eps:
                    print "point =", point
                    warnings.warn("error",'empirical point is not in reference points')
                    
    for Axis_k in range(dimension):
        interpolPoints = interpolationPoints[str(Axis_k)]
        TuckerNodes = TuckerGridNodes[str(Axis_k)]
        refPoints = referencePoints[str(Axis_k)]
        print "===================="
        print "Axis_k = ", Axis_k
        print "interpolPoints =", interpolPoints
        print "refPoints =", refPoints
        print "===================="
  
###===================================================      
### From determined greedy points for each direction, each cross-section chooses from these points for its own points     
###===================================================     
def getListOfInterpolationPointsForEachSection(dimension, listOfCrossSectionNamesWithArgs, listOfBasisFcts, listOfNbOfBasisFcts, listOfFinalEmpiricalInterpolationPoints):
    listOfBasisFctOnFinalPointsForGreedyAlgo = {}
    
    for csName in listOfCrossSectionNamesWithArgs:  
        listOfBasisFctOnFinalPointsForGreedyAlgo_csName = {}
        for Axis_k in range(dimension):
            listOfBasisFctOnFinalPointsForGreedyAlgo_Axis_k = []
            points = listOfFinalEmpiricalInterpolationPoints[str(Axis_k)]
            for fct in listOfBasisFcts[csName][str(Axis_k)]:
                listOfBasisFctOnFinalPointsForGreedyAlgo_Axis_k.append(np.asarray(LI.getInterpolationArr(fct, points)))
                
            listOfBasisFctOnFinalPointsForGreedyAlgo_csName[str(Axis_k)] = listOfBasisFctOnFinalPointsForGreedyAlgo_Axis_k
        listOfBasisFctOnFinalPointsForGreedyAlgo[csName] = listOfBasisFctOnFinalPointsForGreedyAlgo_csName

    
    listOfInterpolationPoints = {}
    listOfPointsForGreedyAlgo =  listOfFinalEmpiricalInterpolationPoints
    for csName in listOfCrossSectionNamesWithArgs:  
        listOfInterpolationPoints_csName = {}
        for Axis_k in range(dimension):
            listOfInterpolationPoints_csName[str(Axis_k)] = greedy.getEmpiricalPoints(listOfBasisFctOnFinalPointsForGreedyAlgo[csName][str(Axis_k)], \
                  listOfPointsForGreedyAlgo[str(Axis_k)], listOfNbOfBasisFcts[csName][Axis_k])[0]  ### [1] is for indexes of empirical points
        listOfInterpolationPoints[csName] =  listOfInterpolationPoints_csName

          
    return listOfInterpolationPoints
 

    
