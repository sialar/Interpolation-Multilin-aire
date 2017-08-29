## -*- coding: utf-8 -*-

import multiDimensIntegralMethods as M
import KLdecomposition as KL
import discretizationChoice as discret
import Tuckerdecomposition as Tucker
import LagrangeInterpolation as LI
import fctMultiVariates as fct
import allCrossSections
 ##===================================================
import array
import datetime
import numpy as np
import sys
import random

def diff(nom, a, b):
    d = ( a/b - 1. ) * 1.E5
    d2 = np.power((a-b)*1.E5,2)
    print 'ECART SECTION ',nom,a,b,d
    return abs(d), d2
###===================================================       

##===================================================     
def getNbPointsForOneTuckerGrid(dimension, axis_i):    
    ### A Tucker grid have ~~ 10 points for the pricipal direction and 2 points for others
    
    listOfNbPointsForOneTuckerGrid = []
    listOfNodesForOneTuckerGrid = []    
    
    n = 2
    N = 25
    for Axis_k in range(dimension):
            listOfNbPointsForOneTuckerGrid.append(n)             
            
    listOfNbPointsForOneTuckerGrid[axis_i] = N 
    ###======= Si x0 est divisé
    #if  axis_i == 0:
        #listOfNbPointsForOneTuckerGrid[axis_i] = 25
    
    ### Si x1 est divisé
    #if  axis_i == 1:
        #listOfNbPointsForOneTuckerGrid[axis_i] = 17 
    
    return listOfNbPointsForOneTuckerGrid  

### For all Tucker decompositions : $d$ dimensions, $d$ Tucker grid where each grid has many points for principal direction and 2 points for others
 ##===================================================
def getInforOfDomain(dimension):

    listOfDomainBorders = {}
    listOfNbPointsOnEachIntervalForFinalDiscretization = {}
    listOfTuckerGridNodes_predefined  = {}
    
    bords = [0.0,1.]
    for Axis_k in range(dimension):
        listOfDomainBorders[str(Axis_k)] = [bords[0], bords[-1]]
        
    for Axis_k in range(dimension):
        listOfNbPointsForOneTuckerGrid_k = getNbPointsForOneTuckerGrid(dimension, Axis_k)
        listOfNbPointsOnEachIntervalForFinalDiscretization[str(Axis_k)] =  listOfNbPointsForOneTuckerGrid_k

        ### peusedo comme dklib
        listOfNodesForOneTuckerGrid = []
        for Axis_j in range(dimension):
            ### Choice 1 : random
            ##===================================================
            #rand = random.uniform(listOfDomainBorders[str(Axis_k)][0],listOfDomainBorders[str(Axis_k)][-1])
            #mid = (listOfDomainBorders[str(Axis_k)][0] + listOfDomainBorders[str(Axis_k)][-1])/2
           
            #if rand >= mid:
                #listOfNodesForOneTuckerGrid.append([listOfDomainBorders[str(Axis_k)][0], rand])
            #elif rand < mid:
                #listOfNodesForOneTuckerGrid.append([rand, listOfDomainBorders[str(Axis_k)][-1]])
            ##===================================================
            ### Choice 2 : extremes
            listOfNodesForOneTuckerGrid.append([listOfDomainBorders[str(Axis_k)][0], listOfDomainBorders[str(Axis_k)][-1]])
            
        points_k = np.linspace ( listOfDomainBorders[str(Axis_k)][0],  listOfDomainBorders[str(Axis_k)][1], num=listOfNbPointsOnEachIntervalForFinalDiscretization[str(Axis_k)][Axis_k]) 
        listOfNodesForOneTuckerGrid[Axis_k] = points_k
        
        listOfTuckerGridNodes_predefined[str(Axis_k)] = listOfNodesForOneTuckerGrid
        
        
        ### ============== Put in list the element whose direction is treated by the couplage  ===========      
        listOfNbPointsOnEachIntervalForFinalDiscretization[str(Axis_k)][Axis_k] = [listOfNbPointsOnEachIntervalForFinalDiscretization[str(Axis_k)][Axis_k]] 
    return [listOfDomainBorders, listOfNbPointsOnEachIntervalForFinalDiscretization, listOfTuckerGridNodes_predefined]       
  

##===================================================  
def writeToFile(nameFile, data):
    File = open(nameFile,'w')
    for item in data:
        File.write("%s\n" % item)
    File.close()

##===================================================  
if __name__ == "__main__":
    #argument1 = sys.argv[1]
    ###================crossSections is replaced function ===========
    ### There is only one dklib containing only one cross-section = f
    ###==============================================================
    ###argsF = [1.0, 2.5, -3.0]
    argsF = [0.5, 0.25, 0.5] # Pour baisser les valeurs de f
    csName = "fct3d"
    dimension = len(argsF)
    f = fct.function(argsF, csName, dimension)
    
    ###===================================================
    listOfCrossSectionNamesWithArgs = [csName]
    listOfFinalCrossSectionOnFullGrid = {}
    listOfFinalCrossSectionOnFullGrid_reference = {}    
    listOfRefCrossSectionsToDetermineCoeffs = {}
    listOfCrossSectionsForTuckerGridsFromDklibs = {}
    
    
    listOfFinalCrossSectionOnFullGrid[csName] = f
    listOfFinalCrossSectionOnFullGrid_reference[csName] = f  
    listOfRefCrossSectionsToDetermineCoeffs[csName] = f
    
    listOfCrossSectionsForTuckerGridsFromDklibs_k = {}
    for Axis_k in range(dimension):        
        listOfCrossSectionsForTuckerGridsFromDklibs_k[str(Axis_k)] = f
        
    listOfCrossSectionsForTuckerGridsFromDklibs[csName] =  listOfCrossSectionsForTuckerGridsFromDklibs_k 
    ###===================================================
    
    ##===================================================

    listOfDomainBorders, listOfNbPointsOnEachIntervalForFinalDiscretization, listOfTuckerGridNodes_predefined\
    = getInforOfDomain(dimension)
    
    ##===================================================   
    for Axis_k in range(dimension):
        print "====================="
        print "listOfTuckerGridNodes_predefined[str(Axis_k)]",listOfTuckerGridNodes_predefined[str(Axis_k)]
        print "====================="
    ###===================================================
    ### CHANGE HERE FOR EACH CHOICE OF LIST OF DKLIBS 
    ###===================================================
    listOfCouplageDirection = {}
    listOfValuesForCouplage = {}
    listOfNumberOfPointsForCouplage = {}
    
    ###===================================================
    ### CHANGE HERE FOR EACH DIRECTION, E.X :
    ###===================================================
    #listOfCouplageDirection['0'] = 1   
    #listOfValuesForCouplage['0'] = [0.25, 0.5] 
    #listOfNumberOfPointsForCouplage['0'] = [9,9]
    ###===================================================
    
    ###===================================================
    D = discret.discretizationChoice(listOfDomainBorders, listOfNbPointsOnEachIntervalForFinalDiscretization,\
                    listOfCouplageDirection, listOfValuesForCouplage, listOfNumberOfPointsForCouplage)
    D.getDiscretization() ### Update the listOfDomainBorders with the listOfCouplageDirection     
    listOfMethodsIntegrals = D.listOfMethodsIntegrals
   
    #print "====================="
    #print "listOfMethodsIntegrals", listOfMethodsIntegrals  
    #print "====================="  
    #raw_input()
    

    listOfTuckerGridNodes = {}
    #print listOfCrossSectionNamesWithArgs
    csName = listOfCrossSectionNamesWithArgs[0] 
   
    listOfKLDecompositionToDefineListOfTuckerGridNodes = {}
    for Axis_k in range(dimension):
        #print "=====================" 
        #print " Axis_k =",  Axis_k
        #print "====================="
        #print "listOfDomainBorders =", listOfDomainBorders
        #print "====================="
        #print "listOfNbPointsOnEachIntervalForFinalDiscretization[str(Axis_k)]", listOfNbPointsOnEachIntervalForFinalDiscretization[str(Axis_k)]
        #print "====================="
        #print "listOfTuckerGridNodes_predefined[str(Axis_k)]",listOfTuckerGridNodes_predefined[str(Axis_k)]
        #print "====================="
        #print " listOfMethodsIntegrals[str(Axis_k)]",  listOfMethodsIntegrals[str(Axis_k)]
        #print "====================="
        #print "listOfCrossSectionsForTuckerGridsFromDklibs[csName][str(Axis_k)]",listOfCrossSectionsForTuckerGridsFromDklibs[csName][str(Axis_k)]
        #print "====================="
        #raw_input()
        
        listOfKLDecompositionToDefineListOfTuckerGridNodes[str(Axis_k)] = KL.KL(Axis_k, listOfDomainBorders, listOfNbPointsOnEachIntervalForFinalDiscretization[str(Axis_k)],\
        listOfTuckerGridNodes_predefined[str(Axis_k)], listOfMethodsIntegrals[str(Axis_k)], listOfCrossSectionsForTuckerGridsFromDklibs[csName][str(Axis_k)], csName)
        ### Pour appeler la fct qui fournit toutes les fcts de base (orthonormée) dans la direction k
        listOfTuckerGridNodes[str(Axis_k)] = listOfKLDecompositionToDefineListOfTuckerGridNodes[str(Axis_k)].listOfTuckerGridNodes_Axis_k[Axis_k]
    
  
    print "====================="
    print "listOfTuckerGridNodes", listOfTuckerGridNodes
    print "====================="
    #raw_input()
    ###===================================================
    ### Create a list of KLdecomposition to applie the greedy algo for all basis fct
    listOfBasisFcts, listOfNbOfBasisFcts, listOfKeptEigenValues = allCrossSections.getAllBasisFctsForAllSections(listOfCrossSectionNamesWithArgs, listOfDomainBorders, listOfNbPointsOnEachIntervalForFinalDiscretization,\
            listOfTuckerGridNodes_predefined, listOfMethodsIntegrals, listOfCrossSectionsForTuckerGridsFromDklibs) 
    
    ###============WITH GREEDY============================
    ###===================================================
    listOfMaxNbFcts, listOfAllBasisFcts = allCrossSections.getListOfBasisFctsByAxis(dimension,\
                          listOfCrossSectionNamesWithArgs, listOfBasisFcts,  listOfNbOfBasisFcts)

    ###===================================================   
   
    ###==============Creating the basis fucntions on the finest sample of points=====================================

    listOfPointsForGreedyAlgo  =  listOfTuckerGridNodes 
  
    ###==============FinalInterpolationPoints for all sections =====================================  
    listOfFinalEmpiricalInterpolationPoints = allCrossSections.getFinalEmpiricalInterpolationPoints(dimension,\
                                listOfMaxNbFcts, listOfAllBasisFcts, listOfPointsForGreedyAlgo)                              
 
    
    ###===========Extract fromm FinalInterpolationPoints the interpolation point for each section========================================
    listOfInterpolationPoints = allCrossSections.getListOfInterpolationPointsForEachSection(dimension,\
        listOfCrossSectionNamesWithArgs, listOfBasisFcts, listOfNbOfBasisFcts, listOfFinalEmpiricalInterpolationPoints)
        
    ###============END WITH GREEDY============================
    ###===================================================
    listOfTuckerDecomp = {}
    for csName in listOfCrossSectionNamesWithArgs:       
        listOfTuckerDecomp[csName] = Tucker.Tucker(listOfBasisFcts[csName], listOfNbOfBasisFcts[csName], listOfRefCrossSectionsToDetermineCoeffs[csName],listOfInterpolationPoints[csName], csName, listOfKeptEigenValues[csName])

        listOfTuckerDecomp[csName].getCoefficients()
    ###===================================================
    ref_values = []
    for Axis_k in range(dimension):
        nbP = 20
        values = np.linspace(listOfDomainBorders[str(Axis_k)][0], listOfDomainBorders[str(Axis_k)][-1],nbP)
        ref_values.append(values)
    ###===================================================
    ###=====================================================================
    nameOfTest = "f3d_analytique"
    ###=====================================================================
    ###=====================================================================
    listOfValues_Tucker = {}
    listOfValues_AP2= {}
    ####===================================================

    
    for csName in listOfCrossSectionNamesWithArgs: 
        count = 1
        values_Tucker = []
        values_AP2 = []
        for x1 in ref_values[0]:
            for x2 in ref_values[1]:
                for x3 in ref_values[2]:
                   
                    point = (x1,x2,x3) 
                    value_T = listOfTuckerDecomp[csName].evaluate(point)
           
                    value_A = listOfFinalCrossSectionOnFullGrid_reference[csName].evaluate(point)

                    values_Tucker.append(value_T)
  
                    values_AP2.append(value_A)
                    ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    count = count + 1
                            
        listOfValues_Tucker[csName] = np.asarray(values_Tucker)        
        listOfValues_AP2[csName]= np.asarray(values_AP2)                    
    
    ###======================Kinf============================= 
    
    listOfRelativeError_Tucker = {}
   
    listOfMaxRelativeError_Tucker = {}

     
    listOfMSE_Tucker = {}  


    ###======================Kinf============================= 
    for csName in listOfCrossSectionNamesWithArgs:
        nameFile_Err = "RelativeError_AP2_Tucker" + nameOfTest +"_" + csName       
        File_Err = open(nameFile_Err,'w')
         ### 4/3/2016: normaliser e_relative = |f-f~|/maxf
        maxAP2 = max(abs(listOfValues_AP2[csName]))
        listOfRelativeError_Tucker[csName] = 1E5*(listOfValues_Tucker[csName] - listOfValues_AP2[csName])/maxAP2
        
        
        ###==========Write the results of relative error of two methods in files===========
        count = 1
        for x1 in ref_values[0]:
            for x2 in ref_values[1]:
                for x3 in ref_values[2]:
                            value_A = listOfValues_AP2[csName][count-1]                            
                            value_T = listOfValues_Tucker[csName][count-1]
                            
                           
                            err_T = 1E5*(value_T - value_A)/maxAP2
                           
                            File_Err.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n"\
                            %(count,x1,x2,x3,value_A, value_T.real,  err_T.real))
                            count = count +1
        File_Err.close()
        ###==========Write the results of relative error of two methods in files===========
       
        listOfMaxRelativeError_Tucker[csName] = max(abs(listOfRelativeError_Tucker[csName]))     
      
        AP2_RS = np.sqrt(sum(np.power(listOfValues_AP2[csName],2)))
        
        listOfMSE_Tucker[csName]  = np.sqrt(sum(np.power(abs(listOfValues_Tucker[csName] - listOfValues_AP2[csName])*1E5, 2)))/AP2_RS   
       
        
        print "====================="
        print csName
        print "relativeError_Tucker (pcm)", listOfMaxRelativeError_Tucker[csName]  

        print "MSE_Tucker (pcm)", listOfMSE_Tucker[csName]
        print "====================="



