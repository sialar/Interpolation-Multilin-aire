# -*- coding: utf-8 -*-
import os
import dkmodel
import dktools
import dkbase
import array
import datetime
import numpy as np
import sys
from numpy import linalg as LA
import scipy
import random
import crossSections
import NeutronicLibraryManager
import Symbols
import geometricElements
from   Geometry import createFlatRegionName
import geometricElements.datarodbanks
import Geometry
loadFromDKLib = NeutronicLibraryManager.loadFromDKLib
###===================================================
### Nouveaux programmes:
import multiDimensIntegralMethods as M
import KLdecomposition as KL
import discretizationChoice as discret
import Tuckerdecomposition as Tucker
import LagrangeInterpolation as LI
import allCrossSections


    
##===================================================     
def getInforOfDomainForOneTuckerGrid(crossSection):
    """
  Argument :
    crossSection[une chaine de caractère] : définier le nom d'une sections
  Retour :
    une liste contient :
      + listOfNbPointsForOneTuckerGrid [une liste de "d" entiers]: le nombre de points lu axe par axe à partir un dklib
      + listOfNodesForOneTuckerGrid [une liste de "d" listes des flottants]: la discrétisation sur tous les axes lue à partir un dklib 
      
  ATTENTION : il y a des valeurs nominales dans quelques directions    
    """
    feedbackMesh = crossSection.getMesh()
    dimension  = feedbackMesh.nAxis()
    
    ### A Tucker grid have ~~ 10 points for the pricipal direction and 2 points for the others
    
    listOfNbPointsForOneTuckerGrid = []
    listOfNodesForOneTuckerGrid = []    
    
    
    for Axis_k in range(dimension):
            nodes_Axis_k = feedbackMesh.getRectilinearAxisSteps(Axis_k)           

            listOfNbPointsForOneTuckerGrid.append(len(nodes_Axis_k))
            
            listOfNodesForOneTuckerGrid.append(nodes_Axis_k)
            
    return [listOfNbPointsForOneTuckerGrid,  listOfNodesForOneTuckerGrid]    

 ##===================================================
def getInforOfDomain(listOfDKlibsForKLDecomposition):
    ###
    ### !!!! ATTENTION !!!!!! Fonction absolument pas générique. Il faut l'adapter en focntion du découpage des axes et dsi c'est le domaine étendu ou pas
    ### Le problèeme est de savoir si on garde l point nominal ou pas pour chaque axe (sauf burnup)
    ###
    """
    Argument : 
      - listOfDKlibsForKLDecomposition[une liste de "d" dklibs] : correspond à "d" grille de Tucker_post
    Retour :
      - listOfDomainBorders [un dictionnaire]:
                + clé [une chaine de caractères] : désigne l'axe. ex :'0', '1'
                + valeur [liste de flottants] :  une liste simple de bornes [min, max]
      
      - listOfNbPointsOnEachIntervalForFinalDiscretization[dictionnaire]:
                 + clé [une chaine de caractères] : désigne l'axe. ex :'0', '1'
                 + valeur [liste de listes d'entiers]: ce sont les nombres totals de points dans la discrétisation
                  '1': [2, [5], 3, 3, 3], '0': [[19], 3, 3, 3, 3], '3': [2, 3, 3, [7], 3]
                  (- sans prendre en compte le coupage
                  - prendre en compte si le point nominal est inclus dedans)
                  
      - listOfTuckerGridNodesFromDKlibsn[dictionnaire]:
                 + clé [une chaine de caractères] : désigne l'axe. ex :'0', '1'
                 + valeur [liste des sets]: valeurs de la discrétisation lue, axe par axe, à partir un dklib
                 e.x : '1': [(0.0, 80000.0),(20.0, 280.6749572753906, 550.0, 910.0, 1539.324951171875, 1800.0), ...]
                 
                  ATTENTION : il y a des valeurs nominales dans quelques directions 
    """

    crossSectionName = "macro_totale"
    dimension  = len(listOfDKlibsForKLDecomposition) #feedbackMesh.nAxis()
    
    listOfDomainBorders = {}
    listOfNbPointsOnEachIntervalForFinalDiscretization = {}
    listOfTuckerGridNodesFromDKlibs  = {}
    
    for Axis_k in range(dimension):
      
        crossSectionTotal =  crossSections.getCrossSections(listOfDKlibsForKLDecomposition[Axis_k], crossSectionName)
        crossSection =  crossSectionTotal[0]
        listOfNbPointsForOneTuckerGrid,  listOfNodesForOneTuckerGrid =\
                                                          getInforOfDomainForOneTuckerGrid(crossSection)
                                                         
    
        
        listOfDomainBorders[str(Axis_k)] = [listOfNodesForOneTuckerGrid[Axis_k][0], listOfNodesForOneTuckerGrid[Axis_k][-1]]
        listOfNbPointsOnEachIntervalForFinalDiscretization[str(Axis_k)] =  listOfNbPointsForOneTuckerGrid
        listOfTuckerGridNodesFromDKlibs[str(Axis_k)] = listOfNodesForOneTuckerGrid
        
        nbOfPointsWithNominalValues =  listOfNbPointsOnEachIntervalForFinalDiscretization[str(Axis_k)][Axis_k]
        
        ###  Mettre en commentaire Axis_k si la subdivision est réaliée en point nominal
        if dimension == 5: 
            # cet axe ext l'axe temperature combustibe, on retire le point nominal en Tc
            if Axis_k == 1:
                listOfNbPointsOnEachIntervalForFinalDiscretization[str(Axis_k)][Axis_k] =  nbOfPointsWithNominalValues -1 ### Ne compter pas les points nominals dans la discretization   
                
            # cet axe ext l'axe densité modératuer, et dans le cas (qui est ici) on est dans le domaine étendu et on garde le point nominal
            # car le découpage de l'axe a été fait pour conserver le point nominal. Don con commente
            #if Axis_k == 2:
                #listOfNbPointsOnEachIntervalForFinalDiscretization[str(Axis_k)][Axis_k] =  nbOfPointsWithNominalValues -1 ### Ne compter pas les points nominals dans la discretization
                
            if Axis_k == 3:
                listOfNbPointsOnEachIntervalForFinalDiscretization[str(Axis_k)][Axis_k] =  nbOfPointsWithNominalValues -1 ### Ne compter pas les points nominals dans la discretization
             
            # cet axe qui est l'axe en Xenon, n'est pas coupé en 2 mais dans le cas du domaine étendu [0, 2] la valeur nominale =1, est aussi à cause de la symmetrie, un point de CC
            # donc il ne faut pas le retirer
            #if Axis_k == 4:
                #listOfNbPointsOnEachIntervalForFinalDiscretization[str(Axis_k)][Axis_k] =  nbOfPointsWithNominalValues -1 ### Ne compter pas les points nominals dans la discretization   
     
        listOfNbPointsOnEachIntervalForFinalDiscretization[str(Axis_k)][Axis_k] = [listOfNbPointsOnEachIntervalForFinalDiscretization[str(Axis_k)][Axis_k]] 
    return [listOfDomainBorders, listOfNbPointsOnEachIntervalForFinalDiscretization, listOfTuckerGridNodesFromDKlibs]       
      
      ##=====================================  
  ##===================================================     
def getAllListOfFinalCrossSection(listOfCrossSectionNames,listOfDKlibsForKLDecomposition):
    """
    Fonction qui détermine toutes les sections pour tous types de dklib : multilinéarie, Tucker, référence...
    Argument : 
      - listOfCrossSectionNames[liste de chaine de caractère] : e.x [ "macro_totale", "macro_absorption",  "macro_scattering", "macro_nu*fission",  "macro_fission"]  
      - dklibName_getCrossSectionsWithArgs, ..., listOfDKlibsForKLDecomposition: dklib ou des dklibs sortis de GAB
    Retours:
       - listOfCrossSectionNamesWithArgs [liste de chaine de caractères] : précise les sections avec groupe d'énergy, anisotropy ...
        e.x listOfCrossSectionNamesWithArgs = 
        ['macro_totale0', 'macro_totale1', 'macro_scattering000', 'macro_scattering001', 'macro_scattering010', 'macro_scattering011']
      
      
    """

    
    listOfCrossSectionNamesWithArgs = []
    listOfCrossSectionsForTuckerGridsFromDklibs = {}
    
    dklibName_getCrossSectionsWithArgs = listOfDKlibsForKLDecomposition[0]
    for crossSectionName in listOfCrossSectionNames :   
        listOfArgs = crossSections.getListOfArgs(dklibName_getCrossSectionsWithArgs, crossSectionName)
        for arg in listOfArgs:
            csNameWithArgument = crossSections.getCrossSectionNameWithArg(crossSectionName, arg)
            listOfCrossSectionNamesWithArgs.append(csNameWithArgument)

            listOfTuckerDataFromDklibsForOneCrossSection = {}
            for Axis_k in range(len(listOfDKlibsForKLDecomposition)):
               
                listOfTuckerDataFromDklibsForOneCrossSection[str(Axis_k)] = crossSections.getfinalCrossSection(listOfDKlibsForKLDecomposition[Axis_k], crossSectionName, arg)
                
            listOfCrossSectionsForTuckerGridsFromDklibs[csNameWithArgument] = listOfTuckerDataFromDklibsForOneCrossSection
    return listOfCrossSectionNamesWithArgs, listOfCrossSectionsForTuckerGridsFromDklibs 
           
          
   
##===================================================  
if __name__ == "__main__":
  
    # ne prend rien en entrée, mais il faut modifier les fonction sau-dessus notamment les methodes style etListOfDklibs_MOX qui indique les DKLibs
    # en fait on n'a besoin que d'une seule fonction générique getListOfDklibs mais pour la restitution Hieu a indiqué toutes les bibliothèques utilisées
    
    listOfNominalValues = {'fuel_temperature' : 550.0, 'coolant_density': 0.71692097187, 'b10_concentration': 3.97352914661e-06, 'xenon_level': 1, 'coolant_temperature': 304.600006104}
    listOfAxis_Parameter_5dimension = {'0': 'BURNUP', '1': 'fuel_temperature', '2': 'coolant_density', '3': 'b10_concentration', '4': 'xenon_level'}
    ###listOfAxis_Parameter_6dimension = {'0': 'BURNUP', '1': 'fuel_temperature', '2':'coolant_temperature', '3': 'coolant_density', '4': 'b10_concentration', '5': 'xenon_level'}
 
    ##===================================================
    ### ATTENTION : A MODIFIER SUIVNAT LE CAS : UOX ou UOX-Gd ou MOX, etc.
    ##===================================================     
    
    dklibName0 = "/local00/home/H41070/EspaceTravail/TestCOCAGNE/DKLibs_Gauss_TH/DKLIB_5parametres_extendDomaine_ByDomain_GAB231/burnup_3_9_9_pointsCC_gado_bis"    
    dklibName1 = "/local00/home/H41070/EspaceTravail/TestCOCAGNE/DKLibs_Gauss_TH/DKLIB_5parametres_extendDomaine_ByDomain_GAB231/Tc_5_gado"    
    dklibName2 = "/local00/home/H41070/EspaceTravail/TestCOCAGNE/DKLibs_Gauss_TH/DKLIB_5parametres_extendDomaine_ByDomain_GAB231/Rho_5_5_gado"    
    dklibName3 = "/local00/home/H41070/EspaceTravail/TestCOCAGNE/DKLibs_Gauss_TH/DKLIB_5parametres_extendDomaine_ByDomain_GAB231/B10_7_gado"  
    dklibName4 = "/local00/home/H41070/EspaceTravail/TestCOCAGNE/DKLibs_Gauss_TH/DKLIB_5parametres_extendDomaine_ByDomain_GAB231/Xe_5_gado"
    
    listOfDKlibsForKLDecomposition = [dklibName0, dklibName1, dklibName2, dklibName3, dklibName4]
    
    ##===================================================
    dimension = len(listOfDKlibsForKLDecomposition)

    ##===================================================
    ## Les types de sections dans le test
    ##===================================================
   
    listOfCrossSectionNames = [ "macro_totale", "macro_absorption",  "macro_scattering", "macro_nu*fission",  "macro_fission"]  

    ##===================================================
    ## On récupère les sections efficaces avec ses arguments (goupe d'énergy, anisotropy, etc) pour
    ## + la liste des noms précis de ces sections (avec le numero du groupe et l'anisotropie): listOfCrossSectionNamesWithArgs
    ## + la décomposition de Tucker : listOfCrossSectionsForTuckerGridsFromDklibs
    ##===================================================
    listOfCrossSectionNamesWithArgs,listOfCrossSectionsForTuckerGridsFromDklibs   = getAllListOfFinalCrossSection(listOfCrossSectionNames,  listOfDKlibsForKLDecomposition)

    listOfDomainBorders, listOfNbPointsOnEachIntervalForFinalDiscretization, listOfTuckerGridNodesFromDKlibs\
    = getInforOfDomain(listOfDKlibsForKLDecomposition)
    

    ###===================================================
    ### CHANGE HERE FOR EACH CHOICE OF LIST OF DKLIBS 
    ###===================================================
    ### Définir quel axe est en coupage, les extremes de segments de coupage  aisi que le nombre de points dans (k-1) premières intervalles (sur k intervalles)  
    ###===================================================
    listOfSubdivisionDirection = {}
    listOfValuesForSubdivision = {}
    listOfNumberOfPointsForSubdivision = {}
    
    # on ne remplit les 3 tableaux précédents que pour les axes qui sont découpés
    
    listOfSubdivisionDirection['0'] = 1
    listOfValuesForSubdivision['0'] = [150, 10000]
    listOfNumberOfPointsForSubdivision['0'] = [3, 9]  # 3 points entre [debut, 150], 9 points entre [150, 10000], le reste pour [10000, fin]. A priori c'est inutile et redondant

    
    listOfSubdivisionDirection['2'] = 1
    listOfValuesForSubdivision['2'] = [0.71692097187]
    listOfNumberOfPointsForSubdivision['2'] = [5]
    

    ###===================================================
    ### Pré-définir les discrétisations pour "d" (d = dimension) grilles de Tucker
    ###===================================================
    D = discret.discretizationChoice(listOfDomainBorders, listOfNbPointsOnEachIntervalForFinalDiscretization,\
                    listOfSubdivisionDirection, listOfValuesForSubdivision, listOfNumberOfPointsForSubdivision)

    D.getDiscretization() ### Update the listOfDomainBorders with the listOfSubdivisionDirection     
    listOfIntegralMethods = D.listOfIntegralMethods
   
    ###===================================================  
    ### Toutes les dicrétisations fines sur "d" axes sont calculées et mettre dans "listOfTuckerGridNodes"
    ###=================================================== 
    listOfTuckerGridNodes = {}
    print listOfCrossSectionNamesWithArgs
    csName = listOfCrossSectionNamesWithArgs[0] #'macro_totale0'
   
    listOfKLDecompositionToDefineListOfTuckerGridNodes = {}
    for Axis_k in range(len(listOfDKlibsForKLDecomposition)):
        listOfKLDecompositionToDefineListOfTuckerGridNodes[str(Axis_k)] = KL.KL(Axis_k, listOfDomainBorders, listOfNbPointsOnEachIntervalForFinalDiscretization[str(Axis_k)],\
        listOfTuckerGridNodesFromDKlibs[str(Axis_k)], listOfIntegralMethods[str(Axis_k)], listOfCrossSectionsForTuckerGridsFromDklibs[csName][str(Axis_k)], csName)        
        listOfTuckerGridNodes[str(Axis_k)] = listOfKLDecompositionToDefineListOfTuckerGridNodes[str(Axis_k)].listOfTuckerGridNodes_Axis_k[Axis_k]
        
    #print "  listOfDomainBorders",    listOfDomainBorders
    #print "listOfNbPointsOnEachIntervalForFinalDiscretization", listOfNbPointsOnEachIntervalForFinalDiscretization
    #print "listOfTuckerGridNodesFromDKlibs ", listOfTuckerGridNodesFromDKlibs
    #print "listOfIntegralMethods ", listOfIntegralMethods
    #print "listOfCrossSectionsForTuckerGridsFromDklibs", listOfCrossSectionsForTuckerGridsFromDklibs
    #print "listOfTuckerGridNodes ", listOfTuckerGridNodes
    #raw_input()
    

    ###===================================================
    ### Create a list of KLdecomposition to apply the greedy algo for all basis fcts
    ###===================================================
    listOfBasisFcts, listOfNbOfBasisFcts = allCrossSections.getAllBasisFctsForAllSections(listOfCrossSectionNamesWithArgs, listOfDomainBorders, listOfNbPointsOnEachIntervalForFinalDiscretization,\
            listOfTuckerGridNodesFromDKlibs, listOfIntegralMethods, listOfCrossSectionsForTuckerGridsFromDklibs) 
    
    ###============WITH GREEDY============================
    ### Collecter touts les fonctions de base, axe par axe pour appliquer greedy ensuit
    ###===================================================
    listOfMaxNbFcts, listOfAllBasisFcts = allCrossSections.getListOfBasisFctsByAxis(len(listOfDKlibsForKLDecomposition),\
                          listOfCrossSectionNamesWithArgs, listOfBasisFcts,  listOfNbOfBasisFcts)

    ###===================================================  
    ### Ici, on fait le greedy sur les points CC dans la discrétisation fine 
    ### (on peut proposer un autre échantillon)
    ###===================================================
    listOfPointsForGreedyAlgo  =  listOfTuckerGridNodes 
  
    ###==============FinalInterpolationPoints for all sections =====================================  
    listOfFinalEmpiricalInterpolationPoints = allCrossSections.getFinalEmpiricalInterpolationPoints(len(listOfDKlibsForKLDecomposition),\
                                listOfMaxNbFcts, listOfAllBasisFcts, listOfPointsForGreedyAlgo)
                                
    ###=============Verify the empirical points with Tucker nodes and reference points======================================    
    #allCrossSections.verifyInterpolPointsInTuckerNodesAndReferenceGrid(dimension, listOfFinalEmpiricalInterpolationPoints, listOfTuckerGridNodes, referencePoints)   
    print "========================="
    print "listOfFinalEmpiricalInterpolationPoints", listOfFinalEmpiricalInterpolationPoints
    print "========================="