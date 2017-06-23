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

##=====================Exact==============================      
def computeKinf(listOfValues):
    """
  Fonction qui calcule le k_inf selon une formule analytique
  Argument : 
    "listOfValues" [un dictionnaire] :
      + clé [une chaine de caractères] : nom de section. ex : "macro_nu*fission0"
      + valeurs [liste de flottants] : les valeurs évaluées par Tucker ou multilinéaire ou GAB sur la grille de référence
  Retour :
    + une liste des valeurs de k_inf selon la formule analytique et sur la grille de référence 
    """
 
  
    nufi0 = listOfValues["macro_nu*fission0"]  
    nufi1 = listOfValues["macro_nu*fission1"]      
    tot0 = listOfValues["macro_totale0"]  
    tot1 = listOfValues["macro_totale1"]
    scatt000 = listOfValues["macro_scattering000"]
    scatt001 = listOfValues["macro_scattering001"]
    scatt010 = listOfValues["macro_scattering010"]
    scatt011 = listOfValues["macro_scattering011"]
    return ( nufi0 * (tot1 - scatt011) + nufi1*scatt010) / ( (tot0 - scatt000 ) * (tot1 - scatt011) - scatt001*scatt010) 
  
    
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
def getAllListOfFinalCrossSection(listOfCrossSectionNames, dklibName_fullGrid, dklibName_fullGrid_reference, dklibReference, listOfDKlibsForKLDecomposition):
    """
    Fonction qui détermine toutes les sections pour tous types de dklib : multilinéarie, Tucker, référence...
    Argument : 
      - listOfCrossSectionNames[liste de chaine de caractère] : e.x [ "macro_totale", "macro_absorption",  "macro_scattering", "macro_nu*fission",  "macro_fission"]  
      - dklibName_fullGrid, ..., listOfDKlibsForKLDecomposition: dklib ou des dklibs sortis de GAB
    Retours:
       - listOfCrossSectionNamesWithArgs [liste de chaine de caractères] : précise les sections avec groupe d'énergy, anisotropy ...
        e.x listOfCrossSectionNamesWithArgs = 
        ['macro_totale0', 'macro_totale1', 'macro_scattering000', 'macro_scattering001', 'macro_scattering010', 'macro_scattering011']
      
      - listOfFinalCrossSectionOnFullGrid [dictionnaire] (ici de dklib calculé sur la grille multilinéaire): 
        + clé [une chaine de caractères] : nom précis de section
        + valeur : section lue de dklib (champ aux noeuds)
        
       - listOfFinalCrossSectionOnFullGrid_reference [dictionnaire] (ici de dklib calculé sur la grille de référence pour mesurer les erreurs d'évaluations): 
        + clé [une chaine de caractères] : nom précis de section
        + valeur : section lue de dklib 
        
      - listOfRefCrossSectionsToDetermineCoeffs [dictionnaire]  (ici de dklib calculé sur la grille qui détermine "b" dans Aa=b pour résoudre les coefficients "a" dans la décomposition de Tucker): 
        + clé [une chaine de caractères] : nom précis de section
        + valeur : section lue de dklib 
      - listOfCrossSectionsForTuckerGridsFromDklibs   
      
    """

    
    listOfCrossSectionNamesWithArgs = []
    listOfFinalCrossSectionOnFullGrid = {}
    listOfFinalCrossSectionOnFullGrid_reference = {}    
    listOfRefCrossSectionsToDetermineCoeffs = {}
    listOfCrossSectionsForTuckerGridsFromDklibs = {}
    
    for crossSectionName in listOfCrossSectionNames :   
        listOfArgs = crossSections.getListOfArgs(dklibName_fullGrid, crossSectionName)
        for arg in listOfArgs:
            csNameWithArgument = crossSections.getCrossSectionNameWithArg(crossSectionName, arg)
            listOfCrossSectionNamesWithArgs.append(csNameWithArgument)
        
            listOfFinalCrossSectionOnFullGrid[csNameWithArgument] = crossSections.getfinalCrossSection(dklibName_fullGrid, crossSectionName, arg)
            listOfFinalCrossSectionOnFullGrid_reference[csNameWithArgument] = crossSections.getfinalCrossSection(dklibName_fullGrid_reference, crossSectionName, arg)
            listOfRefCrossSectionsToDetermineCoeffs[csNameWithArgument] = crossSections.getfinalCrossSection(dklibReference, crossSectionName, arg)
            
            listOfTuckerDataFromDklibsForOneCrossSection = {}
            for Axis_k in range(len(listOfDKlibsForKLDecomposition)):
               
                listOfTuckerDataFromDklibsForOneCrossSection[str(Axis_k)] = crossSections.getfinalCrossSection(listOfDKlibsForKLDecomposition[Axis_k], crossSectionName, arg)
                
            listOfCrossSectionsForTuckerGridsFromDklibs[csNameWithArgument] = listOfTuckerDataFromDklibsForOneCrossSection
    return listOfCrossSectionNamesWithArgs, listOfFinalCrossSectionOnFullGrid,  listOfFinalCrossSectionOnFullGrid_reference,\
           listOfRefCrossSectionsToDetermineCoeffs,listOfCrossSectionsForTuckerGridsFromDklibs 
           
          
###===================================================
def getListOfDklibs_UOX():
  ### dklibs pour UOX: grille de référence pour évaluer erreurs, grille pour multilinéaire, 5 grille de Tucker et une grille pour calculer les coefficients de Tucker

    dklibName_fullGrid_reference = "/local00/home/H41070/EspaceTravail/TestCOCAGNE/DKLibs_Gauss_TH/DKLIB_5parametres_extendDomaine_ByDomain_GAB231/reference_pleine"
    dklibName_fullGrid =  "/local00/home/H41070/EspaceTravail/TestCOCAGNE/DKLibs_Gauss_TH/DKLIB_5parametres_extendDomaine_ByDomain_GAB231/grille_pleine"

    dklibName0 = "/local00/home/H41070/EspaceTravail/TestCOCAGNE/DKLibs_Gauss_TH/DKLIB_5parametres_extendDomaine_ByDomain_GAB231/burnup_3_9_9_pointsCC"    
    dklibName1 = "/local00/home/H41070/EspaceTravail/TestCOCAGNE/DKLibs_Gauss_TH/DKLIB_5parametres_extendDomaine_ByDomain_GAB231/Tc_5"    
    dklibName2 = "/local00/home/H41070/EspaceTravail/TestCOCAGNE/DKLibs_Gauss_TH/DKLIB_5parametres_extendDomaine_ByDomain_GAB231/Rho_5_5"    
    dklibName3 = "/local00/home/H41070/EspaceTravail/TestCOCAGNE/DKLibs_Gauss_TH/DKLIB_5parametres_extendDomaine_ByDomain_GAB231/B10_7"  ###  
    dklibName4 = "/local00/home/H41070/EspaceTravail/TestCOCAGNE/DKLibs_Gauss_TH/DKLIB_5parametres_extendDomaine_ByDomain_GAB231/Xe_5"
      
    dklibReference = "/local00/home/H41070/EspaceTravail/TestCOCAGNE/DKLibs_Gauss_TH/DKLIB_5parametres_extendDomaine_ByDomain_GAB231/reference_CCpoints" #reference"### 

    listOfDKlibsForKLDecomposition = [dklibName0, dklibName1, dklibName2, dklibName3, dklibName4]
    
    return  dklibName_fullGrid, dklibName_fullGrid_reference, dklibReference, listOfDKlibsForKLDecomposition
##===================================================  

###============= UOX-Gd cas==========================
def getListOfDklibs_UOXGd():
  
     ### dklibs pour UOX-Gd : grille de référence pour évaluer erreurs, grille pour multilinéaire, 5 grille de Tucker et une grille pour calculer les coefficients de Tucker 
    dklibName_fullGrid_reference = "/local00/home/H41070/EspaceTravail/TestCOCAGNE/DKLibs_Gauss_TH/DKLIB_5parametres_extendDomaine_ByDomain_GAB231/reference_plein_gado"

    dklibName_fullGrid =  "/local00/home/H41070/EspaceTravail/TestCOCAGNE/DKLibs_Gauss_TH/DKLIB_5parametres_extendDomaine_ByDomain_GAB231/grille_pleine_gado" 
    dklibName0 = "/local00/home/H41070/EspaceTravail/TestCOCAGNE/DKLibs_Gauss_TH/DKLIB_5parametres_extendDomaine_ByDomain_GAB231/burnup_3_9_9_pointsCC_gado_bis"    
    dklibName1 = "/local00/home/H41070/EspaceTravail/TestCOCAGNE/DKLibs_Gauss_TH/DKLIB_5parametres_extendDomaine_ByDomain_GAB231/Tc_5_gado"    
    dklibName2 = "/local00/home/H41070/EspaceTravail/TestCOCAGNE/DKLibs_Gauss_TH/DKLIB_5parametres_extendDomaine_ByDomain_GAB231/Rho_5_5_gado"    
    dklibName3 = "/local00/home/H41070/EspaceTravail/TestCOCAGNE/DKLibs_Gauss_TH/DKLIB_5parametres_extendDomaine_ByDomain_GAB231/B10_7_gado"  
    dklibName4 = "/local00/home/H41070/EspaceTravail/TestCOCAGNE/DKLibs_Gauss_TH/DKLIB_5parametres_extendDomaine_ByDomain_GAB231/Xe_5_gado"

    dklibReference = "/local00/home/H41070/EspaceTravail/TestCOCAGNE/DKLibs_Gauss_TH/DKLIB_5parametres_extendDomaine_ByDomain_GAB231/reference_CCpoints_gado_bis" 

    listOfDKlibsForKLDecomposition = [dklibName0, dklibName1, dklibName2, dklibName3, dklibName4]
    
    return  dklibName_fullGrid, dklibName_fullGrid_reference, dklibReference, listOfDKlibsForKLDecomposition
##=================================================== 

###===================================================
###============= MOX  ================================
def getListOfDklibs_MOX():
    ### dklibs pour MOX: grille de référence pour évaluer erreurs, grille pour multilinéaire, 5 grille de Tucker et une grille pour calculer les coefficients de Tucker
   
    dklibName_fullGrid_reference = "/local00/home/H41070/EspaceTravail/TestCOCAGNE/DKLibs_Gauss_TH/DKLIB_EXTEND_MOX/reference_plein_mox"

    dklibName_fullGrid =  "/local00/home/H41070/EspaceTravail/TestCOCAGNE/DKLibs_Gauss_TH/DKLIB_EXTEND_MOX/grille_pleine_noRacineTc_mox"

    dklibName0 = "/local00/home/H41070/EspaceTravail/TestCOCAGNE/DKLibs_Gauss_TH/DKLIB_EXTEND_MOX/burnup_3_9_9_pointsCC_mox"    
    dklibName1 = "/local00/home/H41070/EspaceTravail/TestCOCAGNE/DKLibs_Gauss_TH/DKLIB_EXTEND_MOX/Tc_5p_mox"    
    dklibName2 = "/local00/home/H41070/EspaceTravail/TestCOCAGNE/DKLibs_Gauss_TH/DKLIB_EXTEND_MOX/Rho_5_5_mox"    
    dklibName3 = "/local00/home/H41070/EspaceTravail/TestCOCAGNE/DKLibs_Gauss_TH/DKLIB_EXTEND_MOX/B10_7_mox"  
    dklibName4 = "/local00/home/H41070/EspaceTravail/TestCOCAGNE/DKLibs_Gauss_TH/DKLIB_EXTEND_MOX/Xe_5_mox"


    dklibReference = "/local00/home/H41070/EspaceTravail/TestCOCAGNE/DKLibs_Gauss_TH/DKLIB_EXTEND_MOX/reference_CCpoints_mox" 

    listOfDKlibsForKLDecomposition = [dklibName0, dklibName1, dklibName2, dklibName3, dklibName4]
    
    return  dklibName_fullGrid, dklibName_fullGrid_reference, dklibReference, listOfDKlibsForKLDecomposition
##=================================================== 



   
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
    """
        For Paris 6
    """
    testCase = "MOX"
    nameFile_infor = "GeneralInfor_" +   testCase + ".txt"  
    File_infor = open(nameFile_infor,'w')
    ##===================================================
    dklibName_fullGrid, dklibName_fullGrid_reference, dklibReference, listOfDKlibsForKLDecomposition\
    = getListOfDklibs_MOX() ### other choices: getListOfDklibs_UOX(); getListOfDklibs_MOX()
    ##===================================================
    dimension = len(listOfDKlibsForKLDecomposition)
    ##===================================================
    ##   on recupere les valeurs des points des axes (axe par axe)  grilles multilinéaire et grille de réference
    ##===================================================
    if len(listOfDKlibsForKLDecomposition) == 5 :
        bu_values,  tc_values, density_values, cb_values, xe_values, ct_values  = crossSections.getFeedbackMesh(dklibName_fullGrid_reference)
        ref_values = crossSections.getFeedbackMesh(dklibReference)
        
    elif len(listOfDKlibsForKLDecomposition) == 6 :
        bu_values, tc_values, ct_values, density_values,  cb_values,  xe_values = crossSections.getFeedbackMesh(dklibName_fullGrid_reference)
    ##===================================================
    ## On détermine le nombre de paramètre de CRN
    ## On met les valeurs dans la grille de référence dans une liste avec la clé est l'indice de l'axe en chaine de caractère : '0', '1',...
    ##===================================================
    dimension = len(listOfDKlibsForKLDecomposition)   
    referencePoints = {}
    
    for Axis_k in range(dimension):
        referencePoints[str(Axis_k)] = ref_values[Axis_k]

    ##===================================================
    ## Les types de sections dans le test
    ##===================================================
   
    listOfCrossSectionNames = [ "macro_totale", "macro_absorption",  "macro_scattering", "macro_nu*fission",  "macro_fission"]  

    ##===================================================
    ## On récupère les sections efficaces avec ses arguments (goupe d'énergy, anisotropy, etc) pour
    ## + la liste des noms précis de ces sections (avec le numero du groupe et l'anisotropie): listOfCrossSectionNamesWithArgs
    ## + l'interpolation multilinéaire : listOfFinalCrossSectionOnFullGrid
    ## + la grille de référence pour mesurer erreur d'évaluation : listOfFinalCrossSectionOnFullGrid_reference
    ## + la grille pour les coefficients de Tucker : listOfRefCrossSectionsToDetermineCoeffs
    ## + la décomposition de Tucker : listOfCrossSectionsForTuckerGridsFromDklibs
    ##===================================================
    listOfCrossSectionNamesWithArgs, listOfFinalCrossSectionOnFullGrid,  listOfFinalCrossSectionOnFullGrid_reference,\
    listOfRefCrossSectionsToDetermineCoeffs,listOfCrossSectionsForTuckerGridsFromDklibs\
    = getAllListOfFinalCrossSection(listOfCrossSectionNames, dklibName_fullGrid, dklibName_fullGrid_reference, dklibReference, listOfDKlibsForKLDecomposition)

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
    File_infor.write("===== GENERAL INFORMATION =====\n")
    File_infor.write("\n")
    File_infor.write("Test case: %s\n" %(testCase))
    File_infor.write("\n")
    File_infor.write("listOfCrossSectionNamesWithArgs %s\n" %(listOfCrossSectionNamesWithArgs))
    File_infor.write("\n")
    File_infor.write("1. On the discretization of each grid used for the KL decomposition:\n")
    File_infor.write("\n")
    File_infor.write("listOfSubdivisionDirection = %s\n" %(listOfSubdivisionDirection)) 
    File_infor.write("\n")
    File_infor.write("listOfValuesForSubdivision = %s\n" %(listOfValuesForSubdivision))
    File_infor.write("\n")
    File_infor.write("listOfNumberOfPointsForSubdivision = %s\n" %(listOfNumberOfPointsForSubdivision))
    File_infor.write("\n")
    File_infor.write("listOfDomainBorders = %s\n" %(listOfDomainBorders))  
    File_infor.write("\n")
    File_infor.write("listOfNbPointsOnEachIntervalForFinalDiscretization = %s\n" %(listOfNbPointsOnEachIntervalForFinalDiscretization))
    File_infor.write("\n")
    File_infor.write("listOfTuckerGridNodesFromDKlibs = %s" %(listOfTuckerGridNodesFromDKlibs))
    File_infor.write("\n")
    File_infor.write("listOfIntegralMethods = %s\n" %(listOfIntegralMethods)) 
    File_infor.write("\n")
    #File_infor.write("listOfCrossSectionsForTuckerGridsFromDklibs = %s" %(listOfCrossSectionsForTuckerGridsFromDklibs))
    File_infor.write("\n")
   
    
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
   
    
    #File_infor.close()
    
    ### Recall these files in KLdecomposition.py and Tuckerdecomposition.py for each csName 
    File_infor_cs = {} 
    for csName in listOfCrossSectionNamesWithArgs:
        nameFile_infor_cs = "TuckerInfor_" +   csName + ".txt"  
        File_infor_cs[csName] = open(nameFile_infor_cs,'w')
    #print File_infor_cs
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
    ### C'est facultatif pour verifier qu'on ne s'est pas trompé en rentrant les points dans GAB
    ### si jamais on s'est troompé, le code s'arrête avec un message d'erreur
    allCrossSections.verifyInterpolPointsInTuckerNodesAndReferenceGrid(dimension, listOfFinalEmpiricalInterpolationPoints, listOfTuckerGridNodes, referencePoints)   
  
    
    ###===========Extract fromm FinalInterpolationPoints the interpolation point for each section========================================
    listOfInterpolationPoints = allCrossSections.getListOfInterpolationPointsForEachSection(len(listOfDKlibsForKLDecomposition),\
        listOfCrossSectionNamesWithArgs, listOfBasisFcts, listOfNbOfBasisFcts, listOfFinalEmpiricalInterpolationPoints)
        
    #####============END WITH GREEDY============================
    
    #####===================================================
    for csName in listOfCrossSectionNamesWithArgs:
        print "====================" 
        print "csName = ", csName
        print "listOfInterpolationPoints[csName] ", listOfInterpolationPoints[csName]         
        print "listOfTuckerGridNodes", listOfTuckerGridNodes   
        print "===================="
        #raw_input()
        
    """
        For Paris 6
    """ 
    File_infor.write( "====================\n" )
    File_infor.write("2. On the EIM:\n")
    File_infor.write("\n")
    File_infor.write("listOfTuckerGridNodes = %s" %(listOfTuckerGridNodes))
    File_infor.write("\n")
    File_infor.write("listOfPointsForGreedyAlgo = %s\n" %(listOfPointsForGreedyAlgo))
    File_infor.write("\n")
    File_infor.write("listOfFinalEmpiricalInterpolationPoints = %s\n" %(listOfFinalEmpiricalInterpolationPoints))
    File_infor.write("\n")
    File_infor.write("listOfInterpolationPoints = %s\n" %(listOfInterpolationPoints))
    File_infor.write("\n")
    File_infor.write("====================\n" )
    
    for csName in listOfCrossSectionNamesWithArgs:
        File_infor.write( "====================\n" )
        File_infor.write("csName = %s\n" %(csName))
        File_infor.write("\n")
        File_infor.write("listOfInterpolationPoints[csName] = %s\n"%(listOfInterpolationPoints[csName]))
        File_infor.write("\n")
        File_infor.write( "====================\n" )
    File_infor.close()
    #####===================================================

    ####============Cese for only one cross-section in the "listOfCrossSectionNamesWithArgs" (already determined)=======================================
    #argument1 = sys.argv[1]    
    #iSection = int(argument1)
    #csName = listOfCrossSectionNamesWithArgs[iSection]
    #listOfCrossSectionNamesWithArgs = [csName]
    
    ### ATTENTION!!! ICI On ne calcule les coeefs de la décopostion que poiur les sections définies dans la liste listOfCrossSectionNamesWithArgs
    ### ici on a écrasé cette liste, pour la remplacer par unbe plus courte qui ne comprend qu'une seule section (c'est pour un test
    ### si on ne fait pas de test (mais vrai calcul en production) alors la ligne suivante est commentée
    #listOfCrossSectionNamesWithArgs = ["macro_scattering001"]#["macro_nu*fission0"] #["macro_scattering010"] #["macro_absorption1"] #["macro_scattering001"]
    ##===================================================
    print "==================="
    print "listOfCrossSectionNamesWithArgs = ", listOfCrossSectionNamesWithArgs
    print "==================="
    ####===================================================
    criterion_post = "posteriori"
    eps = 1E-6
    ###===================================================
    ### Déterminer toutes les décomposition de Tucker pour toutes les sections dans "listOfCrossSectionNamesWithArgs"
    ###===================================================
    listOfTuckerDecomp = {}
    listOfRefParameterCoordinates = {}
    for csName in listOfCrossSectionNamesWithArgs:    
        # les 2 ligne suivantes sont indispensables (calcul et stockage des coefficients)
        listOfTuckerDecomp[csName] = Tucker.Tucker(listOfBasisFcts[csName], listOfNbOfBasisFcts[csName], listOfRefCrossSectionsToDetermineCoeffs[csName],listOfInterpolationPoints[csName], csName)
        listOfTuckerDecomp[csName].getCoefficients()
        # les lignes suivantes c'est pour le creusage 
        #listOfTuckerDecomp[csName].getSparseCoefficients_posteriori (criterion_post, eps)
        #listOfTuckerDecomp[csName].getSparseCoefficients_priori_index()
        #listOfTuckerDecomp[csName].getIntersection_DirectEliminatedPoints_TensorStructurePoints()
    ###===================================================
    ###=====================================================================
    ###nameOfTest = "_sparse_eps1E_10_14700p_UOX_extend_critere4.4_60_30Percent_prioriDirect"
    #nameOfTest = "_sparse_eps1E_10_14700p_UOX_Gd_extend_critere4.4_20_40Percent_prioriDirect"
    
    ####=====================================================================
    #### Calculer et mesurer erreur de l'évaluation
    ####=====================================================================
    #listOfValues_Tucker = {}
    #listOfValues_Cocagne = {}
    #listOfValues_AP2= {}
    
    #listOfValues_Tucker_priori = {}
    #listOfValues_Tucker_posteriori = {}
    
    ####===================================================

    #listOfMaxRelativeError_Tucker_kinf = {}
    #listOfMSE_Tucker_kinf = {}
    
    #listOfMaxRelativeError_Cocagne_kinf = {}
    #listOfMSE_Cocagne_kinf = {}
    
    #listOfMaxRelativeError_Tucker_kinf_priori = {}
    #listOfMSE_Tucker_kinf_priori = {}
    
    #listOfMaxRelativeError_Tucker_kinf_posteriori = {}
    #listOfMSE_Tucker_kinf_posteriori = {}

    #print " listOfCrossSectionNamesWithArgs =",  listOfCrossSectionNamesWithArgs
    
    #for csName in listOfCrossSectionNamesWithArgs:        
  
        #count = 1
        #values_Tucker = []
        #values_Cocagne = []
        #values_AP2 = []
        
        #values_Tucker_priori = []
        #values_Tucker_posteriori = []
        ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        #for bu in bu_values:
            #for cb in cb_values:
                #for xe in xe_values:
                    #for density in density_values:
                        #for tc in tc_values:
                            #if ct_values:
                                #for ct in ct_values:
                                    #point = (bu, tc, ct, density ,cb, xe)
                            #else:
                                #point = (bu, tc, density ,cb, xe)
                           
                           
                            #value_T = listOfTuckerDecomp[csName].evaluate(point)
                            #value_C = listOfFinalCrossSectionOnFullGrid[csName].evaluate(point)
                            #value_A = listOfFinalCrossSectionOnFullGrid_reference[csName].evaluate(point)
                            ####=================================
                            #### A changer le type de creusage a priori ici pour tester la précision de creusage a priori
                            ####=================================
                            #value_T_priori_direct, value_T_priori_ls, value_T_priori_ls_tensorStructure = listOfTuckerDecomp[csName].evaluate_sparse_priori(point)  
                            
                            ####=================================
                            #### For a case where we are interessed in ls_tensorStructure
                            ####=================================
                            
                            #### Choix 1
                            #value_T_priori = value_T_priori_direct 
                            
                            #### Choix 2
                            ##value_T_priori = value_T_priori_ls  
                            
                            #### Choix 3 (résultat semble mauvais)
                            ##value_T_priori = value_T_priori_ls_tensorStructure                           
                            
                            #### Choix 4 (résultat semble mauvais)
                            ##value_T_priori = listOfTuckerDecomp[csName].evaluate_sparse_priori_tensor_direct(point)
                            ####=================================
                            
                            #value_T_posteriori = listOfTuckerDecomp[csName].evaluate_sparse_posteriori(point) 

                            #values_Tucker.append(value_T)
                            #values_Cocagne.append(value_C)
                            #values_AP2.append(value_A)
                            
                            #values_Tucker_priori.append(value_T_priori)
                            #values_Tucker_posteriori.append(value_T_posteriori)
                            
                            #####%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            #print "count =", count
                            #print "Tucker = %s, AP2 = %s, Cocagne = %s, Tucker_priori=%s, Tucker_posteriori = %s" %(values_Tucker[count-1].real,values_AP2[count-1], values_Cocagne[count-1],\
                            #values_Tucker_priori[count-1], values_Tucker_posteriori[count-1])
                            #print "Tucker - AP2 = ", values_Tucker[count-1] - values_AP2[count-1]
                            #print "Cocagne - AP2 =",  values_Cocagne[count-1] - values_AP2[count-1]
                            #print "=================="
                            #print "abs((Tucker - AP2)/AP2) = ", abs((values_Tucker[count-1] - values_AP2[count-1])/values_AP2[count-1])
                            #print "abs((Tucker_prio - AP2)/AP2) = ", abs((values_Tucker_priori[count-1] - values_AP2[count-1])/values_AP2[count-1])
                            #print "abs((Tucker_post - AP2)/AP2) = ", abs((values_Tucker_posteriori[count-1] - values_AP2[count-1])/values_AP2[count-1]) 
                            #print"==================="
                            #print "abs((Cocagne - AP2)/AP2) =",  abs((values_Cocagne[count-1] - values_AP2[count-1])/values_AP2[count-1])
                            #raw_input()
                            #####%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            #count = count + 1
                            
        #listOfValues_Tucker[csName] = np.asarray(values_Tucker)
        #listOfValues_Cocagne[csName] = np.asarray(values_Cocagne)
        #listOfValues_AP2[csName]= np.asarray(values_AP2) 
        #listOfValues_Tucker_priori[csName] = np.asarray(values_Tucker_priori) 
        #listOfValues_Tucker_posteriori[csName] = np.asarray(values_Tucker_posteriori)
                      
      
        
    #COUNT = count -1    
    #print "COUNT =", COUNT    
    ####======================Kinf=============================     
    #listOfRelativeError_Tucker = {}
    
    #listOfRelativeError_Tucker_priori = {}
    #listOfRelativeError_Tucker_posteriori = {}
    #listOfRelativeError_Cocagne = {}
    
    #listOfMaxRelativeError_Tucker = {}
    #listOfMaxRelativeError_Tucker_priori = {}
    #listOfMaxRelativeError_Tucker_posteriori = {}
    
    #listOfMaxRelativeError_Cocagne = {}
     
    #listOfMSE_Tucker = {}  
    #listOfMSE_Tucker_priori = {} 
    #listOfMSE_Tucker_posteriori = {} 
    #listOfMSE_Cocagne = {}
    ####======================Kinf============================= 
    #for csName in listOfCrossSectionNamesWithArgs:
        #nameFile_Err = "RelativeError_AP2_Cocagne_Tucker_TuckerPrio_TuckerPosterio" + nameOfTest +"_" + csName        
        #File_Err = open(nameFile_Err,'w')
         #### 4/3/2016: normaliser e_relative = |f-f~|/maxf
        #maxAP2 = max(abs(listOfValues_AP2[csName]))
        #listOfRelativeError_Tucker[csName] = 1E5*(listOfValues_Tucker[csName] - listOfValues_AP2[csName])/maxAP2
        
        #listOfRelativeError_Tucker_priori[csName] = 1E5*(listOfValues_Tucker_priori[csName] - listOfValues_AP2[csName])/maxAP2
        #listOfRelativeError_Tucker_posteriori[csName] = 1E5*(listOfValues_Tucker_posteriori[csName] - listOfValues_AP2[csName])/maxAP2
        
        #listOfRelativeError_Cocagne[csName] = 1E5*(listOfValues_Cocagne[csName] - listOfValues_AP2[csName])/maxAP2
        ####==========Write the results of relative error of two methods in files===========
        #count = 1
        #for bu in bu_values:
            #for cb in cb_values:
                #for xe in xe_values:
                    #for density in density_values:
                        #for tc in tc_values:
                            #value_A = listOfValues_AP2[csName][count-1]
                            #value_C = listOfValues_Cocagne[csName][count-1]
                            #value_T = listOfValues_Tucker[csName][count-1]
                            
                    
                            
                            #value_T_priori = listOfValues_Tucker_priori[csName][count-1]
                            #value_T_posteriori = listOfValues_Tucker_posteriori[csName][count-1]
                            
                            #err_C = 1E5*(value_C - value_A)/maxAP2
                            #err_T = 1E5*(value_T - value_A)/maxAP2
                            #err_T_priori = 1E5*(value_T_priori - value_A)/maxAP2
                            #err_T_posteriori = 1E5*(value_T_posteriori - value_A)/maxAP2
                            
                            #File_Err.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"\
                            #%(count,bu,tc,density,cb,xe,value_A,value_C, value_T.real,value_T_priori.real,value_T_posteriori.real, err_C, err_T.real,err_T_priori.real, err_T_posteriori.real))
                            #count = count +1
        #File_Err.close()
        ####==========Write the results of relative error of two methods in files===========
        
        #####=====================
        
        #listOfMaxRelativeError_Tucker[csName] = max(abs(listOfRelativeError_Tucker[csName])) #max(abs((listOfValues_Tucker[csName]/listOfValues_AP2[csName] - 1)*1E5))
        #listOfMaxRelativeError_Tucker_priori[csName] = max(abs(listOfRelativeError_Tucker_priori[csName])) #max(abs((listOfValues_Tucker[csName]/listOfValues_AP2[csName] - 1)*1E5))
        #listOfMaxRelativeError_Tucker_posteriori[csName] = max(abs(listOfRelativeError_Tucker_posteriori[csName])) #max(abs((listOfValues_Tucker[csName]/listOfValues_AP2[csName] - 1)*1E5))
        
        #listOfMaxRelativeError_Cocagne[csName] = max(abs(listOfRelativeError_Cocagne[csName]))#max(abs((listOfValues_Cocagne[csName]/listOfValues_AP2[csName] - 1)*1E5))
        
        
        #AP2_RS = np.sqrt(sum(np.power(listOfValues_AP2[csName],2)))
        
        #listOfMSE_Tucker[csName]  = np.sqrt(sum(np.power(abs(listOfValues_Tucker[csName] - listOfValues_AP2[csName])*1E5, 2)))/AP2_RS
        #listOfMSE_Tucker_priori[csName]  = np.sqrt(sum(np.power(abs(listOfValues_Tucker_priori[csName] - listOfValues_AP2[csName])*1E5, 2)))/AP2_RS
        #listOfMSE_Tucker_posteriori[csName]  = np.sqrt(sum(np.power(abs(listOfValues_Tucker_posteriori[csName] - listOfValues_AP2[csName])*1E5, 2)))/AP2_RS

        
        #listOfMSE_Cocagne[csName]  = np.sqrt(sum(np.power(abs(listOfValues_Cocagne[csName] - listOfValues_AP2[csName])*1E5, 2)))/AP2_RS
        
        #print "====================="
        #print csName
        #print "relativeError_Tucker", listOfMaxRelativeError_Tucker[csName]
        #print "relativeError_Tucker_priori", listOfMaxRelativeError_Tucker_priori[csName]
        #print "relativeError_Tucker_posteriori", listOfMaxRelativeError_Tucker_posteriori[csName]        
        #print "relativeError_Cocagne", listOfMaxRelativeError_Cocagne[csName]
        
        #print "MSE_Tucker", listOfMSE_Tucker[csName]
        #print "MSE_Tucker_priori", listOfMSE_Tucker_priori[csName]
        #print "MSE_Tucker_posteriori", listOfMSE_Tucker_posteriori[csName]
        #print "MSE_Cocagne", listOfMSE_Cocagne[csName]
        #print "====================="

    ####====================Reactivity exact========================  
    ####============================================
     #### Calcul très pratique mais le résultat est en doute
    
    #kinf_Tucker =  computeKinf(listOfValues_Tucker) 
    #kinf_Tucker_priori =  computeKinf(listOfValues_Tucker_priori)
    #kinf_Tucker_posteriori =  computeKinf(listOfValues_Tucker_posteriori)
                          
    #kinf_AP2 =  computeKinf(listOfValues_AP2) 
                              
    #kinf_Cocagne =  computeKinf(listOfValues_Cocagne) 
    ####============================================     

    
    #reactivityError_Tucker = (1.0/kinf_Tucker - 1.0/kinf_AP2)*1E5
    #reactivityError_Tucker_priori = (1.0/kinf_Tucker_priori - 1.0/kinf_AP2)*1E5
    #reactivityError_Tucker_posteriori = (1.0/kinf_Tucker_posteriori - 1.0/kinf_AP2)*1E5
    
    #reactivityError_Cocagne =(1.0/kinf_Cocagne - 1.0/kinf_AP2)*1E5
    
    ####==========Write the results of the difference for reactivity between two methods into files===========
    
    #nameFile_reactivity = "ReactivityError_AP2_Cocagne_Tucker_TuckerPrio_TuckerPosterio" + nameOfTest
    
    #File_ReactivityErr = open(nameFile_reactivity,'w')
    #count = 1
    
    #sumReactivity_A = 0
    
    #for bu in bu_values:
        #for cb in cb_values:
            #for xe in xe_values:
                #for density in density_values:
                    #for tc in tc_values:
                        #reactivity_A = kinf_AP2[count-1]
                        #reactivity_C = kinf_Cocagne[count-1]
                        #reactivity_T = kinf_Tucker[count-1]
                        #reactivity_T_priori = kinf_Tucker_priori[count-1]
                        #reactivity_T_posteriori = kinf_Tucker_posteriori[count-1]
                        
                        #sumReactivity_A = sumReactivity_A  + reactivity_A*reactivity_A 

                        #err_C = (1.0/reactivity_C - 1.0/reactivity_A)*1E5
                        #err_T = (1.0/reactivity_T - 1.0/reactivity_A)*1E5
                        #err_T_priori = (1.0/reactivity_T_priori - 1.0/reactivity_A)*1E5
                        #err_T_posteriori = (1.0/reactivity_T_posteriori - 1.0/reactivity_A)*1E5
                        
                        #File_ReactivityErr.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
                        #%(count,bu,tc,density,cb,xe,reactivity_A,reactivity_C, reactivity_T.real, reactivity_T_priori.real,reactivity_T_posteriori.real,err_C, err_T.real, err_T_priori.real, err_T_posteriori.real))
                        #count = count +1
    
    #File_ReactivityErr.close()
    ####=====================
    
    
    #maxReactivityAbsError_Tucker = max(abs(reactivityError_Tucker))
    #maxReactivityAbsError_Tucker_priori = max(abs(reactivityError_Tucker_priori))
    #maxReactivityAbsError_Tucker_posteriori = max(abs(reactivityError_Tucker_posteriori))
    
    #maxReactivityAbsError_Cocagne = max(abs(reactivityError_Cocagne))
    
    
    #reactivity_L2 = np.sqrt(sumReactivity_A)
    
    #reactivityMSE_Tucker = np.sqrt((sum(reactivityError_Tucker*reactivityError_Tucker)))/reactivity_L2
    #reactivityMSE_Tucker_priori = np.sqrt((sum(reactivityError_Tucker_priori*reactivityError_Tucker_priori)))/reactivity_L2
    #reactivityMSE_Tucker_posteriori = np.sqrt((sum(reactivityError_Tucker_posteriori*reactivityError_Tucker_posteriori)))/reactivity_L2
    #reactivityMSE_Cocagne = np.sqrt((sum(reactivityError_Cocagne*reactivityError_Cocagne)))/reactivity_L2
    
    #print "====================="
    #print "len(reactivityError_Tucker)", len(reactivityError_Tucker)
 
    #print "maxReactivityAbsError_Tucker = ", maxReactivityAbsError_Tucker
    #print "maxReactivityAbsError_Tucker_priori = ", maxReactivityAbsError_Tucker_priori
    #print "maxReactivityAbsError_Tucker_posteriori = ", maxReactivityAbsError_Tucker_posteriori
    #print "maxReactivityAbsError_Cocagne = ", maxReactivityAbsError_Cocagne
    
    #print "reactivityMSE_Tucker = ", reactivityMSE_Tucker
    #print "reactivityMSE_Tucker_priori = ", reactivityMSE_Tucker_priori
    #print "reactivityMSE_Tucker_posteriori = ", reactivityMSE_Tucker_posteriori
    
    #print "reactivityMSE_Cocagne = ", reactivityMSE_Cocagne    
    #print "====================="
    
    ####============================================