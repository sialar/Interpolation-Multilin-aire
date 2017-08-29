# -*- coding: utf-8 -*-
import os
#import dkmodel
#import dktools
#import dkbase
#importnp.array
import datetime
import numpy as np
import sys
from numpy import linalg as LA
import scipy
import random
#import crossSections
#import NeutronicLibraryManager
#import Symbols
#import geometricElements
#from   Geometry import createFlatRegionName
#import geometricElements.datarodbanks
#import Geometry
#loadFromDKLib = NeutronicLibraryManager.loadFromDKLib
###===================================================
### Nouveaux programmes:
import multiDimensIntegralMethods as M
import KLdecomposition as KL
import discretizationChoice as discret
import Tuckerdecomposition as Tucker
import LagrangeInterpolation as LI
#import allCrossSections
import LagrangePolynomial
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

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx]

"""
Des changements dans cette fct "getListOfInterpolationFcts"
"""
def getListOfInterpolationFcts(listOfBorders, listOfTuckerGridNodes_Axis_k, finalOrthNormalizedEigVects_Axis_k, Axis_k):   
      
        """
        Fonction qui reconstruit les vecteurs propres gardés dans l'axis "Axis_k" par l'interpolation de Lagrange :
        + Une l'interpolation de Lagrange si il n'y a pas de coupage
        +  Des interpolations de Lagrange par morceau si il y a de coupages
        Argument:
             rien
        retourne:
            rien
        """  
               
        listOfBasicFctsUsingLagrangeInterpolation_Axis_k = []
        nbExtremes =  len(listOfBorders[str(Axis_k)])
        nbOfFcts_Axis_k = len(finalOrthNormalizedEigVects_Axis_k[str(Axis_k)])
        #print nbOfFcts_Axis_k
        #raw_input()
        #print "Axis_k =", Axis_k
        for i in range(nbOfFcts_Axis_k):
            ### Define each sub-polynomial/sub-interval if the main interval of "Axis_k" is subdivided 
            #print "i =", i
            if nbExtremes > 2 :
                polLagrange_ki = []
                
                for j in range(nbExtremes - 1) :
                    #print listOfTuckerGridNodes_Axis_k[str(Axis_k)]
                    #print listOfBorders[str(Axis_k)]
                    #raw_input()
                    #j1Arr = listOfTuckerGridNodes_Axis_k[str(Axis_k)].index(listOfBorders[str(Axis_k)][j])
                    
                    #j2Arr = listOfTuckerGridNodes_Axis_k[str(Axis_k)].index(listOfBorders[str(Axis_k)][j+1])
                    ### Finding the nearest value ensures that j1Arr <> [] and j2Arr <> []
                    value1 = find_nearest(listOfTuckerGridNodes_Axis_k[str(Axis_k)], listOfBorders[str(Axis_k)][j])
                    value2 = find_nearest(listOfTuckerGridNodes_Axis_k[str(Axis_k)], listOfBorders[str(Axis_k)][j+1])
                    
                    j1Arr =np.where(listOfTuckerGridNodes_Axis_k[str(Axis_k)] == value1)[0].tolist()
                    
                    j2Arr = np.where(listOfTuckerGridNodes_Axis_k[str(Axis_k)] == value2)[0].tolist()
                    
                    #j1Arr =np.where(listOfTuckerGridNodes_Axis_k[str(Axis_k)] == listOfBorders[str(Axis_k)][j])[0].tolist()
                    #j2Arr = np.where(listOfTuckerGridNodes_Axis_k[str(Axis_k)] == listOfBorders[str(Axis_k)][j+1])[0].tolist()
                    
                    #print j
                    #print j1Arr
                    #print j2Arr
                   
                    
                    
                    j1 = j1Arr[-1]
                    j2 = j2Arr[0]
                    j2 = j2 + 1
                    
                    
                    #print j1
                    #print j2
                    
                    #raw_input()
                    ###because [j1:j2] take the values before j2
                    
                    interp_k_couplage =  LagrangePolynomial.LagrangePolynomial(listOfTuckerGridNodes_Axis_k[str(Axis_k)][j1:j2],\
                                              finalOrthNormalizedEigVects_Axis_k[str(Axis_k)][i][j1:j2].T)
                          
                    polLagrange_ki.append(interp_k_couplage)
                    
            ### Define the polynomial if the main interval of Axis_k is NOT subdivided         
            elif  nbExtremes == 2: 
                #print "i = ", i
                #print finalOrthNormalizedEigVects_Axis_k[str(Axis_k)][i]
                #print "Axis_k =", Axis_k  
                #raw_input()
                polLagrange_ki = [ LagrangePolynomial.LagrangePolynomial(listOfTuckerGridNodes_Axis_k[str(Axis_k)],\
                                        finalOrthNormalizedEigVects_Axis_k[str(Axis_k)][i])] 
                
            listOfBasicFctsUsingLagrangeInterpolation_Axis_k.append(polLagrange_ki)
    
        return listOfBasicFctsUsingLagrangeInterpolation_Axis_k
    
    
def evaluate(listOfBasisFcts, listOfFinalCoefIndexes_arr, FinalTuckerCoeffs, point) :
        """
      Fonction qui évalue la décompostion de Tucker pour la section considérée en un point donné
      Argument :
        "point" [un list de "d" (=dimension) flottants] : point sur le quel on évalue
      Retour :
        "approxValue" [un flottant]: la veleur d'évaluation de section sur le "point" via Tucker

        """
        approxValue = 0.0
   
        for I in range(len(listOfFinalCoefIndexes_arr)) :
            d = 1.0
           
            for Axis_k in range(dimension) :                
                fct = listOfBasisFcts[str(Axis_k)][listOfFinalCoefIndexes_arr[I][Axis_k]]
                d = d*LI.getInterpolation(fct, point[Axis_k])
               
            approxValue = approxValue +  FinalTuckerCoeffs[I]*d       
     
        return approxValue
###===================================================  
if __name__ == "__main__":
    
    dimension = 5
    
    listOfNominalValues = {'fuel_temperature' : 550.0, 'coolant_density': 0.71692097187, 'b10_concentration': 3.97352914661e-06, 'xenon_level': 1, 'coolant_temperature': 304.600006104}
    listOfAxis_Parameter_5dimension = {'0': 'BURNUP', '1': 'fuel_temperature', '2': 'coolant_density', '3': 'b10_concentration', '4': 'xenon_level'}


   
    listOfCrossSectionNames = [ "macro_totale", "macro_absorption",  "macro_scattering", "macro_nu*fission",  "macro_fission"]  

    ##===================================================
    listOfCrossSectionNamesWithArgs = ['macro_totale0', 'macro_totale1', 'macro_absorption0', 'macro_absorption1', 'macro_scattering000', 'macro_scattering001', 'macro_scattering010', 'macro_scattering011', 'macro_nu*fission0', 'macro_nu*fission1', 'macro_fission0', 'macro_fission1']

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
    #listOfIntegralMethods = {'1': ['T', 'CC', 'T', 'T', 'T'], '0': ['CCbis', 'T', 'T', 'T', 'T'], '3': ['T', 'T', 'T', 'CC', 'T'], '2': ['T', 'T', 'CCbis', 'T', 'T'], '4': ['T', 'T', 'T', 'T', 'CC']}
 ###=========================================  
    #listOfTuckerGridNodes = {'1': np.array([20.        ,   280.67496474,   910.        ,  1539.32503526,
        #1800.]), '0': np.array([  0.00000000e+00,   7.50000000e+01,   1.50000000e+02,
         #1.50000000e+02,   5.24893302e+02,   1.59249910e+03,
         #3.19028410e+03,   5.07500000e+03,   6.95971590e+03,
         #8.55750090e+03,   9.62510670e+03,   1.00000000e+04,
         #1.00000000e+04,   1.26642164e+04,   2.02512627e+04,
         #3.16060799e+04,   4.50000000e+04,   5.83939201e+04,
         #6.97487373e+04,   7.73357836e+04,   8.00000000e+04]), '3': np.array([  0.00000000e+00,   1.17117440e-06,   4.37088238e-06,
         #8.74176476e-06,   1.31126471e-05,   1.63123551e-05,
         #1.74835295e-05]), '2': np.array([ 0.40000001,  0.44641201,  0.55846049,  0.67050897,  0.71692097,
        #0.71692097,  0.75837694,  0.85846049,  0.95854404,  1.        ]), '4': np.array([  9.99999997e-07,   2.92894072e-01,   1.00000050e+00,
         #1.70710693e+00,   2.00000000e+00])}
        
    
    listOfBorders = {}
    listOfBorders['0'] = np.array([0.0, 150.0, 10000.0, 80000.0])
    listOfBorders['1'] = np.array([20.0, 1800.0])
    listOfBorders['2'] = np.array([0.4000000059604645, 0.71692097187, 1.0])
    listOfBorders['3'] = np.array([0.0, 1.748352951835841e-05])
    listOfBorders['4'] = np.array([9.999999974752427e-07, 2.0])
    
 ###=========================================      
    listOfTuckerGridNodes_Axis_k = {}
    listOfTuckerGridNodes_Axis_k['0'] = np.array([0.00000000e+00, 7.50000000e+01, 1.50000000e+02, 1.50000000e+02
, 5.24893302e+02, 1.59249910e+03, 3.19028410e+03, 5.07500000e+03
, 6.95971590e+03, 8.55750090e+03, 9.62510670e+03, 1.00000000e+04
, 1.00000000e+04, 1.26642164e+04, 2.02512627e+04, 3.16060799e+04
, 4.50000000e+04, 5.83939201e+04, 6.97487373e+04, 7.73357836e+04
, 8.00000000e+04])
    
    listOfTuckerGridNodes_Axis_k['1'] = np.array([20.,  280.67496474, 910., 1539.32503526,  1800.])
    
    listOfTuckerGridNodes_Axis_k['2'] = np.array([0.40000001,  0.44641201,  0.55846049,  0.67050897,  0.71692097,  0.71692097,
  0.75837694,  0.85846049,  0.95854404,  1.        ])
    
    listOfTuckerGridNodes_Axis_k['3'] = np.array([  0.00000000e+00,   1.17117440e-06,   4.37088238e-06,   8.74176476e-06,
   1.31126471e-05,   1.63123551e-05,   1.74835295e-05])
    
    listOfTuckerGridNodes_Axis_k['4'] =  np.array([  9.99999997e-07,   2.92894072e-01,   1.00000050e+00,   1.70710693e+00,
   2.00000000e+00])
    
    
    #idx = find_nearest(listOfTuckerGridNodes_Axis_k['2'],listOfBorders['2'][2])
    #print "idx =", idx
    
 ###=========================================   
    
    finalOrthNormalizedEigVects_Axis_k = {}
    finalOrthNormalizedEigVects_Axis_k['0'] = np.array([[3.46498263e-03,  3.46532333e-03,  3.46565629e-03,  3.46565629e-03
, 3.46679117e-03,  3.46844034e-03,  3.47128424e-03,  3.47518875e-03
, 3.47929429e-03,  3.48280065e-03,  3.48512776e-03,  3.48594255e-03
, 3.48594255e-03,  3.49164866e-03,  3.50717035e-03,  3.52842780e-03
, 3.54969334e-03,  3.56553898e-03,  3.57477798e-03,  3.57933995e-03
, 3.58073136e-03],
 [6.87301991e-03,  6.86352325e-03,  6.85062966e-03,  6.85062966e-03
, 6.77737960e-03,  6.54526922e-03,  6.15155310e-03,  5.68286523e-03
, 5.24438151e-03,  4.90038391e-03,  4.68262138e-03,  4.60921043e-03
, 4.60921043e-03,  4.11110644e-03,  2.86225569e-03,  1.18027213e-03
,  -7.90945351e-04, -2.81591381e-03, -4.44876909e-03, -5.42834473e-03
,  -5.75044808e-03],
 [6.67938852e-03,  6.74154999e-03,  6.91468772e-03,  6.91468772e-03
, 7.66802003e-03,  7.92296198e-03,  7.32414763e-03,  6.16020942e-03
, 4.81019634e-03,  3.66303630e-03,  2.91649586e-03,  2.67095709e-03
, 2.67095709e-03,  9.72156049e-04, -2.57330377e-03, -4.56124241e-03
,  -2.93180230e-03,  6.35917158e-04,  3.00761086e-03,  3.66274792e-03
, 3.70048465e-03],
 [1.29890730e-02,  1.28479577e-02,  1.25672386e-02,  1.25672386e-02
, 1.06809626e-02,  7.26143746e-03,  4.02536451e-03,  1.26194015e-03
,  -9.08636853e-04, -2.33245963e-03, -3.12724182e-03, -3.37500720e-03
,  -3.37500720e-03, -4.64696579e-03, -4.57082127e-03,  7.47181492e-05
, 4.01798231e-03,  2.50613908e-03, -1.60566753e-03, -4.34284444e-03
,  -5.12086915e-03],
 [ -2.36750282e-02, -2.33397535e-02, -2.20712166e-02, -2.20712166e-02
,  -1.53167489e-02, -7.01586245e-03, -3.49631537e-04,  3.65310984e-03
, 5.06494376e-03,  5.11241858e-03,  4.78893567e-03,  4.66977111e-03
, 4.66977111e-03,  3.10635268e-03, -1.32077520e-03, -2.66618004e-03
, 1.15254402e-03,  3.32437441e-03,  2.17285588e-04, -5.37503284e-03
,  -8.85212993e-03]])
 
    finalOrthNormalizedEigVects_Axis_k['1'] = np.array([[-0.02366172, -0.02367751, -0.02370469, -0.02372473, -0.02373202], [0.04889343, 0.02981867, -0.00270514, -0.02726803, -0.03624861]])
    
    finalOrthNormalizedEigVects_Axis_k['2'] = np.array([[0.9124528, 0.96873276, 1.10422638, 1.23920536, 1.29497711, 1.29497711
, 1.34472976, 1.46460765, 1.58412178, 1.63350581]
, [-2.65011783, -2.21655361, -1.21376975, -0.28704477, 0.06436376, 0.06436376
, 0.35742806, 0.98475257, 1.50382454, 1.67961492]
, [1.61809818, 1.1795946, -0.21132772, -1.2204461, -1.41361137, -1.41361137
, -1.41304145, -0.54829292, 2.10902196, 4.15249373]])
 
    finalOrthNormalizedEigVects_Axis_k['3'] = np.array([[238.12409505, 238.2882373,  238.70480074, 239.20487502, 239.63436791,  239.90883497, 240.00186117],
 [394.56259325, 347.15677329, 210.54671936, 9.61004127, -204.86511925
, -368.86205624, -430.03103087],
 [567.05243791, 350.43992792, -67.21762983, -264.02969217, -70.63019723
,  316.10251273, 529.14638765]])
 
    finalOrthNormalizedEigVects_Axis_k['4'] = np.array([[-0.70627887, -0.70653249, -0.70711869, -0.70766922, -0.70788737]
, [-1.24608698, -0.87131136, 0.00957682, 0.8604891,  1.20381659]])
    
    
    listOfPointsForGreedyAlgo = {'1':np.array([   20.        ,   280.67496474,   910.        ,  1539.32503526,
        1800.        ]), '0':np.array([  0.00000000e+00,   7.50000000e+01,   1.50000000e+02,
         1.50000000e+02,   5.24893302e+02,   1.59249910e+03,
         3.19028410e+03,   5.07500000e+03,   6.95971590e+03,
         8.55750090e+03,   9.62510670e+03,   1.00000000e+04,
         1.00000000e+04,   1.26642164e+04,   2.02512627e+04,
         3.16060799e+04,   4.50000000e+04,   5.83939201e+04,
         6.97487373e+04,   7.73357836e+04,   8.00000000e+04]), '3':np.array([  0.00000000e+00,   1.17117440e-06,   4.37088238e-06,
         8.74176476e-06,   1.31126471e-05,   1.63123551e-05,
         1.74835295e-05]), '2':np.array([ 0.40000001,  0.44641201,  0.55846049,  0.67050897,  0.71692097,
        0.71692097,  0.75837694,  0.85846049,  0.95854404,  1.        ]), '4':np.array([  9.99999997e-07,   2.92894072e-01,   1.00000050e+00,
         1.70710693e+00,   2.00000000e+00])}
###========================================= 
    listOfBasicFctsUsingLagrangeInterpolation = {}
    for Axis_k in range(dimension):
        basisFcts_Axis_k = getListOfInterpolationFcts(listOfBorders, listOfTuckerGridNodes_Axis_k, finalOrthNormalizedEigVects_Axis_k, Axis_k)
        
        listOfBasicFctsUsingLagrangeInterpolation[str(Axis_k)] = basisFcts_Axis_k
    
    #print "listOfBasicFctsUsingLagrangeInterpolation", listOfBasicFctsUsingLagrangeInterpolation
    #raw_input()
    
    
    
  ###=========================================    
    #listOfFinalEmpiricalInterpolationPoints = {'1': [20.0, 1800.0, 280.67496474397274, 910.0, 1539.3250352560271], '0': [0.0, 1592.4991026562534, 80000.0, 45000.0, 20251.262658470834, 69748.737341529166, 8557.5008973437471], '3': [0.0, 1.7483529518358409e-05, 8.7417647591792047e-06], '2': [1.0, 0.40000000596046448, 0.9585440361362102, 0.71692097187000003, 0.55846048891523226, 0.85846048593500002], '4': [9.9999999747524271e-07, 2.0, 1.0000004999999987]}

    #listOfInterpolationPoints = {'macro_totale1': {'1': [20.0, 1800.0, 910.0], '0': [0.0, 8557.5008973437471, 80000.0, 45000.0], '3': [0.0, 1.7483529518358409e-05], '2': [1.0, 0.40000000596046448, 0.85846048593500002, 0.71692097187000003], '4': [9.9999999747524271e-07, 2.0]}, 'macro_totale0': {'1': [20.0, 1800.0], '0': [0.0, 80000.0, 8557.5008973437471, 45000.0, 20251.262658470834], '3': [0.0, 1.7483529518358409e-05, 8.7417647591792047e-06], '2': [1.0, 0.40000000596046448, 0.71692097187000003], '4': [9.9999999747524271e-07, 2.0]}, 'macro_scattering000': {'1': [20.0, 1800.0], '0': [0.0, 80000.0, 8557.5008973437471, 45000.0, 20251.262658470834], '3': [0.0, 1.7483529518358409e-05, 8.7417647591792047e-06], '2': [1.0, 0.40000000596046448, 0.71692097187000003], '4': [9.9999999747524271e-07, 2.0]}, 'macro_scattering001': {'1': [20.0, 1800.0, 280.67496474397274, 910.0, 1539.3250352560271], '0': [0.0, 1592.4991026562534, 80000.0, 45000.0, 20251.262658470834, 69748.737341529166], '3': [0.0, 1.7483529518358409e-05, 8.7417647591792047e-06], '2': [0.40000000596046448, 1.0, 0.55846048891523226, 0.9585440361362102, 0.71692097187000003, 0.85846048593500002], '4': [9.9999999747524271e-07, 2.0, 1.0000004999999987]}, 'macro_scattering011': {'1': [20.0, 1800.0], '0': [0.0, 80000.0, 8557.5008973437471, 45000.0], '3': [0.0, 1.7483529518358409e-05], '2': [1.0, 0.40000000596046448, 0.85846048593500002, 0.55846048891523226], '4': [9.9999999747524271e-07, 2.0]}, 'macro_scattering010': {'1': [20.0, 1800.0, 910.0], '0': [0.0, 1592.4991026562534, 80000.0, 45000.0, 20251.262658470834, 69748.737341529166, 8557.5008973437471], '3': [0.0, 1.7483529518358409e-05, 8.7417647591792047e-06], '2': [1.0, 0.40000000596046448, 0.55846048891523226, 0.85846048593500002, 0.71692097187000003], '4': [9.9999999747524271e-07, 2.0]}, 'macro_nu*fission0': {'1': [20.0, 1800.0, 910.0], '0': [0.0, 80000.0, 8557.5008973437471, 45000.0, 20251.262658470834], '3': [1.7483529518358409e-05, 0.0, 8.7417647591792047e-06], '2': [1.0, 0.40000000596046448, 0.71692097187000003, 0.85846048593500002], '4': [9.9999999747524271e-07, 2.0]}, 'macro_nu*fission1': {'1': [20.0, 1800.0, 910.0], '0': [0.0, 20251.262658470834, 80000.0, 45000.0, 8557.5008973437471], '3': [0.0, 1.7483529518358409e-05, 8.7417647591792047e-06], '2': [1.0, 0.40000000596046448, 0.9585440361362102, 0.71692097187000003, 0.55846048891523226], '4': [9.9999999747524271e-07, 2.0]}, 'macro_fission0': {'1': [20.0, 1800.0, 910.0], '0': [0.0, 80000.0, 8557.5008973437471, 45000.0, 20251.262658470834], '3': [1.7483529518358409e-05, 0.0, 8.7417647591792047e-06], '2': [1.0, 0.40000000596046448, 0.85846048593500002, 0.55846048891523226], '4': [9.9999999747524271e-07, 2.0]}, 'macro_fission1': {'1': [20.0, 1800.0, 910.0], '0': [0.0, 20251.262658470834, 80000.0, 45000.0, 8557.5008973437471], '3': [0.0, 1.7483529518358409e-05, 8.7417647591792047e-06], '2': [1.0, 0.40000000596046448, 0.9585440361362102, 0.71692097187000003, 0.55846048891523226], '4': [9.9999999747524271e-07, 2.0]}, 'macro_absorption1': {'1': [20.0, 1800.0, 910.0], '0': [0.0, 80000.0, 8557.5008973437471, 45000.0, 20251.262658470834], '3': [0.0, 1.7483529518358409e-05, 8.7417647591792047e-06], '2': [1.0, 0.40000000596046448, 0.9585440361362102, 0.71692097187000003, 0.55846048891523226], '4': [9.9999999747524271e-07, 2.0, 1.0000004999999987]}, 'macro_absorption0': {'1': [20.0, 1800.0, 910.0], '0': [0.0, 1592.4991026562534, 80000.0, 45000.0, 20251.262658470834, 69748.737341529166], '3': [0.0, 1.7483529518358409e-05, 8.7417647591792047e-06], '2': [1.0, 0.40000000596046448, 0.85846048593500002, 0.55846048891523226, 0.71692097187000003], '4': [9.9999999747524271e-07, 2.0]}}

    FinalTuckerCoeffs =[2.96521471e+01, -2.08379223e-03,  1.83564712e-03, -9.03070210e-04
, -8.01356492e-04,  2.68801744e-03,  2.95618579e-03, -1.35250576e-03
, -2.15676440e-04, -1.07551343e-04, -5.48047832e-03,  2.72816542e-02
,  1.29386264e-03, -1.57009011e-04, -3.29310471e-05, -5.78211629e-03
,  2.05217576e-04,  8.38714251e-05,  2.37454258e-05,  2.68471342e-07
,  2.91809748e-03, -4.19445684e-04, -9.30790020e-04,  4.45878006e-04
,  3.26171254e-04, -3.02997835e-04,  1.16022408e-04, -1.82255117e-05
, -5.17929321e-06,  8.49713967e-06,  1.47282116e-03, -1.18473274e-02
,  1.79535724e-04, -4.56286079e-04,  2.21658431e-04,  1.44201590e-03
, -3.69470742e-04, -5.18859190e-05,  2.14868048e-05,  4.34348708e-06
, -9.37090792e-04, -2.32490647e-03,  6.31757240e-04,  6.08121585e-05
, -9.98714045e-05,  4.63936149e-05, -5.01465505e-05,  1.91014499e-05
, -4.73767100e-06, -3.94090583e-06, -3.23102066e-05, -1.07973322e-03
,  3.98905901e-04,  3.01200005e-05, -6.47591652e-05,  2.25533827e-05
, -1.82869368e-05,  8.46894690e-07,  1.79170985e-06,  4.66935000e-07
,  9.05637758e-04,  8.76721372e-05, -6.41062589e-05, -9.96532569e-05
,  2.19836235e-05,  6.63603065e-05, -2.21853509e-05, -4.72259449e-06
,  8.22079088e-07,  1.51424242e-06, -6.50831470e-04,  7.96849455e-05
,  1.09415123e-05,  4.76313303e-06, -2.23606028e-06, -1.33402156e-05
,  9.54722829e-06, -5.36563064e-06,  3.02656024e-06,  1.61617510e-06
, -2.87087782e-04,  3.17177534e-05,  2.19793307e-05, -2.87719033e-06
, -5.41261344e-06, -1.33836070e-06, -9.38778349e-07,  7.93158823e-07
, -4.88163948e-07, -2.80198557e-07,  1.58852413e-03, -5.58856348e-03
, -1.83567192e-04,  2.90575908e-04,  1.21254175e-04,  4.72775421e-04
, -1.70451000e-05, -2.74763771e-05, -8.24325443e-07,  1.90294650e-06
, -1.34976573e-03, -9.16572213e-04,  1.10641949e-04,  8.63818714e-05
, -1.70624492e-05,  2.17183980e-05, -1.19961142e-05,  7.50990343e-07
, -1.57166410e-06, -5.43114681e-07, -1.05115485e-04, -2.94945106e-04
,  1.06300597e-04,  3.50611044e-05, -1.64352230e-05, -3.23820242e-06
, -6.43347969e-06, -9.14619241e-07,  5.97459952e-07,  2.03462548e-07
, -2.65421884e-03, -6.23868430e-05,  1.41459378e-05, -1.80156977e-05
,  3.98297082e-06,  2.77554213e-05, -6.65290594e-06, -2.58139220e-06
, -2.49064048e-07,  7.86335655e-07, -4.01406233e-04,  5.84969366e-05
,  5.25797088e-06,  5.98900290e-06, -8.84444418e-07, -1.12674607e-05
,  7.20460524e-06, -4.08706743e-06,  2.13656894e-06,  1.18212208e-06
, -1.55857158e-04,  3.45160299e-05,  1.39653214e-05, -2.70505613e-07
, -3.01181853e-06, -3.15897688e-06, -1.05987542e-06,  5.75458961e-07
, -4.36970402e-07, -3.20938376e-07, -2.51363731e-04,  1.28198717e-05
,  2.85199964e-06, -3.63982293e-06, -4.78337315e-07,  2.17624679e-06
, -2.87795910e-07, -1.23156414e-06,  2.19873459e-07,  1.62344246e-07
, -2.85099397e-05, -5.93959470e-07,  6.77453298e-06, -2.50843448e-06
, -1.60543513e-06,  1.32026564e-06, -2.60458466e-06,  1.57377172e-06
, -8.36202546e-07, -5.44275766e-07, -1.77053922e-05,  5.20650878e-06
,  8.17131900e-07,  3.67089270e-08, -3.84129520e-07, -1.04848340e-06
,  6.34492207e-07, -4.63601441e-07,  1.43109863e-07, -5.13624562e-08]
    
    
    
    listOfFinalCoefIndexes_arr = [[0, 0, 0, 0, 0], [1, 0, 0, 0, 0], [2, 0, 0, 0, 0], [3, 0, 0, 0, 0], [4, 0, 0, 0, 0], [0, 1, 0, 0, 0], [1, 1, 0, 0, 0], [2, 1, 0, 0, 0], [3, 1, 0, 0, 0], [4, 1, 0, 0, 0], [0, 0, 1, 0, 0], [1, 0, 1, 0, 0], [2, 0, 1, 0, 0], [3, 0, 1, 0, 0], [4, 0, 1, 0, 0], [0, 1, 1, 0, 0], [1, 1, 1, 0, 0], [2, 1, 1, 0, 0], [3, 1, 1, 0, 0], [4, 1, 1, 0, 0], [0, 0, 2, 0, 0], [1, 0, 2, 0, 0], [2, 0, 2, 0, 0], [3, 0, 2, 0, 0], [4, 0, 2, 0, 0], [0, 1, 2, 0, 0], [1, 1, 2, 0, 0], [2, 1, 2, 0, 0], [3, 1, 2, 0, 0], [4, 1, 2, 0, 0], [0, 0, 0, 1, 0], [1, 0, 0, 1, 0], [2, 0, 0, 1, 0], [3, 0, 0, 1, 0], [4, 0, 0, 1, 0], [0, 1, 0, 1, 0], [1, 1, 0, 1, 0], [2, 1, 0, 1, 0], [3, 1, 0, 1, 0], [4, 1, 0, 1, 0], [0, 0, 1, 1, 0], [1, 0, 1, 1, 0], [2, 0, 1, 1, 0], [3, 0, 1, 1, 0], [4, 0, 1, 1, 0], [0, 1, 1, 1, 0], [1, 1, 1, 1, 0], [2, 1, 1, 1, 0], [3, 1, 1, 1, 0], [4, 1, 1, 1, 0], [0, 0, 2, 1, 0], [1, 0, 2, 1, 0], [2, 0, 2, 1, 0], [3, 0, 2, 1, 0], [4, 0, 2, 1, 0], [0, 1, 2, 1, 0], [1, 1, 2, 1, 0], [2, 1, 2, 1, 0], [3, 1, 2, 1, 0], [4, 1, 2, 1, 0], [0, 0, 0, 2, 0], [1, 0, 0, 2, 0], [2, 0, 0, 2, 0], [3, 0, 0, 2, 0], [4, 0, 0, 2, 0], [0, 1, 0, 2, 0], [1, 1, 0, 2, 0], [2, 1, 0, 2, 0], [3, 1, 0, 2, 0], [4, 1, 0, 2, 0], [0, 0, 1, 2, 0], [1, 0, 1, 2, 0], [2, 0, 1, 2, 0], [3, 0, 1, 2, 0], [4, 0, 1, 2, 0], [0, 1, 1, 2, 0], [1, 1, 1, 2, 0], [2, 1, 1, 2, 0], [3, 1, 1, 2, 0], [4, 1, 1, 2, 0], [0, 0, 2, 2, 0], [1, 0, 2, 2, 0], [2, 0, 2, 2, 0], [3, 0, 2, 2, 0], [4, 0, 2, 2, 0], [0, 1, 2, 2, 0], [1, 1, 2, 2, 0], [2, 1, 2, 2, 0], [3, 1, 2, 2, 0], [4, 1, 2, 2, 0], [0, 0, 0, 0, 1], [1, 0, 0, 0, 1], [2, 0, 0, 0, 1], [3, 0, 0, 0, 1], [4, 0, 0, 0, 1], [0, 1, 0, 0, 1], [1, 1, 0, 0, 1], [2, 1, 0, 0, 1], [3, 1, 0, 0, 1], [4, 1, 0, 0, 1], [0, 0, 1, 0, 1], [1, 0, 1, 0, 1], [2, 0, 1, 0, 1], [3, 0, 1, 0, 1], [4, 0, 1, 0, 1], [0, 1, 1, 0, 1], [1, 1, 1, 0, 1], [2, 1, 1, 0, 1], [3, 1, 1, 0, 1], [4, 1, 1, 0, 1], [0, 0, 2, 0, 1], [1, 0, 2, 0, 1], [2, 0, 2, 0, 1], [3, 0, 2, 0, 1], [4, 0, 2, 0, 1], [0, 1, 2, 0, 1], [1, 1, 2, 0, 1], [2, 1, 2, 0, 1], [3, 1, 2, 0, 1], [4, 1, 2, 0, 1], [0, 0, 0, 1, 1], [1, 0, 0, 1, 1], [2, 0, 0, 1, 1], [3, 0, 0, 1, 1], [4, 0, 0, 1, 1], [0, 1, 0, 1, 1], [1, 1, 0, 1, 1], [2, 1, 0, 1, 1], [3, 1, 0, 1, 1], [4, 1, 0, 1, 1], [0, 0, 1, 1, 1], [1, 0, 1, 1, 1], [2, 0, 1, 1, 1], [3, 0, 1, 1, 1], [4, 0, 1, 1, 1], [0, 1, 1, 1, 1], [1, 1, 1, 1, 1], [2, 1, 1, 1, 1], [3, 1, 1, 1, 1], [4, 1, 1, 1, 1], [0, 0, 2, 1, 1], [1, 0, 2, 1, 1], [2, 0, 2, 1, 1], [3, 0, 2, 1, 1], [4, 0, 2, 1, 1], [0, 1, 2, 1, 1], [1, 1, 2, 1, 1], [2, 1, 2, 1, 1], [3, 1, 2, 1, 1], [4, 1, 2, 1, 1], [0, 0, 0, 2, 1], [1, 0, 0, 2, 1], [2, 0, 0, 2, 1], [3, 0, 0, 2, 1], [4, 0, 0, 2, 1], [0, 1, 0, 2, 1], [1, 1, 0, 2, 1], [2, 1, 0, 2, 1], [3, 1, 0, 2, 1], [4, 1, 0, 2, 1], [0, 0, 1, 2, 1], [1, 0, 1, 2, 1], [2, 0, 1, 2, 1], [3, 0, 1, 2, 1], [4, 0, 1, 2, 1], [0, 1, 1, 2, 1], [1, 1, 1, 2, 1], [2, 1, 1, 2, 1], [3, 1, 1, 2, 1], [4, 1, 1, 2, 1], [0, 0, 2, 2, 1], [1, 0, 2, 2, 1], [2, 0, 2, 2, 1], [3, 0, 2, 2, 1], [4, 0, 2, 2, 1], [0, 1, 2, 2, 1], [1, 1, 2, 2, 1], [2, 1, 2, 2, 1], [3, 1, 2, 2, 1], [4, 1, 2, 2, 1]]
    
    ###===================================================
    #point = (150.0, 280.0, 0.6,  8.74176476e-06, 1.0)
    #point = (24000.0,	20.0,	0.40000000596,	1.74835295184e-05,	0.75)
    point = (17000.0,	1800.0,	0.620000004768,	1.74835295184e-05,	0.75)
    value_eva = evaluate(listOfBasicFctsUsingLagrangeInterpolation, listOfFinalCoefIndexes_arr, FinalTuckerCoeffs, point) 
    print value_eva
    
    #24000.0	20.0	0.40000000596	1.74835295184e-05	0.75	0.380883157253	0.38087658725	0.380885006086
    
    #17000.0	1800.0	0.620000004768	1.74835295184e-05	0.75	0.492738544941	0.492710966701	0.492733514241	-3.96303467188	-0.722919167469
    
    
    #print "==================="
    #print "listOfCrossSectionNamesWithArgs = ", listOfCrossSectionNamesWithArgs
    #print "==================="
    #####===================================================
    #criterion_post = "posteriori"
    #eps = 1E-6
    ####===================================================
    #### Déterminer toutes les décomposition de Tucker pour toutes les sections dans "listOfCrossSectionNamesWithArgs"
    ####===================================================
    #listOfTuckerDecomp = {}
    #listOfRefParameterCoordinates = {}
    #for csName in listOfCrossSectionNamesWithArgs:    
        ## les 2 ligne suivantes sont indispensables (calcul et stockage des coefficients)
        #listOfTuckerDecomp[csName] = Tucker.Tucker(listOfBasisFcts[csName], listOfNbOfBasisFcts[csName], listOfRefCrossSectionsToDetermineCoeffs[csName],listOfInterpolationPoints[csName], csName)
        #listOfTuckerDecomp[csName].getCoefficients()
        ## les lignes suivantes c'est pour le creusage 
        ##listOfTuckerDecomp[csName].getSparseCoefficients_posteriori (criterion_post, eps)
        ##listOfTuckerDecomp[csName].getSparseCoefficients_priori_index()
        ##listOfTuckerDecomp[csName].getIntersection_DirectEliminatedPoints_TensorStructurePoints()
    ####===================================================
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