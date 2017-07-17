# -*- coding: utf-8 -*-
import KLdecomposition as KL
import numpy as np
import greedyAlgo as greedy
import LagrangeInterpolation as LI
import crossSections

##=====================================
###  Next function is used to gather all basis functions of all cross-sections, the main goal is to
### + Determine all common points per axis for all basis functions by the greedy algorithm 
### + From these common points, determine all points used in the linear system for the coefficients of Tucker, again by the greedy algorithm 
##=====================================

def getAllBasisFctsForAllSections(listOfCrossSectionNamesWithArgs, listOfDomainBorders, listOfNbPointsOnEachIntervalForFinalDiscretization,\
            listOfTuckerGridNodesFromDKlibs, listOfMethodsIntegrals, listOfCrossSectionsForTuckerGridsFromDklibs):
    
    """
    Fonction qui détaille toutes les fonctions de bases et le cardinal de ces fonctions/axe pour toutes les sections, axe par axe.
    On regroupe tous les polynomes de Lagrange calculés à partir des modes propres mais on en les traite  pas encore.
    Arguments :
      - listOfCrossSectionNamesWithArgs [une liste des chaines de caractères] : nom de toutes les sections testées avec les groupes d'énergy, anisotropy ...
        
        e.x listOfCrossSectionNamesWithArgs = 
        ['macro_totale0', 'macro_totale1', 'macro_scattering000', 'macro_scattering001', 'macro_scattering010', 'macro_scattering011', 'macro_nu*fission0', 'macro_nu*fission1']

      -  listOfDomainBorders [dictionnaire]:      
                + clé [une chaine de caractères] : désigne l'axe. ex :'0', '1'
                + valeur [liste de flottants] : soit une liste simple de bornes [min, max], soit une liste de listes (decoupa l'axe)]
                e.x : listOfDomainBorders {'1': [20.0, 1800.0], '0': [0.0, 150, 10000, 80000.0]}
                
      - listOfNbPointsOnEachIntervalForFinalDiscretization[dictionnaire]:
                 + clé [une chaine de caractères] : désigne l'axe. ex :'0', '1'
                 + valeur [liste de listes d'entiers]: ce sont les nombres de points dans la discrétisation
        e.x {'1': [2, [10], 2], '0': [[3, 9, 11], 2, 2], '2': [2, 2, [10]]}
        
      -  listOfTuckerGridNodesFromDKlibs [dictionnaire]:
                + clé [chaine de caractères] : désigne l'axe. ex :'0', '1'
                + valeur [liste de "d" arrays ((d= dimension)]: chaque array contient les coordonnées dans la discrétisation d'un axe. 
                "d" arrays correspondent à la discrétisation axiale pour une grille de Tucker 
                               
          e.x : '4': [array([     0.,  80000.]), array([   20.,  1800.]), array([ 0.40000001,  1.        ]), array([  0.00000000e+00,   1.74835295e-05]), array([  9.99999997e-07,   2.92894072e-01,   1.00000050e+00,
         1.70710693e+00,   2.00000000e+00])]
         (pour la direction de Xenon, où on a 5 points CC sur l'axe xenon, l'autres axes ont 2 points)
          
      
      - listOfMethodsIntegrals [dictionnaire]: 
              + clé [chaine de caractères] : désigne l'axe. ex :'0', '1'
              + valeur [une liste de "d" composants, chaque composant est une chaine de caractères] : désigne quelle méthode 
              d'intégration utilisée pour quel axe.
         e.x: {'1': ['T', 'CC', 'T', 'T', 'T'], '0': ['CCbis', 'T', 'T', 'T', 'T'], ...}
         
      - listOfCrossSectionsForTuckerGridsFromDklibs [dictionnaire]:
              + clé [chaine de caractères] : désigne le nom de section efficace. e.x 'macro_totale1'
              + valeur [dictionnaire] :
                + clé [chaine de caractères] : désigne l'axe. ex :'0', '1'
                + valeur [la même section mais lue à partir de "d" dklibs sur "d" grilles de Tucker]
                
              e.x : 'macro_totale1': {'1': <dkmodel.TFParameterXsLinear_Ptr; proxy of <Swig Object of type 'std::shared_ptr< Descartes::Tools::TFParameterXsLinear > *' at 0x7f6b7b51e4e0> >, 
              '0': <dkmodel.TFParameterXsLinear_Ptr; proxy of <Swig Object of type 'std::shared_ptr< Descartes::Tools::TFParameterXsLinear > *' at 0x7f6b7b51e450> >, ...}
     Retour :
      - listOfBasisFcts  [dictionnaire]: 
        + clé [chaine de caractères] : désigne le nom de section efficace. e.x 'macro_totale1'
        + valeur[dictionnaire]:
          + clé [chaine de caractères] : désigne l'axe. ex :'0', '1'
          + valeur [une liste des objets de la classe LagrangePolynomial] : chaque objet est une fonction de base = une polynome de Lagrange
            ou polynome de Lagrange par morceau (si le coupage). Le nombre de composant dans la liste est égale à nombre de mode propre 
            dans une direction pour la section dans la clé
            
        e.x : 
        
          'macro_totale1': {'1': [[<LagrangePolynomial.LagrangePolynomial instance at 0x7f6b7ae45710>], [<LagrangePolynomial.LagrangePolynomial instance at 0x7f6b7ae456c8>]],
          '0': [[<LagrangePolynomial.LagrangePolynomial instance at 0x7f6b7ae453f8>, <LagrangePolynomial.LagrangePolynomial instance at 0x7f6b7ae45440>, <LagrangePolynomial.LagrangePolynomial instance at 0x7f6b7ae45488>], 
          [<LagrangePolynomial.LagrangePolynomial instance at 0x7f6b7ae454d0>, <LagrangePolynomial.LagrangePolynomial instance at 0x7f6b7ae45518>, <LagrangePolynomial.LagrangePolynomial instance at 0x7f6b7ae45560>], 
          [<LagrangePolynomial.LagrangePolynomial instance at 0x7f6b7ae455a8>, <LagrangePolynomial.LagrangePolynomial instance at 0x7f6b7ae455f0>, <LagrangePolynomial.LagrangePolynomial instance at 0x7f6b7ae45638>]], 
          , ...}
          
          
      - listOfNbOfBasisFcts [dictionnaire]: 
        + clé [chaine de caractères] : désigne le nom de section efficace. e.x 'macro_totale1'
        + valeur [liste de "d" entiers]: le nombre de modes propres gardé par direction. 
            e.x [3, 2, 4, 2, 2] 
    """
    listOfBasisFcts = {}
    listOfNbOfBasisFcts = {}
    dimension = len(listOfDomainBorders)
    for csName in listOfCrossSectionNamesWithArgs:

      
        listOfBasisFcts_csName = {}
        listOfNbOfBasisFcts_csName = []
        for Axis_k in range(dimension):
            KLDecomposition_csName_Axis_k = KL.KL(Axis_k, listOfDomainBorders, listOfNbPointsOnEachIntervalForFinalDiscretization[str(Axis_k)],\
            listOfTuckerGridNodesFromDKlibs[str(Axis_k)], listOfMethodsIntegrals[str(Axis_k)], listOfCrossSectionsForTuckerGridsFromDklibs[csName][str(Axis_k)], csName)
            ### Pour appeler la fct qui fournit toutes les fcts de base (orthonormée) dans la direction k
            KLDecomposition_csName_Axis_k.finalOrthoNormalEigens_Axis_k()
            KLDecomposition_csName_Axis_k.getListOfInterpolationFcts()
            
            listOfBasisFcts_csName[str(Axis_k)] = KLDecomposition_csName_Axis_k.listOfBasicFctsUsingLagrangeInterpolation_Axis_k
            listOfNbOfBasisFcts_csName.append(KLDecomposition_csName_Axis_k.nbOfFcts_Axis_k)

            
        listOfBasisFcts[csName] =  listOfBasisFcts_csName
        listOfNbOfBasisFcts[csName] = listOfNbOfBasisFcts_csName

        CardinalOfBasisFcts = 1
        for Axis_k in range(dimension):
            CardinalOfBasisFcts = CardinalOfBasisFcts*listOfNbOfBasisFcts[csName][Axis_k]  
   
    return listOfBasisFcts,  listOfNbOfBasisFcts
    
###===================================================    
def getListOfBasisFctsByAxis(dimension, listOfCrossSectionNamesWithArgs, listOfBasisFcts,  listOfNbOfBasisFcts):
    """
  Fonction qui 
    + collecte et mélange toutes les fonctions de base, axe par axe pour la méthode de Greey après
    + définit le max du nombre de fonctions de base par axe sur toutes les sections efficaces -> pour déterminer le
      nombre de point commun entre toutes les sections
  Argument :
    - dimension [un entier] : le nombre de paramètre CRN
    - listOfCrossSectionNamesWithArgs [une liste des chaines de caractères] : nom de toutes les sections testées avec les groupes d'énergy, anisotropy ...
        
        e.x listOfCrossSectionNamesWithArgs = 
        ['macro_totale0', 'macro_totale1', 'macro_scattering000', 'macro_scattering001', 'macro_scattering010', 'macro_scattering011']

    - listOfBasisFcts  [dictionnaire]: 
    + clé [chaine de caractères] : désigne le nom de section efficace. e.x 'macro_totale1'
    + valeur[dictionnaire]:
      + clé [chaine de caractères] : désigne l'axe. ex :'0', '1'
      + valeur [une liste des objets de la classe LagrangePolynomial] : chaque objet est une fonction de base = une polynome de Lagrange
        ou polynome de Lagrange par morceau (si le coupage). Le nombre de composant dans la liste est égale à nombre de mode propre 
        dans une direction pour la section dans la clé
        
    e.x : 
    
      'macro_totale1': {'1': [[<LagrangePolynomial.LagrangePolynomial instance at 0x7f6b7ae45710>], [<LagrangePolynomial.LagrangePolynomial instance at 0x7f6b7ae456c8>]],
      '0': [[<LagrangePolynomial.LagrangePolynomial instance at 0x7f6b7ae453f8>, <LagrangePolynomial.LagrangePolynomial instance at 0x7f6b7ae45440>, <LagrangePolynomial.LagrangePolynomial instance at 0x7f6b7ae45488>], 
      [<LagrangePolynomial.LagrangePolynomial instance at 0x7f6b7ae454d0>, <LagrangePolynomial.LagrangePolynomial instance at 0x7f6b7ae45518>, <LagrangePolynomial.LagrangePolynomial instance at 0x7f6b7ae45560>], 
      [<LagrangePolynomial.LagrangePolynomial instance at 0x7f6b7ae455a8>, <LagrangePolynomial.LagrangePolynomial instance at 0x7f6b7ae455f0>, <LagrangePolynomial.LagrangePolynomial instance at 0x7f6b7ae45638>]], 
      , ...}
      
      
    - listOfNbOfBasisFcts [dictionnaire]: 
      + clé [chaine de caractères] : désigne le nom de section efficace. e.x 'macro_totale1'
      + valeur [liste de "d" entiers]: le nombre de modes propres respectivement gardé pour "d" directions 
          e.x [3, 2, 4, 2, 2] 
          
  Retour :
    - listOfMaxNbFcts [dictionnaire]:
       + clé [chaine de caractères] : désigne l'axe. ex :'0', '1'
       + valeur [liste de "d" entiers] : chaque entier = le max du nombre de fonctions de base par axe sur toutes les sections efficaces
    
    -  listOfAllBasisFcts [dictionnaire]: collecter toutes les fonctions de base/axe sur toutes les sections
       + clé [chaine de caractères] : désigne l'axe. ex :'0', '1'
       + valeur [liste des sous-listes, sous-liste contient une fonction de base/axe  = un polynome de Lagrange]
        e.x '1': [[<LagrangePolynomial.LagrangePolynomial instance at 0x7f21266d1950>], ...]
       
    """
    
   
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
###def getListOfPointsForGreedyAlgo(listOfDomainBorders, listOfNbPointsForGreedyAlgo):
  ###"""
  ###Fontion qui détermine/retrouve la liste de points dans les discrétisations fines des grille de Tucker, axe par axe.
  ###---> L'idée est de construire un échantillon des points sur le quel on sélecte les points par greedyAlgo
  
  
  ###"""
    ###listOfPointsForGreedyAlgo = {}    
   
    ###dimension = len(listOfDomainBorders)
    ###for Axis_k in range(dimension):
        ###nbExtremes = len(listOfDomainBorders[str(Axis_k)])

        ###listOfPointsForGreedyAlgo_Axis_k = np.asarray([])
        
        ###for indexOfBorder in range(nbExtremes-1):
            ###leftBorder = listOfDomainBorders[str(Axis_k)][indexOfBorder]
            ###rightBorder = listOfDomainBorders[str(Axis_k)][indexOfBorder + 1]
            ###listOfPointsForGreedyAlgo_Axis_k_OnSubInterval =  np.linspace(leftBorder,rightBorder,num=listOfNbPointsForGreedyAlgo[Axis_k])
            ###if indexOfBorder > 0 :
                ###listOfPointsForGreedyAlgo_Axis_k_OnSubInterval = np.delete(listOfPointsForGreedyAlgo_Axis_k_OnSubInterval, 0)
            ###listOfPointsForGreedyAlgo_Axis_k = np.concatenate((listOfPointsForGreedyAlgo_Axis_k, listOfPointsForGreedyAlgo_Axis_k_OnSubInterval), axis=0)
        
        ###listOfPointsForGreedyAlgo[str(Axis_k)] =  listOfPointsForGreedyAlgo_Axis_k
        
    ###print "========================="      
    ###print "listOfPointsForGreedyAlgo", listOfPointsForGreedyAlgo   
    ###print "========================="  
    ####raw_input()
    ###return listOfPointsForGreedyAlgo

###===================================================  
### Construct a list of points/Axis_k, used for all basis functions in this Axis_k 
###=================================================== 
def getFinalEmpiricalInterpolationPoints(dimension, listOfMaxNbFcts, listOfAllBasisFcts, listOfPointsForGreedyAlgo):
    listOfBasisFctOnPointsForGreedyAlgo = {}
    """
    Cette fonction retourne les points calculés par la méthode de greeedy (ici on parle plutot de  Empirical Interpolation Method)
    Ces poihnts seront utilisés pour la résoltution des coeff de Tucker (partie "b" dans Aa = b 
     où A est la matrice obtenue par les modes propres en des points choisis par greedy,, "a" est le vecteur des coefficients de  la decomposition de Tucker)
    Arguments :
    - dimension [un entier] : le nombre de paramètre CRN
      
      - listOfMaxNbFcts [dictionnaire]:
       + clé [chaine de caractères] : désigne l'axe. ex :'0', '1'
       + valeur [liste de "d" entiers] : chaque entier = le max du nombre de fonctions de base par axe sur toutes les sections efficaces
      
      -  listOfAllBasisFcts [dictionnaire]: collecter toutes les fonctions de base/axe sur toutes les sections
       + clé [chaine de caractères] : désigne l'axe. ex :'0', '1'
       + valeur [liste des sous-listes, sous-liste contient une fonction de base/axe  = un polynome de Lagrange]
        e.x '1': [[<LagrangePolynomial.LagrangePolynomial instance at 0x7f21266d1950>], ...]
        
      - listOfPointsForGreedyAlgo [dictionnaire]:
        + clé [chaine de caractères] : désigne l'axe. ex :'0', '1'
        + valeur [array des flottants] : valeurs de CC points sur la discrétisation fine d'un axes
          (contenant tous les extremes communs entre les segments dans le découpage)
          
        e.x:  
        '4': array([  9.99999997e-07,   2.92894072e-01,   1.00000050e+00,
         1.70710693e+00,   2.00000000e+00])
         
         '2': array([ 0.40000001,  0.44641201,  0.55846049,  0.67050897,  0.71692097,
        0.71692097,  0.75837694,  0.85846049,  0.95854404,  1.        ])
        
    Retours :
        - listOfFinalEmpiricalInterpolationPoints[dictionnaire]:
          + clé [chaine de caractères] : désigne l'axe. ex :'0', '1'
          + valeur [liste de flottants] : les points communs choisit par greedy (axe par axe) pour constituer les points dans "b"
           la résoultion des coefficients de Tucker "a" : Aa=b. 
           Size(liste) = le  max du nombre de fonctions de base sur un axe dans "clé" et sur toutes les sections 
           
        e.x : '1': [20.0, 1800.0, 280.67496474397274, 910.0]
    
    """
    
    
    ###ici , listOfPointsForGreedyAlgo = listOfTuckerGridNodes
    
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
    """
    Fonction qui vérifie si points choisits par le greedy sont bien ceux calculés dans la bibiothèque sortis de GAB
    (car sinon, des valeurs de "b" dans Aa =b ne sont que des valeurs approximatives)
    
    Argument :
      - dimension [un entier] : le nombre de paramètre CRN
      
      - interpolationPoints [dictionnaire]: : c'est exacte "listOfFinalEmpiricalInterpolationPoints" (voir précédence)
      
      - TuckerGridNodes [dictionnaire]:
        + clé [chaine de caractères] : désigne l'axe. ex :'0', '1'
        + valeur [array des flottants] : valeurs de CC points sur la discrétisation fine d'un axes
          (contenant tous les extremes communs entre les segments dans le découpage)
        e.x:  
        '4': array([  9.99999997e-07,   2.92894072e-01,   1.00000050e+00,
         1.70710693e+00,   2.00000000e+00])
         
         '2': array([ 0.40000001,  0.44641201,  0.55846049,  0.67050897,  0.71692097,
        0.71692097,  0.75837694,  0.85846049,  0.95854404,  1.        ])  
        
       - referencePoints [dictionnaire]: valeurs lues de dklib qui calcule les valeurs des sections dans "b" pour résoudre les coefficients "a" dans Aa =b
        +  clé [chaine de caractères] : désigne l'axe. ex :'0', '1'
        + valeur [liste de flottants] : valeurs lues d'un dklip, en général, elles ont des valeurs nominals 
          e.x : '1': [20.0, 280.6749572753906, 550.0, 910.0, 1539.324951171875, 1800.0], 
                '0': [0.0, 1592.489990234375, 5075.0, 20251.259765625, 45000.0, 69748.0, 80000.0]
    Retour:
      rien (il y a des avertisements)
    """
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
    """
     From determined greedy points for each direction, each cross-section chooses from these points for its own points    
     Argument:
      - Voir les descritptions précédentes
     Retour :
      - "listOfInterpolationPoints" [dictionnaire]:
        + clé[une chaine de caractères] : le nom de section (avec groupe d'énergy, anisotropy...)
        + valeur[dictionnaire]:
          + clé [chaine de caractères] : désigne l'axe. ex :'0', '1'
          + valeur[une liste de flottants] : les points CC choisits pour la section considérée 
          
        e.x : 'macro_fission1': {'1': [20.0, 1800.0, 910.0], '0': [0.0, 80000.0, 5075.0, 45000.0], ...}
    
    """
    
    
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
 

    
