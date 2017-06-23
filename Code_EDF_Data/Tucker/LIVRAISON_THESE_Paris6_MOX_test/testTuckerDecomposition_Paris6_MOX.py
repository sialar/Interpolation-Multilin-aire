# -*- coding: utf-8 -*-
import numpy as np
import sys
from numpy import linalg as LA
import json, ast ### To convert string to dictionnary
###===================================================
### Nouveaux programmes:
import LagrangeInterpolation as LI
import LagrangePolynomial
import string
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

###===================================================
def find_nearest(myList,value):
    return min(myList, key=lambda x:abs(x-value)) 

###=================================================== 
def getListOfInterpolationFcts(listOfDomainBorders, listOfTuckerGridNodes, finalOrthNormalizedEigVects, Axis_k):   
      
        """
        Fonction qui reconstruit les vecteurs propres gardés dans l'axis "Axis_k" par l'interpolation de Lagrange :
        + Une l'interpolation de Lagrange si il n'y a pas de coupage
        +  Des interpolations de Lagrange par morceau si il y a de coupages
        Argument:
            + listOfDomainBorders : liste des extrêmes (après avoir effectué le coupage)
            + listOfTuckerGridNodes : dictionnaire, contenant les points dans tous les discrétisations fines
            + finalOrthNormalizedEigVects : dictionnaire, contenant tous les vecteurs de base orthonormés par direction 
            Axis_k : axis considéré
        retourne:
            les fonctions de base (les polynômes de Lagrange) pour la direction considérée (Axis_k)
        """  
               
        listOfBasicFctsUsingLagrangeInterpolation_Axis_k = []
        nbExtremes =  len(listOfDomainBorders[str(Axis_k)])
        nbOfFcts_Axis_k = len(finalOrthNormalizedEigVects[str(Axis_k)])
        
        listOfTuckerGridNodes[str(Axis_k)] = np.asarray(listOfTuckerGridNodes[str(Axis_k)])
        for i in range(nbOfFcts_Axis_k):
            ### Define each sub-polynomial/sub-interval if the main interval of "Axis_k" is subdivided 
       
            if nbExtremes > 2 :
                polLagrange_ki = []
                
                for j in range(nbExtremes - 1) :
                    
                    #print listOfTuckerGridNodes[str(Axis_k)]
                    #print listOfDomainBorders[str(Axis_k)][j]
                    #raw_input()
                    value1 = find_nearest(listOfTuckerGridNodes[str(Axis_k)], listOfDomainBorders[str(Axis_k)][j])
                    value2 = find_nearest(listOfTuckerGridNodes[str(Axis_k)], listOfDomainBorders[str(Axis_k)][j+1])

                    j1Arr =np.where(listOfTuckerGridNodes[str(Axis_k)] == value1)[0].tolist()
                    
                    j2Arr = np.where(listOfTuckerGridNodes[str(Axis_k)] == value2)[0].tolist()
    
                    
                    j1 = j1Arr[-1]
                    j2 = j2Arr[0]
                    j2 = j2 + 1

                    interp_k_couplage =  LagrangePolynomial.LagrangePolynomial(listOfTuckerGridNodes[str(Axis_k)][j1:j2],\
                                              finalOrthNormalizedEigVects[str(Axis_k)][i][j1:j2])
                          
                    polLagrange_ki.append(interp_k_couplage)
                    
            ### Define the polynomial if the main interval of Axis_k is NOT subdivided         
            elif  nbExtremes == 2: 
                polLagrange_ki = [ LagrangePolynomial.LagrangePolynomial(listOfTuckerGridNodes[str(Axis_k)],\
                                        finalOrthNormalizedEigVects[str(Axis_k)][i])] 
                
            listOfBasicFctsUsingLagrangeInterpolation_Axis_k.append(polLagrange_ki)
    
        return listOfBasicFctsUsingLagrangeInterpolation_Axis_k
    
###===================================================     
def evaluate(listOfBasisFcts, listOfFinalCoefIndexes_arr, FinalTuckerCoeffs, dimension, point) :
        """
      Fonction qui évalue la décompostion de Tucker pour la section considérée en un point donné
      Argument :
        + listOfBasisFcts : dictionnaire, contenant tous les fcts de base de toutes les directions.
        + FinalTuckerCoeffs : liste des coefficients de Tucker
        + listOfFinalCoefIndexes_arr : liste des indices des fcts de base associées aux coefficients de Tucker.
        + point : [un list de "d" (=dimension) flottants] : point sur le quel on évalue
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
def check_string(strs_begin, strs_end, NameFile):
    """
    Chercher une seule fois les strings dans un fichier
    """
    with open(NameFile) as f:
        found = False
        found_strs = ""
        for line in f:             
            if (strs_begin in line) and (strs_end in line):                
                return line
            elif (strs_begin in line) and (strs_end not in line):
                found_strs = found_strs + line
                while strs_end not in line:                    
                    line = next(f)
                    found_strs = found_strs + line

                found = True     
                return found_strs
        if found == False :
            print('The translation cannot be found!')

###=================================================== 
def get_list_check_string(strs_begin, strs_end, NameFile):
    """
    Chercher plusieurs fois les strings dans un fichier
    
    """
    list_found_strs = []
    with open(NameFile) as f:
        for line in f: 
            found_strs = ""
         
            if (strs_begin in line) and (strs_end not in line):
                found_strs = found_strs + line
                while strs_end not in line:                    
                    line = next(f)
                    found_strs = found_strs + line
                list_found_strs.append(found_strs)               
    return list_found_strs
###===================================================   

def convert_multiLines_oneLine(lines):
    converted_line =  ' '.join(' '.join(line.split()) for line in lines.split('\n')) #re.split(r'(\s+)', line)
    return converted_line
###===================================================  
def replace_str(strs, str_old, str_new):
    if str_old in strs:        
        new_strs = string.replace(strs, str_old, str_new)
        #print "new_strs", new_strs
        return new_strs
    else:
        warnings.warn("warning",'strs can not replaced')
"""
malformed string

https://stackoverflow.com/questions/14236708/how-to-convert-a-malformed-string-to-a-dictionary
best reponse:

https://stackoverflow.com/questions/32695699/valueerror-malformed-string-using-ast-literal-eval 
"""
###=================================================== 
def convert_str_dic(strs, strs_split):
    """
    Convertir strings à dictionnaire (ou liste) avec le split = strs_split
    """
    if strs_split in strs:
        strs_dic = strs.split(strs_split)[1]   ### 2 parties, on récupére la 2ème     
        dic = ast.literal_eval(strs_dic)
        return dic
    else:
        warnings.warn("warning",'strs_split is not in strs')
###===================================================  

###===================================================
def getFinalOrthNormalizedEigVects(NameFile):
    strs_begin = "Orthonormalized eigenvectors"
    strs_end = "]]"   
    
    list_check_string = get_list_check_string(strs_begin, strs_end, NameFile)
    
    #print "======================"
    #print list_check_string
    #print "len(list_check_string)", len(list_check_string)
    #print "======================"
    
    finalOrthNormalizedEigVects = {}
    
    for i in range(dimension) :
        lines = list_check_string[i]
        #print "==========lines==========="
        #print lines
        #print "====================="
        #raw_input()
        converted_line = convert_multiLines_oneLine(lines)
       
        
        str_old = "[ " ### there are many cases : "[[ ", [[-"
        str_new = "["
        if str_old in converted_line:
            converted_line = replace_str(converted_line, str_old, str_new)
        
        str_old = "]] "
        str_new = "]]"
        converted_line = replace_str(converted_line, str_old, str_new)
        
        str_old = " "
        str_new = ", "
        converted_line = replace_str(converted_line, str_old, str_new)
        #print "==========converted_line==========="
        #print converted_line
        #print "====================="
        #raw_input()
        #print "converted_line\n", converted_line
        #raw_input
        strs_split = "self.finalOrthNormalizedEigVects_Axis_k:, "
        finalOrthNormalizedEigVects[str(i)] = convert_str_dic(converted_line, strs_split)
    
    #print "========finalOrthNormalizedEigVects============"
    #print finalOrthNormalizedEigVects
    #raw_input()
    return finalOrthNormalizedEigVects
###===================================================
def getFinalTuckerCoeffs(NameFile):
    strs_begin = "Coefficient values" #self.finalOrthNormalizedEigVects: [[ "
    strs_end = "]"
    lines = check_string(strs_begin, strs_end, NameFile)
    
    converted_line = convert_multiLines_oneLine(lines)
    
    str_old = "[ " ### there are many cases : "[[ ", [[-"
    str_new = "["
    if str_old in converted_line:
        converted_line = replace_str(converted_line, str_old, str_new)
        
    str_old = "]: " ### there are many cases : "[[ ", [[-"
    str_new = "]"
    if str_old in converted_line:
        converted_line = replace_str(converted_line, str_old, str_new)    
    
    str_old = " "
    str_new = ", "
    converted_line = replace_str(converted_line, str_old, str_new)
    
    #print "==========================="
    #print "converted_line", converted_line
    #print "==========================="
    #raw_input()
    strs_split = "self.FinalTuckerCoeffs, ="
    FinalTuckerCoeffs = convert_str_dic(converted_line, strs_split)
    return FinalTuckerCoeffs

###===================================================
def getListOfFinalCoefIndexes_arr(NameFile):
    strs_begin = "self.listOfFinalCoefIndexes_arr"
    strs_end = "]]"
    lines = check_string(strs_begin, strs_end, NameFile)
    
    converted_line = convert_multiLines_oneLine(lines)
    
    #print "==========================="
    #print "converted_line", converted_line
    #print "==========================="
    #raw_input()

    strs_split = "self.listOfFinalCoefIndexes_arr = "
    listOfFinalCoefIndexes_arr = convert_str_dic(converted_line, strs_split)
    return listOfFinalCoefIndexes_arr
###===================================================
def getListOfBasicFcts(listOfDomainBorders, listOfTuckerGridNodes, finalOrthNormalizedEigVects, dimension ):

    listOfBasicFctsUsingLagrangeInterpolation = {}
    for Axis_k in range(dimension):
        basisFcts_Axis_k = getListOfInterpolationFcts(listOfDomainBorders, listOfTuckerGridNodes, finalOrthNormalizedEigVects, Axis_k)
        
        listOfBasicFctsUsingLagrangeInterpolation[str(Axis_k)] = basisFcts_Axis_k
    return listOfBasicFctsUsingLagrangeInterpolation
###===================================================


"""
Les fichiers à lire :

GeneralInfor_UOXGd.txt (ou GeneralInfor_UOX.txt, GeneralInfor_MOX.txt)

TuckerInfor_macro_totale1.txt
TuckerInfor_macro_totale0.txt
TuckerInfor_macro_scattering011.txt
TuckerInfor_macro_scattering010.txt
TuckerInfor_macro_scattering001.txt
TuckerInfor_macro_scattering000.txt
TuckerInfor_macro_nu_fission1.txt
TuckerInfor_macro_nu_fission0.txt
TuckerInfor_macro_fission1.txt
TuckerInfor_macro_fission0.txt
TuckerInfor_macro_absorption1.txt
TuckerInfor_macro_absorption0.txt

"""
if __name__ == "__main__":
    
    dimension = 5
   
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
    """
    BEGIN : récupérer des données généralees pour toutes les sections efficaces
    """
    ###===================================================
    """
    - listOfDomainBorders : list des extrêmes après la définition de coupage
    - Ici, listOfDomainBorders est lu à partir de fichier GeneralInfor_UOX.txt 
    - (3 cas tests : UOX, MOX, UOXGd mais ils ont la même discrétisation)
    
    """
    ###===================================================    
    """
    -Récupérer listOfDomainBorders de "GeneralInfor_UOX.txt " sous forme string puis transformer sous forme dictionnaire
    - Assigner à la variable listOfDomainBorders
    - listOfDomainBorders : list des extrêmes après la définition de coupage
    - Ici, listOfDomainBorders est lu à partir de fichier GeneralInfor_UOX.txt 
    - (3 cas tests : UOX, MOX, UOXGd mais ils ont la même discrétisation)
    """
    ###===================================================  
    """
    CHANGE ICI POUR CHAQUE CAS TEST : UOX, MOX, UOXGd
    """
    NameFile = "GeneralInfor_MOX.txt" 
    ###===================================================  
    
    strs_begin = "listOfDomainBorders = {"
    strs_end = "}"    
    
    lines = check_string(strs_begin, strs_end, NameFile)
    #line = check_string(strs, NameFile)
    converted_line = convert_multiLines_oneLine(lines)

    strs_split = " = "
    listOfDomainBorders = convert_str_dic(converted_line, strs_split)
 
    
    #print "listOfDomainBorders =", listOfDomainBorders
    #print "type(listOfDomainBorders) = ", type(listOfDomainBorders)
    
 ###=========================================   
    """
    -Récupérer listOfTuckerGridNodes de "GeneralInfor_UOX.txt " sous forme string puis transformer sous forme dictionnaire
    """
    strs_begin = "listOfTuckerGridNodes = {"
    strs_end = "}"

    lines = check_string(strs_begin, strs_end, NameFile)

    converted_line = convert_multiLines_oneLine(lines)
    """
    Eliminate all "array(" and  ")" - caused malformed strings
    """
    str_old = "array("
    str_new = ""
    converted_line = replace_str(converted_line, str_old, str_new)
    
    str_old = ")"
    str_new = ""
    converted_line = replace_str(converted_line, str_old, str_new)

    strs_split = " = "
    listOfTuckerGridNodes_dic = convert_str_dic( converted_line, strs_split)
    
    #print "============="
    #print "listOfTuckerGridNodes_dic = ", listOfTuckerGridNodes_dic    
    #print "type(listOfTuckerGridNodes_dic) = ", type(listOfTuckerGridNodes_dic)
    #print "============="
    #raw_input()
 
    listOfTuckerGridNodes = listOfTuckerGridNodes_dic 
    ###===================================================
    """
    END : récupérer des données généralees pour toutes les sections efficaces
    """
    ###===================================================
 ###========================================= 
    """
    Chaque fichier contient des informations de la décompostion de Tucker pour UNE section efficace 
    
    """
    listOfCrossSectionNamesWithArgs = ['macro_totale0', 'macro_totale1', 'macro_absorption0', 'macro_absorption1', 'macro_scattering000', 'macro_scattering001', 'macro_scattering010', 'macro_scattering011', 'macro_nu*fission0', 'macro_nu*fission1', 'macro_fission0', 'macro_fission1']

    listOfFiles = {'macro_totale1' : "TuckerInfor_macro_totale1.txt", 'macro_absorption0': "TuckerInfor_macro_absorption0.txt", 'macro_absorption1': "TuckerInfor_macro_absorption1.txt", 'macro_scattering000': "TuckerInfor_macro_scattering000.txt", 'macro_scattering001':"TuckerInfor_macro_scattering001.txt", 'macro_scattering010': "TuckerInfor_macro_scattering010.txt", 'macro_scattering011' :  "TuckerInfor_macro_scattering011.txt", 'macro_nu*fission0': "TuckerInfor_macro_nu_fission0.txt", 'macro_nu*fission1':"TuckerInfor_macro_nu_fission1.txt", 'macro_fission0': "TuckerInfor_macro_fission0.txt", 'macro_fission1': "TuckerInfor_macro_fission1.txt", 'macro_totale0': "TuckerInfor_macro_totale0.txt"}

###============================================
    """
   Définir tous les données de Tucker (fcts de base, coefficients) pour toutes les sections  
   """
###============================================
    finalOrthNormalizedEigVects = {}
    FinalTuckerCoeffs = {}
    listOfFinalCoefIndexes_arr = {}
    listOfBasicFctsUsingLagrangeInterpolation = {}
    
    for csName in listOfCrossSectionNamesWithArgs:
     
        NameFile = listOfFiles[csName]
        finalOrthNormalizedEigVects[csName]  = getFinalOrthNormalizedEigVects(NameFile)
        
        FinalTuckerCoeffs[csName]  = getFinalTuckerCoeffs(NameFile)
        
        listOfFinalCoefIndexes_arr[csName]  = getListOfFinalCoefIndexes_arr(NameFile)
        
        listOfBasicFctsUsingLagrangeInterpolation[csName]  = getListOfBasicFcts(listOfDomainBorders, listOfTuckerGridNodes, finalOrthNormalizedEigVects[csName], dimension)
    
###=========================================   

    """
    Récupérer les résultats de référence
    
    ordre bu tc rho b10 xe AP2 Cocagne Tucker
    => point = (bu tc rho b10 xe)
    """
    ###===================================================
    ###===================================================
    """
    CHANGE ICI POUR CHAQUE CAS TEST : UOX, MOX, UOXGd
    """
    repertoire = "FinalResults_eps1E_10_14000p_MOX_extendedD_Paris6/"
    nameTest="RelativeError_AP2_Cocagne_Tucker_withGreedy_eps1E_10_14000p_mox_"
    ###===================================================
    NameFileTest = repertoire + nameTest

    listOfRefResultsFiles = {}
    for csName in listOfCrossSectionNamesWithArgs:    
        listOfRefResultsFiles[csName] = NameFileTest + csName
        if csName == 'macro_nu*fission0':
            listOfRefResultsFiles[csName] = NameFileTest + 'macro_nu_fission0'
        if csName == 'macro_nu*fission1':
            listOfRefResultsFiles[csName] = NameFileTest + 'macro_nu_fission1'
    
    ###===================================================
    """
    Récupérer les points dans la grille de référence
    """
    ###===================================================
    listOfReferencePoints = []
    NameFile = listOfRefResultsFiles['macro_totale0']
    f=open(NameFile,"r")
    lines=f.readlines()

    for x in lines:
        point_str = []
        for i in range(1, dimension + 1): 
            point_str.append(float(x.split('\t')[i]))
        listOfReferencePoints.append(point_str)        
    f.close()
    
    ###===================================================
    """
    Récupérer les valeurs exactes (AP2) par section 
    """
    ###===================================================
    listOfAP2 = {}    
    for csName in listOfCrossSectionNamesWithArgs:
        NameFile = listOfRefResultsFiles[csName]
        f=open(NameFile,"r")
        lines=f.readlines()
        AP2_idx = 6
        AP2_csName = []
        for x in lines:
            AP2_csName.append(float(x.split('\t')[AP2_idx]))
          
        f.close()
        listOfAP2[csName] =  AP2_csName
        
        
    ###==================1. TESTER AVEC UNE SEULE SECTION=================================
    """
    -Si on veut tester avec une seule section, 
    - Sinon, mettre en commentaire 3 lignes :
    
    iSection = ...
    csName = listOfCrossSectionNamesWithArgs[iSection]
    listOfCrossSectionNamesWithArgs = [csName]
    
    Ici, on a :
        listOfCrossSectionNamesWithArgs = ['macro_totale0', 'macro_totale1', 'macro_absorption0', 'macro_absorption1', 'macro_scattering000', 'macro_scattering001', 'macro_scattering010', 'macro_scattering011', 'macro_nu*fission0', 'macro_nu*fission1', 'macro_fission0', 'macro_fission1']
    """
    iSection = 11
    csName = listOfCrossSectionNamesWithArgs[iSection]
    listOfCrossSectionNamesWithArgs = [csName]
    ###===================================================  
 
    ###==================2. TESTER AVEC LES SECTIONS DANS REACTIVITE=================================
    """
    Si on veut tester rapidement pour avoir la reactivité, change  listOfCrossSectionNamesWithArgs :
    Sinon, mettre en commentaire la ligne :
     listOfCrossSectionNamesWithArgs = ['macro_totale0', 'macro_totale1', 'macro_scattering000', 'macro_scattering001', 'macro_scattering010', 'macro_scattering011', 'macro_nu*fission0', 'macro_nu*fission1']
    """
    #listOfCrossSectionNamesWithArgs = ['macro_totale0', 'macro_totale1', 'macro_scattering000', 'macro_scattering001', 'macro_scattering010', 'macro_scattering011', 'macro_nu*fission0', 'macro_nu*fission1']
  
    ###=================================================== 
    
    ###==================3. TESTER AVEC TOUTES (=12) SECTIONS =================================
    """
    Si on veut tester pour avoir la reactivité et toutes les sections
    Sinon, mettre en commentaire la ligne :
   
    """
    #listOfCrossSectionNamesWithArgs = ['macro_totale0', 'macro_totale1', 'macro_absorption0', 'macro_absorption1', 'macro_scattering000', 'macro_scattering001', 'macro_scattering010', 'macro_scattering011', 'macro_nu*fission0', 'macro_nu*fission1', 'macro_fission0', 'macro_fission1']
  
    ###=================================================== 
 
    ###=====================================================================
    ### Calculer et mesurer erreur de l'évaluation
    ###=====================================================================
    listOfValues_Tucker = {}
    listOfValues_AP2= {}

    ###===================================================

    listOfMaxRelativeError_Tucker_kinf = {}
    listOfMSE_Tucker_kinf = {}   

    ###===================================================
            
    for csName in listOfCrossSectionNamesWithArgs:        
  
        count = 1
        values_Tucker = []        
        values_AP2 = listOfAP2[csName] 
        """
        point = (bu, tc, ct, density ,cb, xe)
        """
        for i in range(len(listOfReferencePoints)):
            point = listOfReferencePoints[i]
            value_T = evaluate(listOfBasicFctsUsingLagrangeInterpolation[csName], listOfFinalCoefIndexes_arr[csName], FinalTuckerCoeffs[csName], dimension, point) 
            value_A = values_AP2[i]
      

            values_Tucker.append(value_T) 
            #values_AP2.append(value_A)
  
            
            ####%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            print "count =", count
            print "point = ", point
            print "AP2 = %s, Tucker = %s" %(values_AP2[count-1], values_Tucker[count-1].real)
            print "Tucker - AP2 = ", values_Tucker[count-1] - values_AP2[count-1]
      
            print "=================="
            print "abs((Tucker - AP2)/AP2) = ", abs((values_Tucker[count-1] - values_AP2[count-1])/values_AP2[count-1])
            print"==================="
          
            raw_input()
            ####%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            count = count + 1
            
        listOfValues_Tucker[csName] = np.asarray(values_Tucker) 
        listOfValues_AP2[csName]= np.asarray(values_AP2) 

                      
    ###======================Relative error and RMS error=============================     
    listOfRelativeError_Tucker = {}    
    listOfMaxRelativeError_Tucker = {}     
    listOfMSE_Tucker = {}  

    ###======================Kinf============================= 
    for csName in listOfCrossSectionNamesWithArgs:
        ### 4/3/2016: normaliser e_relative = |f-f~|/maxf
        maxAP2 = max(abs(listOfValues_AP2[csName]))
        listOfRelativeError_Tucker[csName] = 1E5*(listOfValues_Tucker[csName] - listOfValues_AP2[csName])/maxAP2

        listOfMaxRelativeError_Tucker[csName] = max(abs(listOfRelativeError_Tucker[csName])) 
  
        
        AP2_RS = np.sqrt(sum(np.power(listOfValues_AP2[csName],2)))
        
        listOfMSE_Tucker[csName]  = np.sqrt(sum(np.power(abs(listOfValues_Tucker[csName] - listOfValues_AP2[csName])*1E5, 2)))/AP2_RS

        print "====================="
        print csName
        print "relativeError_Tucker", listOfMaxRelativeError_Tucker[csName]
        
        print "MSE_Tucker", listOfMSE_Tucker[csName] 
        print "=====================" 
        
        
    ###====================Reactivity ========================       
    
    kinf_Tucker =  computeKinf(listOfValues_Tucker)                           
    kinf_AP2 =  computeKinf(listOfValues_AP2)                               

    ###============================================   

    
    reactivityError_Tucker = (1.0/kinf_Tucker - 1.0/kinf_AP2)*1E5
    
    count = 1    
    sumReactivity_A = 0
    
    for point in listOfReferencePoints: 
        reactivity_A = kinf_AP2[count-1]
        reactivity_T = kinf_Tucker[count-1]
        sumReactivity_A = sumReactivity_A  + reactivity_A*reactivity_A 
        
        err_T = (1.0/reactivity_T - 1.0/reactivity_A)*1E5
        count = count +1        

    ###=====================   
    
    maxReactivityAbsError_Tucker = max(abs(reactivityError_Tucker))

    
    
    reactivity_L2 = np.sqrt(sumReactivity_A)    
    reactivityMSE_Tucker = np.sqrt((sum(reactivityError_Tucker*reactivityError_Tucker)))/reactivity_L2
 
    
    print "====================="
    print "len(reactivityError_Tucker)", len(reactivityError_Tucker)
 
    print "maxReactivityAbsError_Tucker = ", maxReactivityAbsError_Tucker
    print "maxReactivityAbsError_Tucker_priori = ", maxReactivityAbsError_Tucker_priori

    
    print "reactivityMSE_Tucker = ", reactivityMSE_Tucker
    
    ###============================================
