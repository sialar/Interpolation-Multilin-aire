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

def find_nearest(myList,value):
    return min(myList, key=lambda x:abs(x-value)) 

"""
Des changements dans cette fct "getListOfInterpolationFcts"
"""
def getListOfInterpolationFcts(listOfDomainBorders, listOfTuckerGridNodes, finalOrthNormalizedEigVects, Axis_k):   
      
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
def check_string(strs_begin, strs_end, NameFile):
    """
    Data which we want to recover from a file can be: a line or multilines 
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

###=================================================== 
def check_string_bis(strs_begin, strs_end, NameFile):
    """
    Data which we want to recover from a file can be: a line or multilines 
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
                return found_strs, line
        if found == False :
            print('The translation cannot be found!')
###=================================================== 


def get_list_check_string(strs_begin, strs_end, NameFile):
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
    if strs_split in strs:
        strs_dic = strs.split(strs_split)[1]        
        dic = ast.literal_eval(strs_dic)
        return dic
    else:
        warnings.warn("warning",'strs_split is not in strs')
###===================================================  
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
    - listOfDomainBorders : list des extrêmes après la définition de coupage
    - Ici, listOfDomainBorders est lu à partir de fichier GeneralInfor_MOX.txt 
    - (3 cas tests : UOX, MOX, UOXGd mais ils ont la même discrétisation)
    
    """
    ###===================================================    
    """
    -Récupérer listOfDomainBorders de "GeneralInfor_MOX.txt " sous forme string puis transformer sous forme dictionnaire
    - Assigner à la variable listOfDomainBorders
    """
    strs_begin = "listOfDomainBorders = {"
    strs_end = "}"
    NameFile = "GeneralInfor_MOX.txt" 
    lines = check_string(strs_begin, strs_end, NameFile)
    #line = check_string(strs, NameFile)
    converted_line = convert_multiLines_oneLine(lines)

    strs_split = " = "
    listOfDomainBorders = convert_str_dic(converted_line, strs_split)
 
    
    #print "listOfDomainBorders =", listOfDomainBorders
    #print "type(listOfDomainBorders) = ", type(listOfDomainBorders)
    
 ###=========================================   
 
     ###===================================================
    """
    - listOfDomainBorders : list des extrêmes après la définition de coupage
    - Ici, listOfDomainBorders est lu à partir de fichier GeneralInfor_MOX.txt 
    - (3 cas tests : UOX, MOX, UOXGd mais ils ont la même discrétisation)
    
    """
    ###===================================================
    
    """
    -Récupérer listOfDomainBorders de "GeneralInfor_MOX.txt " sous forme string puis transformer sous forme dictionnaire
    - Assigner à la variable listOfDomainBorders
    """
    strs_begin = "listOfTuckerGridNodes = {"
    strs_end = "}"
    NameFile = "GeneralInfor_MOX.txt" 
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
 ###=========================================   
    listOfFiles = {'0' : "TuckerInfor_macro_totale1.txt", '1': "TuckerInfor_macro_absorption0.txt", '2': "TuckerInfor_macro_absorption1.txt", '3': "TuckerInfor_macro_scattering000.txt", '4':"TuckerInfor_macro_scattering001.txt", '5': "TuckerInfor_macro_scattering010.txt", '6' :  "TuckerInfor_macro_scattering011.txt", '7': "TuckerInfor_macro_nu_fission0.txt", '8':"TuckerInfor_macro_nu_fission1.txt", '9': "TuckerInfor_macro_fission0.txt", '10': "TuckerInfor_macro_fission1.txt", '11': "TuckerInfor_macro_totale0.txt"}

###============================================
    strs_begin = "Orthonormalized eigenvectors" #self.finalOrthNormalizedEigVects: [[ "
    strs_end = "]]"
    NameFile = listOfFiles['11']
    
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
    
###======================================================    
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
    
    
###====================================================== 

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
    
###========================================= 
    listOfBasicFctsUsingLagrangeInterpolation = {}
    for Axis_k in range(dimension):
        basisFcts_Axis_k = getListOfInterpolationFcts(listOfDomainBorders, listOfTuckerGridNodes, finalOrthNormalizedEigVects, Axis_k)
        
        listOfBasicFctsUsingLagrangeInterpolation[str(Axis_k)] = basisFcts_Axis_k
    
    #print "listOfBasicFctsUsingLagrangeInterpolation", listOfBasicFctsUsingLagrangeInterpolation
    #raw_input()

###===================================================
    #point = (150.0, 280.0, 0.6,  8.74176476e-06, 1.0)
    #point = (24000.0,	20.0,	0.40000000596,	1.74835295184e-05,	0.75)
    point = (17000.0,	1800.0,	0.620000004768,	1.74835295184e-05,	0.75)
    value_eva = evaluate(listOfBasicFctsUsingLagrangeInterpolation, listOfFinalCoefIndexes_arr, FinalTuckerCoeffs, point) 
    print value_eva
    
   