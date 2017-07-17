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
        listOfBasicFctsUsingLagrangeInterpolation_Axis_k = []
        nbExtremes =  len(listOfDomainBorders[str(Axis_k)])
        nbOfFcts_Axis_k = len(finalOrthNormalizedEigVects[str(Axis_k)])

        listOfTuckerGridNodes[str(Axis_k)] = np.asarray(listOfTuckerGridNodes[str(Axis_k)])
        for i in range(nbOfFcts_Axis_k):

            if nbExtremes > 2 :
                polLagrange_ki = []

                for j in range(nbExtremes - 1) :
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

            elif  nbExtremes == 2:
                polLagrange_ki = [ LagrangePolynomial.LagrangePolynomial(listOfTuckerGridNodes[str(Axis_k)],\
                                        finalOrthNormalizedEigVects[str(Axis_k)][i])]

            listOfBasicFctsUsingLagrangeInterpolation_Axis_k.append(polLagrange_ki)

        return listOfBasicFctsUsingLagrangeInterpolation_Axis_k

###===================================================
def evaluate(listOfBasisFcts, listOfFinalCoefIndexes_arr, FinalTuckerCoeffs, dimension, point) :
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

def convert_str_dic(strs, strs_split):
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
    finalOrthNormalizedEigVects = {}

    for i in range(dimension) :
        lines = list_check_string[i]
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
        strs_split = "self.finalOrthNormalizedEigVects_Axis_k:, "
        finalOrthNormalizedEigVects[str(i)] = convert_str_dic(converted_line, strs_split)
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
    strs_split = "self.FinalTuckerCoeffs, ="
    FinalTuckerCoeffs = convert_str_dic(converted_line, strs_split)
    return FinalTuckerCoeffs

###===================================================
def getListOfFinalCoefIndexes_arr(NameFile):
    strs_begin = "self.listOfFinalCoefIndexes_arr"
    strs_end = "]]"
    lines = check_string(strs_begin, strs_end, NameFile)

    converted_line = convert_multiLines_oneLine(lines)
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

    listOfCrossSectionNamesWithArgs = ['macro_totale0', 'macro_totale1', 'macro_absorption0', 'macro_absorption1', 'macro_scattering000', 'macro_scattering001', 'macro_scattering010', 'macro_scattering011', 'macro_nu*fission0', 'macro_nu*fission1', 'macro_fission0', 'macro_fission1']

    listOfSubdivisionDirection = {}
    listOfValuesForSubdivision = {}
    listOfNumberOfPointsForSubdivision = {}


    listOfSubdivisionDirection['0'] = 1
    listOfValuesForSubdivision['0'] = [150, 10000]
    listOfNumberOfPointsForSubdivision['0'] = [3, 9]  # 3 points entre [debut, 150], 9 points entre [150, 10000], le reste pour [10000, fin]. A priori c'est inutile et redondant


    listOfSubdivisionDirection['2'] = 1
    listOfValuesForSubdivision['2'] = [0.71692097187]
    listOfNumberOfPointsForSubdivision['2'] = [5]

    NameFile = "GeneralInfor_UOXGd.txt"
    strs_begin = "listOfDomainBorders = {"
    strs_end = "}"

    lines = check_string(strs_begin, strs_end, NameFile)
    converted_line = convert_multiLines_oneLine(lines)

    strs_split = " = "
    listOfDomainBorders = convert_str_dic(converted_line, strs_split)

    strs_begin = "listOfTuckerGridNodes = {"
    strs_end = "}"

    lines = check_string(strs_begin, strs_end, NameFile)
    converted_line = convert_multiLines_oneLine(lines)
    str_old = "array("
    str_new = ""
    converted_line = replace_str(converted_line, str_old, str_new)
    str_old = ")"
    str_new = ""
    converted_line = replace_str(converted_line, str_old, str_new)

    strs_split = " = "
    listOfTuckerGridNodes_dic = convert_str_dic( converted_line, strs_split)

    listOfTuckerGridNodes = listOfTuckerGridNodes_dic
    listOfCrossSectionNamesWithArgs = ['macro_totale0', 'macro_totale1', 'macro_absorption0', 'macro_absorption1', 'macro_scattering000', 'macro_scattering001', 'macro_scattering010', 'macro_scattering011', 'macro_nu*fission0', 'macro_nu*fission1', 'macro_fission0', 'macro_fission1']
    listOfFiles = {'macro_totale1' : "TuckerInfor_macro_totale1.txt", 'macro_absorption0': "TuckerInfor_macro_absorption0.txt", 'macro_absorption1': "TuckerInfor_macro_absorption1.txt", 'macro_scattering000': "TuckerInfor_macro_scattering000.txt", 'macro_scattering001':"TuckerInfor_macro_scattering001.txt", 'macro_scattering010': "TuckerInfor_macro_scattering010.txt", 'macro_scattering011' :  "TuckerInfor_macro_scattering011.txt", 'macro_nu*fission0': "TuckerInfor_macro_nu_fission0.txt", 'macro_nu*fission1':"TuckerInfor_macro_nu_fission1.txt", 'macro_fission0': "TuckerInfor_macro_fission0.txt", 'macro_fission1': "TuckerInfor_macro_fission1.txt", 'macro_totale0': "TuckerInfor_macro_totale0.txt"}
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


    #/**************************************************************************/
    csName = sys.argv[1]
    point = []
    for i in range(2,len(sys.argv)):
        point.append(float(sys.argv[i]))

    print evaluate(listOfBasicFctsUsingLagrangeInterpolation[csName], listOfFinalCoefIndexes_arr[csName], FinalTuckerCoeffs[csName], dimension, point)
