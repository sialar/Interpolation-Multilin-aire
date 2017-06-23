# -*- coding: utf-8 -*-
import multiDimensIntegralMethods as M
import multiDimensionOrthoNormalL2 as orth
import greedyAlgo as greedy
import numpy as np
from numpy import linalg as LA
import scipy
import warnings
import LagrangePolynomial
##=====================================
###=====================================
### This class is used to provide basis functions:
### + In a given parameter/direction : Axis_k
### + For a given cross-section, for example : total0
###=====================================
### INPUTS:
### + Discretizations on Axis_k
### + Cross-section values (computed by GAB/AP2)
###=====================================
### METHODS:
### Solve the intergral equation of the given cross-section on the given direction ( Axis_k)
###=====================================
### OUTPUTS:
### Basis functions (Lagrange Polynomial(s)) for Axis_k and used in the Tucker decomposition
###=====================================
### TO CHANGE (following the test case):
### eps value in method : "def finalEigens_Axis_k(self, M)"
###=====================================


##=====================================
def separator():
    print '='*206

debug = 0
##=====================================

class KL(M.domainD) :
    def __init__(self, Axis_k, listOfDomainBorders, listOfNbPointsOnEachIntervalForFinalDiscretization_Axis_k, \
    listOfTuckerGridNodes_Axis_k, listOfIntegralMethods_Axis_k, crossSection, crossSectionName): # HERE : crossSection = f, crossSectionName = fNom

        M.domainD.__init__(self, Axis_k, listOfDomainBorders, listOfNbPointsOnEachIntervalForFinalDiscretization_Axis_k, listOfTuckerGridNodes_Axis_k)
        M.domainD.discretizeDomainEntire(self, listOfIntegralMethods_Axis_k)

        self.crossSection =  crossSection
        self.crossSectionName = crossSectionName

##=====================================

    def getIntegralOperator_MatrixForm_Axis_k(self, MkI, WkI):

        MW= MkI*WkI
        MW_T =MW.T
        M_Axis_k = np.dot(MkI,MW_T)

        return M_Axis_k

##=====================================

    def eigsSort(self, M):
        eigs = LA.eig(M) # contain eigenvalues [0] (increasing order) and eigenvectors [1]
        eigVals = eigs[0]
        eigVects = eigs[1].T

        sortedEigVals = np.sort(eigVals )
        argSort = np.argsort(eigVals)
        sortedEigVects =  eigVects[argSort]
        return [sortedEigVals, sortedEigVects]
##=====================================
    def finalEigens_Axis_k(self, M):
        eigs = self.eigsSort(M)

        eigenVal = eigs[0] # increasing order eigenVal[0]<eigenVal[1]<....
        eigenVect = eigs[1]


        decreasedEigenVals = eigenVal[::-1] # to get decreasing order
        decreasedEigenVects = eigenVect[::-1]

        if debug:
            separator()
            print "direction =", self.Axis_k
            print "decreasedEigenVals", decreasedEigenVals
            separator()
        #raw_input()
        ###=====================================
        ### Criterion to select "j" eigenvectors via their corresponding eigenvalues : lambda_i/lambda_0 >= eps
        ###=====================================
        n = len(eigenVal)
        j = 0
        #eps = np.power(10.0,-17)
        eps = np.power(10.0,-10)
        for i in range(0,n):
            if decreasedEigenVals[i]/decreasedEigenVals[0] >= eps:
                j = j+1
            else :
                break
        ###==========================
        if debug:
            separator()
            print "in KL, eps = ", eps
            separator()
        ###==========================
        #j = 5
        #if self.Axis_k == 2:
	    #j = 5

        #if self.Axis_k == 2:
	    #j = 4
        #if self.Axis_k == 1 | self.Axis_k == 2 :
            #j=5
        ###==========================
        finalEigVals =   decreasedEigenVals[:j] #decreasedEigenVals[:j] give:  decreasedEigenVals[0], ...,decreasedEigenVals[j-1]
        finalEigVects = decreasedEigenVects[:j]
        ##===================================
        if debug:
            separator()
            print "finalEigVals", finalEigVals
            #raw_input()
            separator()
        return finalEigVals,  finalEigVects
##=====================================
    def finalOrthoNormalEigens_Axis_k(self):

        MkI = self.getMatrixOfCrossSectionMkI(self.crossSection)
        WkI = self.getWeights()
        M_Axis_k = self.getIntegralOperator_MatrixForm_Axis_k(MkI , WkI)
        self.finalEigVals_Axis_k, self.finalEigVects_Axis_k  = self.finalEigens_Axis_k(M_Axis_k)

        ort =  orth.orthoNormalL2(self.weights[self.Axis_k])

        self.finalOrthNormalizedEigVects_Axis_k =  ort.orthonormalVecteurs(self.finalEigVects_Axis_k)
        self.nbOfFcts_Axis_k = len(self.finalOrthNormalizedEigVects_Axis_k)
##=====================================

###=====================================
    def getListOfInterpolationFcts(self):

        self.listOfBasicFctsUsingLagrangeInterpolation_Axis_k = []
        nbExtremes =  len(self.listOfBorders[str(self.Axis_k)])

        for i in range(self.nbOfFcts_Axis_k):
            ### Define each sub-polynomial/sub-interval if the main interval of Axis_k is subdivided
            if nbExtremes > 2 :
                polLagrange_ki = []

                for j in range(nbExtremes - 1) :
                    j1Arr = np.where(self.listOfTuckerGridNodes_Axis_k[self.Axis_k] == self.listOfBorders[str(self.Axis_k)][j])[0].tolist()
                    j2Arr = np.where(self.listOfTuckerGridNodes_Axis_k[self.Axis_k] == self.listOfBorders[str(self.Axis_k)][j+1])[0].tolist()
                    if debug:
                        print "j1=", j1Arr
                        print "j2=", j2Arr
                        print self.listOfBorders[str(self.Axis_k)]
                    j1 = j1Arr[-1]
                    j2 = j2Arr[0]
                    j2 = j2 + 1
                    ###because [j1:j2] take the values before j2

                    interp_k_couplage =  LagrangePolynomial.LagrangePolynomial(self.listOfTuckerGridNodes_Axis_k[self.Axis_k][j1:j2],\
                                              self.finalOrthNormalizedEigVects_Axis_k[i][j1:j2].T)

                    polLagrange_ki.append(interp_k_couplage)

            ### Define the polynomial if the main interval of Axis_k is NOT subdivided
            elif  nbExtremes == 2:

                polLagrange_ki = [ LagrangePolynomial.LagrangePolynomial(self.listOfTuckerGridNodes_Axis_k[self.Axis_k],\
                                        self.finalOrthNormalizedEigVects_Axis_k[i].T)]

            self.listOfBasicFctsUsingLagrangeInterpolation_Axis_k.append(polLagrange_ki)
