# -*- coding: utf-8 -*-
#import methodsIntegralMultiDimens_TuckerHackbusch_bis as M
import multiDimensIntegralMethods as M
import multiDimensionOrthoNormalL2 as orth
import greedyAlgo as greedy
import numpy as np
from numpy import linalg as LA
import scipy
#from scipy import interpolate
import warnings
import LagrangePolynomial

class KL(M.domainD) :
    ###=====================================
    ### This class is a sub-class of the class "multiDimensIntegralMethods"
    ### It is associated with the KL decomposition to solve the eigenproblem: A*v=lambda*v 
    ### It is used to provide basis functions from the discretization of a Tucker grid (defined in "multiDimensIntegralMethods") :
    ### + basis functions for only one direction indicated here : Axis_k (There are "d" axes = "d" parameters")
    ### + For a given cross-section, for example : total0
    ###=====================================
    def __init__(self, Axis_k, listOfDomainBorders, listOfNbPointsOnEachIntervalForFinalDiscretization_Axis_k, \
    listOfTuckerGridNodes_Axis_k, listOfIntegralMethods_Axis_k, crossSection, crossSectionName): 
      #"""          
          #Argumentes:
            #"Axis_k"  [entier] : represente le paramètre qui est finement discretizé.
             #"listOfBorders" [dictionnaire] clé [chaine], valeur [liste soit une liste simple de bornes [min, max], soit une liste de listes (decoupa l'axe)]
             #"listOfNbPointsOnEachIntervalForFinalDiscretization_Axis_k" [liste de listes d'entiers]: ce sont les nombres de points dans la discrétisation
             #"listOfTuckerGridNodes_Axis_k" [un dictionnaire]: 
             #+ clé : c'est la chain de caractère correspondant à un numéro d'axe  ex '0', '1' etc
             #+ valeur : les coordonnées de points dans la discrétisation (pour une grille de Tucker) pour tous les axes
             #+ crossSection : lue à partir dklib
             #+ crossSectionName : une chain de caractère précise la section, groupe d'énergy, ect. ex : macro_totale0
          #Retourn:
            #rien (constructeur)
      #""" 
    
        M.domainD.__init__(self, Axis_k, listOfDomainBorders, listOfNbPointsOnEachIntervalForFinalDiscretization_Axis_k, listOfTuckerGridNodes_Axis_k) 
        M.domainD.discretizeDomainEntire(self, listOfIntegralMethods_Axis_k)        
      
        self.crossSection =  crossSection      
        self.crossSectionName = crossSectionName
  
    
##=====================================
        
    def getIntegralOperator_MatrixForm_Axis_k(self, MkI, WkI):  
        """
        focntion qui calcule la matrice A dans la décomposition de KL A*v = lambda*v, pour un seul axe "Axis_k"
        
        Argument:
             MkI [matrice] : valeur de section efficace sur la grille de Tucker
            WkI [matrice] : poids d'intégration 
        retourne:
           M_Axis_k [matrice] : la matrice dans le problème de valeurs propres pour la décomposition de KL
        """
    
        MW= MkI*WkI       
        MW_T =MW.T       
        M_Axis_k = np.dot(MkI,MW_T)    
        
        return M_Axis_k

##=====================================

    def eigsSort(self, M):
        """
      
        Argument:
             M [matrice] 
        retourne:
        une liste de :
           sortedEigVals [un array] : un array des valeurs propres en ordre croissant de matrice "M"
           sortedEigVects [un array de array] :  un array des vecteurs propres de "M" (en ordre croissant de valeurs propres)
        """
        eigs = LA.eig(M) # contain eigenvalues [0] (increasing order) and eigenvectors [1] 
        eigVals = eigs[0] 
        eigVects = eigs[1].T

        sortedEigVals = np.sort(eigVals )
        argSort = np.argsort(eigVals)
        sortedEigVects =  eigVects[argSort]       
        return [sortedEigVals, sortedEigVects]
##=====================================
    def finalEigens_Axis_k(self, M):
        """
        Fonction qui retourne les valeurs propres et vecteurs propres gardés selon un critère de convergence
        Argument:
             M [matrice] 
        retourne:
        une liste de :
           finalEigVals [un array] : un array des valeurs propres en ordre décroissant de matrice "M"
           finalEigVects [un array de array] :  un array des vecteurs propres de "M" (en ordre décroissant de valeurs propres)
        ATTENTION :  A changer ici la valeur de eps (coir corps du texte) pour le critère de sélection des modes propres
        """
        eigs = self.eigsSort(M)
          
        eigenVal = eigs[0] # increasing order eigenVal[0]<eigenVal[1]<....
        eigenVect = eigs[1]


        decreasedEigenVals = eigenVal[::-1] # to get decreasing order
        decreasedEigenVects = eigenVect[::-1]
        
        
        ###=====================================
        ### Criterion to select "j" eigenvectors via their corresponding eigenvalues : lambda_i/lambda_0 >= eps
        ###=====================================  
    
        n = len(eigenVal)       
        j = 0      
        
        ## !!! CECI EST LE CRITERE DE CONVERGENCE
        ## IL DEVRAIT etre en argumetn de la fonction au même titre que M
        eps = np.power(10.0,-9) ### UOX/MOX: eps = np.power(10.0,-10); UOX-Gd np.power(10.0,-9)
        for i in range(0,n):
            if decreasedEigenVals[i]/decreasedEigenVals[0] >= eps: 
                j = j+1
            else :
                break
       
        ###==========================
    
        
        
        finalEigVals =   decreasedEigenVals[:j] #decreasedEigenVals[:j] give:  decreasedEigenVals[0], ...,decreasedEigenVals[j-1]
        finalEigVects = decreasedEigenVects[:j]  
        ##===================================
              
        """
        For Paris 6
        """
        nameFile_infor_cs = "TuckerInfor_" +   self.crossSectionName + ".txt"  
        File_infor_cs = open(nameFile_infor_cs,'a')
        File_infor_cs.write("For cross section = %s\n" %(self.crossSectionName))
        
        """
        For Paris 6        """          
        
        File_infor_cs.write("=================\n")
        File_infor_cs.write("%s. BEGIN for information of KL in the DIRECTION = %s\n" %(self.Axis_k, self.Axis_k))
        File_infor_cs.write("In KL, to select the first eigenvectors, eps = %s\n" %(eps))
        File_infor_cs.write("=================\n")
        ###==========================
        File_infor_cs.write("=================\n")
        File_infor_cs.write("Total number of eigenvalues: %s\n" %(len(decreasedEigenVals)))
        File_infor_cs.write("=================\n")
        File_infor_cs.write("Eigenvalues in decreasing values, decreasedEigenVals = : %s\n" %(decreasedEigenVals))

        File_infor_cs.write("=================\n")
        File_infor_cs.write("Eigenvectors corresponding to decreasing eigenvalues, decreasedEigenVects = : %s\n" %(decreasedEigenVects))
    
        File_infor_cs.write("=================\n")
        File_infor_cs.write("Number of kept eigenvectors: %s/%s\n" %(len(finalEigVals),len(decreasedEigenVals)))
     
        File_infor_cs.write("=================\n")
        File_infor_cs.write("Kepts eigenvectors, finalEigVects: %s\n" %(finalEigVects))
      
        File_infor_cs.write("=================\n")
        File_infor_cs.close()
        #print "=================\n"
        #print "BEGIN for information of KL in the direction = ", self.Axis_k
        #print "In KL, to select the first eigenvectors, eps = ", eps
        #print "=================\n"
        ####==========================
        #print "=================\n"
        #print "Eigenvalues in decreasing values:\n"
        #print decreasedEigenVals
        #print "=================\n"
        #print "Eigenvectors corresponding to decreasing eigenvalues:\n"
        #print decreasedEigenVects
        #print "=================\n"
        #print "Number of kept eigenvectors:\n"
        #print len(finalEigVals)
        #print "=================\n"
        #print "Kepts eigenvectors:\n"
        #print finalEigVects
        #print "=================\n"
        
        return finalEigVals,  finalEigVects
##=====================================   
    def finalOrthoNormalEigens_Axis_k(self):  
    
        """
        Fonction qui orthonormise les vecteurs propres gardés dans l'axis "Axis_k" et compte le nombre de ces fonctions propres
        Argument:
             rien
        retourne:
            rien
        """
        MkI = self.getMatrixOfCrossSectionMkI(self.crossSection)      
        WkI = self.getWeights() 
        M_Axis_k = self.getIntegralOperator_MatrixForm_Axis_k(MkI , WkI)        
        self.finalEigVals_Axis_k, self.finalEigVects_Axis_k  = self.finalEigens_Axis_k(M_Axis_k) 
        
        ort =  orth.orthoNormalL2(self.weights[self.Axis_k])      
       
        self.finalOrthNormalizedEigVects_Axis_k =  ort.orthonormalVecteurs(self.finalEigVects_Axis_k)
        self.nbOfFcts_Axis_k = len(self.finalOrthNormalizedEigVects_Axis_k)
        
        
        """
        For Paris 6
        """
        nameFile_infor_cs = "TuckerInfor_" +   self.crossSectionName + ".txt"  
        File_infor_cs = open(nameFile_infor_cs,'a')
        
        File_infor_cs.write("=================\n")
        File_infor_cs.write("Recall: infor of the fine discretization:\n")
        File_infor_cs.write("\n")
        File_infor_cs.write("listOfBorders, self.listOfBorders[str(self.Axis_k)] = %s\n" %( self.listOfBorders[str(self.Axis_k)]))
        File_infor_cs.write("\n")
        File_infor_cs.write("listOfPoints, self.listOfTuckerGridNodes_Axis_k[self.Axis_k] = %s\n"%( self.listOfTuckerGridNodes_Axis_k[self.Axis_k]))
        File_infor_cs.write("=================\n")
        File_infor_cs.write("Orthonormalized eigenvectors, self.finalOrthNormalizedEigVects_Axis_k: %s\n" %( self.finalOrthNormalizedEigVects_Axis_k))
   
        File_infor_cs.write("\n")
        File_infor_cs.write("END for information of KL in the direction = %s\n" %(self.Axis_k))
        File_infor_cs.write("\n")
        File_infor_cs.write("=================\n")
        
        
        
        #print "=================\n"
        #print "Recall: infor of the fine discretization:\n"
        #print "listOfBorders = ", self.listOfBorders[str(self.Axis_k)]
        #print "listOfPoints = ", self.listOfTuckerGridNodes_Axis_k[self.Axis_k]
        #print "=================\n"
        #print "Orthonormalized eigenvectors:"
        #print self.finalOrthNormalizedEigVects_Axis_k
        #print "\n"
        #print "END for information of KL in the direction = ", self.Axis_k
        #print "=================\n"
        
        #print "=================\n"
        #print "Recall: infor of the fine discretization:\n"
        #print "listOfBorders = ", self.listOfBorders[str(self.Axis_k)]
        #print "listOfPoints = ", self.listOfTuckerGridNodes_Axis_k[self.Axis_k]
        #print "=================\n"
        #print "Orthonormalized eigenvectors:"
        #print self.finalOrthNormalizedEigVects_Axis_k
        #print "\n"
        #print "END for information of KL in the direction = ", self.Axis_k
        #print "=================\n"
        #raw_input()
##=====================================  

###===================================== 
    def getListOfInterpolationFcts(self):   
      
        """
        Fonction qui reconstruit les vecteurs propres gardés dans l'axis "Axis_k" par l'interpolation de Lagrange :
        + Une l'interpolation de Lagrange si il n'y a pas de coupage
        +  Des interpolations de Lagrange par morceau si il y a de coupages
        Argument:
             rien
        retourne:
            rien
        """  
               
        self.listOfBasicFctsUsingLagrangeInterpolation_Axis_k = []
        nbExtremes =  len(self.listOfBorders[str(self.Axis_k)])
      
        for i in range(self.nbOfFcts_Axis_k):
            ### Define each sub-polynomial/sub-interval if the main interval of "Axis_k" is subdivided 
            if nbExtremes > 2 :
                polLagrange_ki = []
                
                for j in range(nbExtremes - 1) :
                    j1Arr = np.where(self.listOfTuckerGridNodes_Axis_k[self.Axis_k] == self.listOfBorders[str(self.Axis_k)][j])[0].tolist()
                    j2Arr = np.where(self.listOfTuckerGridNodes_Axis_k[self.Axis_k] == self.listOfBorders[str(self.Axis_k)][j+1])[0].tolist()
                    
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
                      
                   

 