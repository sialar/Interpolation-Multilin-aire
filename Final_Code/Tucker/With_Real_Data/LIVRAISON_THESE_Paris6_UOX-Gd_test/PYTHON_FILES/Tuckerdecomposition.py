# -*- coding: utf-8 -*-
import numpy as np
import KLdecomposition as KL
import LagrangeInterpolation as LI
import random

### INPUTS:
### + List of all basis functions for "d" directions (defined by KLdecomposition class)
### + Reference points (determined by Greedy algorithm) to solve system of coefficients
##=====================================
### OUTPUTS:
### + Coefficients of the Tucker decomposition of the given "crossSectionName"
### + Tucker decomposition which can evaluate any given point
##=====================================
###=====================================
### TO CHANGE (following the test case): 
### listOfKeptPercent for a priori sparsity
###=====================================

###########################################=======================
###  Formule pour trouver les coeeficients de la décomposition de Tucker:   somme (ai  produit Phij) = b  = Aa, A = produit Phij, a =ai 
############################################



class Tucker:
##=====================================
### This class is used to define the Tucker decomposition for a given cross-section (via crossSectionName, e.g: macro_totale0)
##=====================================
    def __init__(self, listOfBasisFcts, listOfNbOfBasisFcts, refCrossSectionToDetermineCoeffs, listOfInterpolationPoints, crossSectionName):
        """
        Data for the given cross-section via its "crossSectionName"
        "listOfBasisFcts" [un dictionnaire]: 
                  + clé [une chain de caractère d'un entier] : ex '0', '1' ---> désigne l'axe
                  + Valeur [une liste de sous-listes] : chaque sous-liste contient un objet correspondant à une fonction de base (une interpolation de Lagrange)
                  ou plusieurs objets correspondant à fonction de base par morceau si on fait le coupage 
                  e.x : [<LagrangePolynomial.LagrangePolynomial instance at 0x2eb3098>]
                  ou [<LagrangePolynomial.LagrangePolynomial instance at 0x2eb35a8>, <LagrangePolynomial.LagrangePolynomial instance at 0x2eb35f0>, <LagrangePolynomial.LagrangePolynomial instance at 0x2eb3638>]
        
        "listOfNbOfBasisFcts" [une liste de "d" entiers]: 
                  + chaque entier est le nombre de fonctions de base prises après la KL, e.x [7,3,5,2,3]
        "refCrossSectionToDetermineCoeffs" [une section lue de dkzip] : cette section calcule les points dans "b" du système des coefficients a: Aa = b
        
        "listOfInterpolationPoints" [un dictionnaire]: 
                  + clé [une chain de caractère d'un entier] : ex '0', '1' ---> désigne l'axe
                  + Valeur [une liste de flottants] : 
                      + liste des points choisis par axe en fonction de modes propres gardés (et par greedy) pour constituer les points utilisés dans "b"
                      + liste est utilisée pour une section déterminée par "crossSectionName"
                      e.x [0.0, 0.96984631039295421, 0.11697777844051105, 0.58682408883346515, 0.75] 
                      
        "crossSectionName" [une chaine de caractères] : désigne le nom de section efficace que on l'approach. ex :  "macro_totale0"                 
                  
        """
     

        
        self.listOfBasisFcts = listOfBasisFcts
        self.listOfNbOfBasisFcts  = listOfNbOfBasisFcts 
        self.refCrossSectionToDetermineCoeffs = refCrossSectionToDetermineCoeffs      
        self.listOfInterpolationPoints = listOfInterpolationPoints
        self.crossSectionName = crossSectionName
        
        """
          "dimension" [un entier] : le nombre de paramètre de CRN
          "CardinalOfBasisFcts" [un entier] : le produit de tous les nombres de modes propres gardés pour la section considérée ici.
                    e.x : nombre de modes propres par axe pour dimension = 5: [5,2,3,2,3] => CardinalOfBasisFcts = 5*2*3*2*3 = 180
        """
        self.dimension = len(self.listOfBasisFcts)
      
        ###=====================================    
        CardinalOfBasisFcts = 1
        for Axis_k in range(self.dimension):
            CardinalOfBasisFcts = CardinalOfBasisFcts*listOfNbOfBasisFcts[Axis_k]
            
        ###===================================== 
        self.CardinalOfBasisFcts = CardinalOfBasisFcts
 ##===================================== 
    def produitArr(self, u): 
        """
        rgument:
            "u" : [liste ou array de flottants]
            
        retourn:
            un flottant qui est le produit des elements de le tableau "u"
        """
        n = len(u)   
        s = 1.0
        for i in range(0,n):
            s = s*u[i]
        return s
##=====================================    
    def transfIndexIntToArray(self, I, d, listOfNbPoints): 
        """
        
        focntion qui transforme un indice entier à une liste I[entier] ---> [i_0, ...,i_(d-1)]
         Argument:
         I [un entier] : indice de points sur une grille en "d" dimensions
         d [un entier] : la dimension = le nombre de paramètres CRN
         listOfNbPoints [une liste de entiers] : la base pour transformer l'indice
        
            
        retourn:
        "i" : une liste de "d" entier qui désigne la transformation de I dans la base "listOfNbPoints"
            
        """ 
    
        I = int(I)           
        i = []
        ii = [0]
        for k in range(0,d):
            i.append(ii)
                
        index = range(0,d)
        indexInverse = index[::-1] 
            
        for j in indexInverse:
            if j > 0 :            
                i[j] = int(I/self.produitArr(listOfNbPoints[:j] ))
                I = I - i[j]*self.produitArr(listOfNbPoints[:j] )
                I = int(I)
            
            else :
                i[j] = int(I) 
            
        return i

 ##=====================================
    def lastSquare(self, A, b):
        """
       Argument: 
          A[une matrice] : taille : mxn 
          b[un array] : taille m (m>n)
      retourn:
          sol_n [un array] : solution "a" de Aa = b par le moindres carrées
        """
        sol_np = np.linalg.lstsq(A, b)[0]
        return sol_np
 ##=====================================
    def getCoefficients(self):
        """
      Fonction qui détermine :
          + "self.FinalTuckerCoeffs"[une liste de flottants] : les coefficients de la décomposition de Tucker pour la section considérée (via crossSectionName).
            
          + " self.listOfFinalCoefIndexes_arr" [une liste de liste de entiers]: 
                    chaque liste de entiers contient "d" entiers qui désigne les indices des modes propres associés à un coefficient. 
                    e.x a_I*phi_i1*phi_i2...*phi_i5, => I <=> [i1,i2, ...,i5]            
       
       Argument: 
          rien
      retourn:
          rien 
        """
        ### Matrice "A" et vecteur "b" dans Aa=b
        self.matrixA_full = np.zeros(self.CardinalOfBasisFcts*self.CardinalOfBasisFcts).reshape(self.CardinalOfBasisFcts,self.CardinalOfBasisFcts)
        self.vectorOfExactValues = np.zeros(self.CardinalOfBasisFcts) 
        ###==============================================
        ### self.listOfRefPoint used to determine cross-section values in b (on the right hand size of the system A*coeffs = b)
        self.listOfRefPoint = []
        ###==============================================
        for i in range(self.CardinalOfBasisFcts):
            ### ==================== Determine exact value :f(x1,x2,...xd)====================
            iArr = self.transfIndexIntToArray(i, self.dimension, self.listOfNbOfBasisFcts)           
            refPoint = []
           
            for Axis_k in range(self.dimension):
               
                refValue_Axis_k = self.listOfInterpolationPoints[str(Axis_k)][iArr[Axis_k]]
                refPoint.append(refValue_Axis_k)
            
            self.vectorOfExactValues[i]  = self.refCrossSectionToDetermineCoeffs.evaluate(refPoint)
            
            self.listOfRefPoint.append(refPoint)
           
            ### ========================================
            ### Determine the elements of matrix "A" via the evaluation of eigenvectors on reference points 
            ### ========================================
           
            for j in range(self.CardinalOfBasisFcts):
                p = 1.0
                jArr = self.transfIndexIntToArray(j, self.dimension, self.listOfNbOfBasisFcts)   
                for Axis_k in range(self.dimension):         
     
                    basisFct =  self.listOfBasisFcts[str(Axis_k)][jArr[Axis_k]]
                    valueOfEigenVector_Axis_k_jk = LI.getInterpolation(basisFct, refPoint[Axis_k])
                   
                    p = p*valueOfEigenVector_Axis_k_jk
                   
                self.matrixA_full[i][j] = p
        ###===================================================
        ### self.FinalTuckerCoeffs is an array, each index [I] associated with an array from indexes of basis functions: I = [i1,i2,...,id]
                                                                                                                                
        self.FinalTuckerCoeffs = np.linalg.solve(self.matrixA_full,self.vectorOfExactValues)
        ### Transform list in I to list in [i1,i2,...,id]
        self.listOfFinalCoefIndexes_arr = []
        for j in range(self.CardinalOfBasisFcts):
            jArr = self.transfIndexIntToArray(j, self.dimension, self.listOfNbOfBasisFcts)         
            self.listOfFinalCoefIndexes_arr.append(jArr)
        
        
        """
        For Paris 6
        """
        nameFile_infor_cs = "TuckerInfor_" +   self.crossSectionName + ".txt"  
        File_infor_cs = open(nameFile_infor_cs,'a')
        
        #nameFile_infor_cs = "TuckerInfor_" +   self.crossSectionName + ".txt"  
        #File_infor_cs = open(nameFile_infor_cs,'r+')
        
        File_infor_cs.write("For cross section = %s\n" %(self.crossSectionName))
        File_infor_cs.write("=================\n")
        File_infor_cs.write("BEGIN for information of coefficients in Tucker decomposition")
        File_infor_cs.write("\n")
        File_infor_cs.write("self.listOfNbOfBasisFcts = %s\n"%(self.listOfNbOfBasisFcts))
        File_infor_cs.write("\n")
        
        File_infor_cs.write("self.CardinalOfBasisFcts =%s\n" %(self.CardinalOfBasisFcts))
        File_infor_cs.write("\n")
        File_infor_cs.write("Number of coefficients, len(self.FinalTuckerCoeffs) = %s\n" %(len(self.FinalTuckerCoeffs)))
        File_infor_cs.write("=================\n")
        File_infor_cs.write("Coefficient values, self.FinalTuckerCoeffs =%s:\n" %(self.FinalTuckerCoeffs))
        File_infor_cs.write("\n")
        
        File_infor_cs.write("=================\n")
        File_infor_cs.write("Coefficient indexes under array format :\n")
        File_infor_cs.write("\n")
        File_infor_cs.write("self.listOfFinalCoefIndexes_arr = %s\n"%(self.listOfFinalCoefIndexes_arr))
        File_infor_cs.write("\n")
        File_infor_cs.write("END for information of coefficients in Tucker decomposition\n")
        File_infor_cs.write("=================\n")
        File_infor_cs.close()
        
        #print "For cross section = %s\n" %(self.crossSectionName)
        #raw_input()

###===================================================
    ###==================================================
    ### For a posteriori sparse representation
    ### Criterion : KEEP a_I if |a_I|/sum(|a_I|) > eps
    ### eps ~  1E-6
    ###===================================================
    def getSparseCoefficients_posteriori (self, criterion, eps): 
        """
      Fonction qui détermine les coefficients ainsi que les indices associés après un creusage a posteriori.
      Ici, le critère testé est : garder a_I si |a_I|/sum(|a_I|) > eps
      
      Les nouveaux coefficients et indices sont dans :
          + "self.FinalTuckerCoeffs_sparse_posteriori" [une liste de flottants] 
          + "self.listOfFinalCoefIndexes_arr_posteriori" [une liste de liste de entiers]
       Argument: 
            "criterion" [une chaine de caractères] 
            "eps" [un flottant]  
      retourn:
          rien 
        """

        keptIndexesAfterSparsity_posteriori = []
        if criterion == "posteriori":

            sumOfCoef = np.sum(abs(self.FinalTuckerCoeffs))
           
            for j in range(self.CardinalOfBasisFcts):                             
                if abs(self.FinalTuckerCoeffs[j])/sumOfCoef > eps :
                    keptIndexesAfterSparsity_posteriori.append(j)
        ### List of a posteriori array indexes             
        self.listOfFinalCoefIndexes_arr_posteriori = [self.listOfFinalCoefIndexes_arr[i] for i in keptIndexesAfterSparsity_posteriori]
        self.FinalTuckerCoeffs_sparse_posteriori = self.FinalTuckerCoeffs[keptIndexesAfterSparsity_posteriori]
        
        print "===========A POSTERIORI================"
        print "eps =", eps
        print "======================================="
        print "crossSectionName = ", self.crossSectionName
        print "listOfNbOfBasisFcts = ", self.listOfNbOfBasisFcts
        print "Total = ", self.CardinalOfBasisFcts
        reducedNb = len(self.listOfFinalCoefIndexes_arr) - len(self.listOfFinalCoefIndexes_arr_posteriori)
        print "Percentage of a A POSTERIORI reduction = %s %%" %(100.0*reducedNb/len(self.listOfFinalCoefIndexes_arr))
        print "=========END A POSTERIORI=============="
###=============================================================
    ###==================================================
    ### For a priori sparse representation
    ### Criterion : KEEP a_I if I = [i1,..,id] s.t: (i1*...*id)/(r1*...*rd) < eps <=1
    ### rk : number of basis function in the direction k
    ###===================================================
    def f_iArr_indexCriterion(self,iArr): ### an operator applied to iArr
        """
      Fonction qui effectue un opérateur sur l'indice sous forme array d'un coefficient.
      Opérateur est pour calculer les valeurs selon un critère de creusage
      Arguments: 
        iArr [une liste ou un array de entiers] : le multi-indice d'un coefficient provenant des indices des modes propres associés
      Retour:
        productOfIndexes [un entier] : le produit des composants dans "iArr"
      
        """
        productOfIndexes = 1.0
        for Axis_k in range(self.dimension):                   
            productOfIndexes = productOfIndexes*float((iArr[Axis_k]+1))/self.listOfNbOfBasisFcts[Axis_k]    
        return productOfIndexes
 ##===================================== 
  ##=====================================  
    def getKeptCoefIndexes_PrioriCriterion_Index(self,listOfFinalCoefIndexes, keptPercent):
        """
      Fonction qui détermine les INDICES des coefficients gardés en fonction de pourcentage proposé par un critère de creusage
      Arguments: 
        listOfFinalCoefIndexes [une liste de entiers] : les indices en entier des coefficients
        keptPercent [un flottant entre 0.0 et 1.0] : le pourcentage des coefficients gardés après un creusage
      Retour:
        keptIndexesAfterSparsity [une liste de entiers] : les indices en entier des coefficients après le creusage
      
        """
        list_f_iArr = [] 
        for i in range(len(listOfFinalCoefIndexes)):
            iArr = listOfFinalCoefIndexes[i]
            f_iArr = self.f_iArr_indexCriterion(iArr)
            list_f_iArr.append(f_iArr)
        
        arr_f_iArr = np.asarray(list_f_iArr)       
        #sort_list = np.sort(arr_f_iArr) # order croissan
        sort_indices =  np.argsort(arr_f_iArr) # order croissant, f_i small -> i important -> will be kept

        ### N : number of elements which are kept
        N = int(len(listOfFinalCoefIndexes)*keptPercent)
        ### Keep only the first N coefs
        keptIndexesAfterSparsity = sort_indices[:N]        
        return keptIndexesAfterSparsity
 ##=====================================  
   
    def getKeptRefPointIndexes_PrioriSparsity(self, keptPercent):
        """
        Fonction qui détermine les indices gardés  dans "b" du système Aa=b 
        (en fonction du pourcentage des coefficients gardé)
        Cette fonction est utilisée pour tester le creusage tensoriel         
        
        Arguments:           
          keptPercent [un flottant entre 0.0 et 1.0] : le pourcentage des coefficients gardés après un creusage
        Retour:
          keptRefPointIndexes_tensorStructure [une liste des entiers] : les indices gardés dans "b"
         
        
        """
        iArrMax = self.listOfNbOfBasisFcts 

        nbFct_k_max = np.max(iArrMax)
        Axis_k_max = np.argmax(iArrMax)
        
        #nbFct_k_max = np.min(iArrMax)
        #Axis_k_max = np.argmin(iArrMax)
        
        #nbFct_k_max = iArrMax[self.dimension-1]
        #Axis_k_max = self.dimension-1 #np.argmin(iArrMax)
        
        index_k_max_tensorStructure = nbFct_k_max - 1 ### In reality, index is scaled by 1 because it begins by 0
        
                  
        ### First test, reduce only 1 point on the axe having the maximal number of basis fcts, such that : R(Initial) > R' > R" (after the a priori sparsity)
        if nbFct_k_max*keptPercent >= 1.0: ### Subtract by 1 but still having basis functions in this direction k
            index_k_max_tensorStructure = index_k_max_tensorStructure - 1 ###  R' can be defined
            print "================="
            print "Case argmax"
            print "nbFct =", nbFct_k_max
            print "Axis =", Axis_k_max
            print "================="
            #raw_input()
        else:
            ###  R' can not be defined, R'=R
            print "================="
            print "Can not reduce by tensor structure"
            print "================="
            #raw_input()
        ### List of integer indexes
        keptRefPointIndexes_tensorStructure = []
        
        for i in range(self.CardinalOfBasisFcts):
            iArr = self.transfIndexIntToArray(i, self.dimension, self.listOfNbOfBasisFcts)
            if iArr[Axis_k_max] <= index_k_max_tensorStructure:
                keptRefPointIndexes_tensorStructure.append(i)   

        return keptRefPointIndexes_tensorStructure

###=============================================================
    def getSparseCoefficients_priori_index (self):
        """
        Fonction qui détermine les coefficients gardés (a-> a') en fonction de pourcentage proposé par un critère de creusage.
        Cette fonction détermine aussi les indices gardés associés avec les coefficients gardés
        Après, des technique de creusage dans "b" sont proposées pour obtenir b', ainsi que la résolution du syxtème A'a' = b':
          + Si b'=b  : utiliser le moindres carrées pour A'a'=b'
          + Si size(a') < size(b') < size(b) et b' est obtenue par un creusage tensoriel: utiliser le moindres carrées pour A'a'=b'
          + Si size(a') = size(b') : résoudre directement  A'a'=b'
        
      Arguments: 
        rien
      Retour:
        rien
      
        """
        """
        Les dictionnaires suivants proposent un pourcentage de creusage des coefficients "a"
        A MODIFIER POUR CHAQUE CAS 
      
        """
        # There are percentages of kept coefs after a priori sparsenesse
        #E.g listOfKeptPercentPriori_IndexCriterion['macro_totale0'] = 0.60 <=> 60% of coefs are kept <=> Reduced 40% of coefs
        ###listOfKeptPercentPriori_IndexCriterion = {}
        ###listOfKeptPercentPriori_IndexCriterion['macro_totale0'] = 0.60 
        ###listOfKeptPercentPriori_IndexCriterion['macro_totale1'] = 0.50 
        ###listOfKeptPercentPriori_IndexCriterion['macro_absorption0']= 0.4 
        ###listOfKeptPercentPriori_IndexCriterion['macro_absorption1']= 0.60 
        ###listOfKeptPercentPriori_IndexCriterion['macro_scattering000'] = 0.60  
        ###listOfKeptPercentPriori_IndexCriterion['macro_scattering001']= 0.4 
        ###listOfKeptPercentPriori_IndexCriterion['macro_scattering010']= 0.4
        ###listOfKeptPercentPriori_IndexCriterion['macro_scattering011']= 0.4 
        ###listOfKeptPercentPriori_IndexCriterion['macro_nu*fission0']= 0.4 
        ###listOfKeptPercentPriori_IndexCriterion['macro_nu*fission1']= 0.4
        ###listOfKeptPercentPriori_IndexCriterion['macro_fission0']= 0.4
        ###listOfKeptPercentPriori_IndexCriterion['macro_fission1']= 0.4 
        
        ###=======================UOX cas for direct sparsity ======================================
        #listOfKeptPercentPriori_IndexCriterion = {}
        #listOfKeptPercentPriori_IndexCriterion['macro_totale0'] = 0.70 
        #listOfKeptPercentPriori_IndexCriterion['macro_totale1'] = 0.70 
        #listOfKeptPercentPriori_IndexCriterion['macro_absorption0']= 0.7 
        #listOfKeptPercentPriori_IndexCriterion['macro_absorption1']= 0.60 
        #listOfKeptPercentPriori_IndexCriterion['macro_scattering000'] = 0.7  
        #listOfKeptPercentPriori_IndexCriterion['macro_scattering001']= 0.4 
        #listOfKeptPercentPriori_IndexCriterion['macro_scattering010']= 0.7
        #listOfKeptPercentPriori_IndexCriterion['macro_scattering011']= 0.7 
        #listOfKeptPercentPriori_IndexCriterion['macro_nu*fission0']= 0.4 
        #listOfKeptPercentPriori_IndexCriterion['macro_nu*fission1']= 0.4
        #listOfKeptPercentPriori_IndexCriterion['macro_fission0']= 0.4
        #listOfKeptPercentPriori_IndexCriterion['macro_fission1']= 0.4
       
        ### =======================MOX cas for direct sparsity ========================================
        
        #listOfKeptPercentPriori_IndexCriterion = {}
        #listOfKeptPercentPriori_IndexCriterion['macro_totale0'] = 0.70 
        #listOfKeptPercentPriori_IndexCriterion['macro_totale1'] = 0.70 
        #listOfKeptPercentPriori_IndexCriterion['macro_absorption0']= 0.70 
        #listOfKeptPercentPriori_IndexCriterion['macro_absorption1']= 0.70 
        #listOfKeptPercentPriori_IndexCriterion['macro_scattering000'] = 0.70  
        #listOfKeptPercentPriori_IndexCriterion['macro_scattering001']= 0.70 
        #listOfKeptPercentPriori_IndexCriterion['macro_scattering010']= 0.70
        #listOfKeptPercentPriori_IndexCriterion['macro_scattering011']= 0.70 
        #listOfKeptPercentPriori_IndexCriterion['macro_nu*fission0']= 0.80 
        #listOfKeptPercentPriori_IndexCriterion['macro_nu*fission1']= 0.80
        #listOfKeptPercentPriori_IndexCriterion['macro_fission0']= 0.80
        #listOfKeptPercentPriori_IndexCriterion['macro_fission1']= 0.80
        ### =========================================================================================
        ### =======================UOX-Gd cas for direct sparsity ===================================
        listOfKeptPercentPriori_IndexCriterion = {}
        listOfKeptPercentPriori_IndexCriterion['macro_totale0'] = 0.6 #0.96
        listOfKeptPercentPriori_IndexCriterion['macro_totale1'] = 0.8 #0.85
        listOfKeptPercentPriori_IndexCriterion['macro_absorption0']= 0.6 #0.96
        listOfKeptPercentPriori_IndexCriterion['macro_absorption1']= 0.6 #0.85
        listOfKeptPercentPriori_IndexCriterion['macro_scattering000'] = 0.6 #0.96       
        listOfKeptPercentPriori_IndexCriterion['macro_scattering001']= 0.8 #0.94
        listOfKeptPercentPriori_IndexCriterion['macro_scattering010']= 0.6 #0.85
        listOfKeptPercentPriori_IndexCriterion['macro_scattering011']= 0.7 #0.85
        listOfKeptPercentPriori_IndexCriterion['macro_nu*fission0']= 0.6 #0.96
        listOfKeptPercentPriori_IndexCriterion['macro_nu*fission1']= 0.6 #0.85
        listOfKeptPercentPriori_IndexCriterion['macro_fission0']= 0.6 #0.96
        listOfKeptPercentPriori_IndexCriterion['macro_fission1']= 0.6 #0.85
        
        
        
        keptPercent= listOfKeptPercentPriori_IndexCriterion[self.crossSectionName]

        self.keptCoefIndexes_Priori = self.getKeptCoefIndexes_PrioriCriterion_Index(self.listOfFinalCoefIndexes_arr, keptPercent)
        

        self.listOfFinalCoefIndexes_arr_priori=  [self.listOfFinalCoefIndexes_arr[i] for i in self.keptCoefIndexes_Priori]     
       
        
        matrixA = self.matrixA_full[:, self.keptCoefIndexes_Priori]### take only columns with indexes in keptIndexesAfterSparsity <==> put 0 to small coefficient
                                                            ### Used for least square to solve matrixA*a = b (over detemined system car dim(b) > dim(a),dim( matrixA) = dim(b)*dim(a)                                                         
                                                            
        ### ==============Extract from exact values====================
        self.FinalTuckerCoeffs_ViaPriori = self.FinalTuckerCoeffs[self.keptCoefIndexes_Priori] 
        ### ===========================================================  
        
        
        ###=================================================
        ### CHOICE 1 : Full least square (A'x = b, b : initial )
        ###=================================================
        self.A_sparse_priori = matrixA 
        self.FinalTuckerCoeffs_sparse_priori_fullLS = np.linalg.lstsq(self.A_sparse_priori,self.vectorOfExactValues)[0]
        
        ###=================================================
        ### CHOICE 2 : Reduced directly by a square system (A'x = b', size(x) =size(b')  )
        ### Not using least square
        ###=================================================
        b_sparse_priori_direct = self.vectorOfExactValues[self.keptCoefIndexes_Priori]       
        
        self.A_sparse_priori_direct = matrixA[self.keptCoefIndexes_Priori, :]  ### take only rows with indexes in self.keptCoefIndexes_Priori
        self.FinalTuckerCoeffs_sparse_priori_direct = np.linalg.solve(self.A_sparse_priori_direct, b_sparse_priori_direct)
        ###==================================
        
        
        ###=================================================
        ### CHOICE 3 : Reduced  by intermediate way via a tensor structure applied to reference points used for b 
        ### Using least square
        ###=================================================

        self.keptRefPointIndexes_tensorStructure = self.getKeptRefPointIndexes_PrioriSparsity(keptPercent) ### These are not indexes of coeffs
        self.A_sparse_priori_tensorStructure = matrixA[self.keptRefPointIndexes_tensorStructure, :]   ### take only rows with indexes in keptRefPointIndexes_tensorStructure 
        self.b_sparse_priori_tensorStructure = self.vectorOfExactValues[self.keptRefPointIndexes_tensorStructure]                                                                                                        ###==> eliminate the point not used in the new system with tensor structure
        print "len(self.b_sparse_priori_tensorStructure)", len(self.b_sparse_priori_tensorStructure)
        #raw_input()
        
        self.FinalTuckerCoeffs_sparse_priori_tensorStructure = np.linalg.lstsq(self.A_sparse_priori_tensorStructure, self.b_sparse_priori_tensorStructure)[0]
        

        print "===========A PRIORI================"
        reducedNb = len(self.listOfFinalCoefIndexes_arr) - len(self.listOfFinalCoefIndexes_arr_priori)
        print "Percentage of a PRIORI reduction = %s %%" %(100.0*reducedNb/len(self.listOfFinalCoefIndexes_arr))
        print "=========END A PRIORI=================="
        
##===================================================   
    def evaluate(self, point) :
        """
      Fonction qui évalue la décompostion de Tucker pour la section considérée en un point donné
      Argument :
        "point" [un list de "d" (=dimension) flottants] : point sur le quel on évalue
      Retour :
        "approxValue" [un flottant]: la veleur d'évaluation de section sur le "point" via Tucker

        """
        approxValue = 0.0
   
        for I in range(len(self.listOfFinalCoefIndexes_arr)) :
            d = 1.0
           
            for Axis_k in range(self.dimension) :                
                fct = self.listOfBasisFcts[str(Axis_k)][self.listOfFinalCoefIndexes_arr[I][Axis_k]]
                d = d*LI.getInterpolation(fct, point[Axis_k])
               
            approxValue = approxValue +  self.FinalTuckerCoeffs[I]*d       
     
        return approxValue
 
##===================================================   
    def evaluate_sparse_priori(self, point) :
        """
      Fonction qui évalue la décompostion de Tucker en un point après avoir creusé a priori les coefficients 
      Argument :
        point [une liste de "d" (=dimension) flottants] : point sur le quel on évalue
      Retour :
        + approxValue_direct [un flottant]: la veleur d'évaluation où les coefficients dans A'a' = b' 
          sont directement déterminés (sans moindres carrées)
        +  approxValue_ls [un flottant]: la veleur d'évaluation où les coefficients dans A'a' = b' 
          sont déterminés avec la méthode de moindres carrées et b' = b
        +  approxValue_ls_tensorStructure [un flottant]: la veleur d'évaluation où les coefficients dans A'a' = b' 
          sont déterminés avec la méthode de moindres carrées et b' <> b, size(b')<size(b), les indices gardés dans
          b sont tensoriels.        

        """
        approxValue_direct = 0.0
        approxValue_ls = 0.0
        approxValue_ls_tensorStructure = 0.0

        for I in range(len(self.listOfFinalCoefIndexes_arr_priori)) :
            d = 1.0
           
            for Axis_k in range(self.dimension) :
                fct = self.listOfBasisFcts[str(Axis_k)][self.listOfFinalCoefIndexes_arr_priori[I][Axis_k]]
                d = d*LI.getInterpolation(fct, point[Axis_k])
               
            approxValue_direct = approxValue_direct +  self.FinalTuckerCoeffs_sparse_priori_direct[I]*d   
            approxValue_ls = approxValue_ls +  self.FinalTuckerCoeffs_sparse_priori_fullLS[I]*d 
            approxValue_ls_tensorStructure = approxValue_ls_tensorStructure + self.FinalTuckerCoeffs_sparse_priori_tensorStructure[I]*d
        #return approxValue_ls_tensorStructure    
        return approxValue_direct, approxValue_ls, approxValue_ls_tensorStructure
##===================================================   
    def evaluate_sparse_posteriori(self, point) :   
        """
      Fonction qui évalue la décompostion de Tucker en un point après avoir creusé a posteriori les coefficients 
      Argument :
        point [un list de "d" (=dimension) flottants] : point sur le quel on évalue
      Retour :
         "approxValue" [un flottant]: la veleur d'évaluation de section sur le "point" via Tucker creusé a posteriori  

        """
        approxValue = 0.0
        for I in range(len(self.listOfFinalCoefIndexes_arr_posteriori)) :
            d = 1.0
           
            for Axis_k in range(self.dimension):
                fct = self.listOfBasisFcts[str(Axis_k)][self.listOfFinalCoefIndexes_arr_posteriori[I][Axis_k]]
                d = d*LI.getInterpolation(fct, point[Axis_k])
               
            approxValue = approxValue +  self.FinalTuckerCoeffs_sparse_posteriori[I]*d    
     
        return approxValue



##=================================================== 
    def getIntersection_DirectEliminatedPoints_TensorStructurePoints(self):
        ### Refind the eliminated points/indexes by a priori sparsity
        allIndexes = np.arange(self.CardinalOfBasisFcts)  
        mask = np.ones(len(allIndexes), dtype=bool) # all elements included/True.
        mask[self.keptCoefIndexes_Priori] = False    
        
        i_eliminated_Priori = allIndexes[mask] 
        
        print "len(i_eliminated_Priori)", len(i_eliminated_Priori)
        print "len(self.keptRefPointIndexes_tensorStructure)", len(self.keptRefPointIndexes_tensorStructure)
        
        #raw_input()
        
        commonIndexes = list(set(i_eliminated_Priori).intersection(self.keptRefPointIndexes_tensorStructure))
        
        
        nb_difference_tensorStruct_DirectPrio = len(self.keptRefPointIndexes_tensorStructure) - len(self.keptCoefIndexes_Priori)
        
        print "len(commonIndexes)", len(commonIndexes)
        print "nb_difference_tensorStruct_DirectPrio = ", nb_difference_tensorStruct_DirectPrio
        #raw_input()
        
        if len(commonIndexes) >= nb_difference_tensorStruct_DirectPrio:
          ### Choose the first elements in common points
            #new_eliminated_indexes = commonIndexes[:nb_difference_tensorStruct_DirectPrio] ### Take 702 points from the common points
          ### Choose randomly 
            new_eliminated_indexes =random.sample(set(commonIndexes), nb_difference_tensorStruct_DirectPrio)
          
            self.new_kept_indexes_tensor_direct = [i for i in self.keptRefPointIndexes_tensorStructure if i not in new_eliminated_indexes]
            
            print "len(self.new_kept_indexes_tensor_direct)", len(self.new_kept_indexes_tensor_direct)
            print self.new_kept_indexes_tensor_direct
            #raw_input()
            
            self.A_sparse_priori_tensorStructure_new = self.A_sparse_priori[self.new_kept_indexes_tensor_direct, :]   ### take only rows with indexes in keptRefPointIndexes_tensorStructure 
            self.b_sparse_priori_tensorStructure_new = self.vectorOfExactValues[self.new_kept_indexes_tensor_direct]                                                                 
  
            
            self.FinalTuckerCoeffs_sparse_priori_tensorStructure_new = np.linalg.solve(self.A_sparse_priori_tensorStructure_new, self.b_sparse_priori_tensorStructure_new)
            
            
            print "len(self.FinalTuckerCoeffs_sparse_priori_tensorStructure_new) =", len(self.FinalTuckerCoeffs_sparse_priori_tensorStructure_new)
            print "len(new_eliminated_indexes)", len(new_eliminated_indexes)
            print "self.A_sparse_priori_tensorStructure_new.shape =", self.A_sparse_priori_tensorStructure_new.shape
            print "len(self.b_sparse_priori_tensorStructure_new) = ", len(self.b_sparse_priori_tensorStructure_new)
            print "len(self.FinalTuckerCoeffs_sparse_priori_tensorStructure_new) =", len(self.FinalTuckerCoeffs_sparse_priori_tensorStructure_new)
            #raw_input()
            
            
        else:
            print "Not possible to take all points in direct elimination for reducing the size of tensor structure"
###===================================================
    def evaluate_sparse_priori_tensor_direct(self, point) :
        approxValue_tensor_direct = 0.0
        

        for I in range(len(self.FinalTuckerCoeffs_sparse_priori_tensorStructure_new)) :
            d = 1.0
           
            for Axis_k in range(self.dimension) :
                fct = self.listOfBasisFcts[str(Axis_k)][self.listOfFinalCoefIndexes_arr_priori[I][Axis_k]]
                d = d*LI.getInterpolation(fct, point[Axis_k])
               
            approxValue_tensor_direct = approxValue_tensor_direct +  self.FinalTuckerCoeffs_sparse_priori_tensorStructure_new[I]*d   
   
        return approxValue_tensor_direct 
        
        
        
####===================================================   
#### Partie de vérification sur le creusage a priori avec le structire tensoriel 
###===================================================   
    #def verifyLeastSquare(self, A, coefs, b):
        #deltaError = np.dot(A,coefs)-b
        #deltaError_L2 = np.sqrt(sum (np.power(deltaError, 2)))
        #return deltaError_L2

    #def normInf(self, coefs1, coefs2):
        #return max(abs(coefs1 - coefs2))
    
    #def normL2(self, coefs1, coefs2):
        #return np.sqrt(sum(np.power(coefs1 - coefs2,2))) /len(coefs1)
        
        
    #def getVerification(self):
        #print "self.FinalTuckerCoeffs.size = ", len(self.FinalTuckerCoeffs)
        #print "self.FinalTuckerCoeffs_sparse_priori.size = ", len(self.FinalTuckerCoeffs_sparse_priori)

      
        #diffCoefInf_LS = self.normInf(self.FinalTuckerCoeffs_ViaPriori, self.FinalTuckerCoeffs_sparse_priori)
        #diffCoefInf_LS_reduit = self.normInf(self.FinalTuckerCoeffs_ViaPriori, self.FinalTuckerCoeffs_sparse_priori_tensorStructure)
        #diffCoefInf_direct = self.normInf(self.FinalTuckerCoeffs_ViaPriori, self.FinalTuckerCoeffs_sparse_priori_direct)
        
        #diffCoefL2_LS = self.normL2(self.FinalTuckerCoeffs_ViaPriori, self.FinalTuckerCoeffs_sparse_priori)
        #diffCoefL2_LS_reduit = self.normL2(self.FinalTuckerCoeffs_ViaPriori, self.FinalTuckerCoeffs_sparse_priori_tensorStructure)
        #diffCoefL2_direct = self.normL2(self.FinalTuckerCoeffs_ViaPriori, self.FinalTuckerCoeffs_sparse_priori_direct)
        
        #print "========================"
        #print "diffCoefInf_LS =", diffCoefInf_LS
        #print "diffCoefInf_LS_reduit =", diffCoefInf_LS_reduit
        #print "diffCoefInf_direct = ", diffCoefInf_direct
        #print "diffCoefL2_LS = ", diffCoefL2_LS
        #print "diffCoefL2_LS_reduit = ", diffCoefL2_LS_reduit
        #print "diffCoefL2_direct = ", diffCoefL2_direct
        #print "========================"
        #raw_input()
        #fit1_1620 = self.verifyLeastSquare(self.A_sparse_priori,self.FinalTuckerCoeffs_sparse_priori, self.vectorOfExactValues)
        #fit2_1350 = self.verifyLeastSquare(self.A_sparse_priori_tensorStructure, self.FinalTuckerCoeffs_sparse_priori_tensorStructure, self.b_sparse_priori_tensorStructure)
        
        #fit1_withSol1350 = self.verifyLeastSquare(self.A_sparse_priori,self.FinalTuckerCoeffs_sparse_priori_tensorStructure, self.vectorOfExactValues)
        #fit2_withSol1620 = self.verifyLeastSquare(self.A_sparse_priori_tensorStructure, self.FinalTuckerCoeffs_sparse_priori, self.b_sparse_priori_tensorStructure)
        
        
        ##fit2_1350 = self.verifyLeastSquare(self.A_sparse_priori_tensorStructure, self.FinalTuckerCoeffs_sparse_priori_tensorStructure, self.b_sparse_priori_tensorStructure)
        
        
        
        
        #print "========================"
        #print "We have to have: fit1_1620 < fit1_withSol1350"
        #print "fit1_1620 = ", fit1_1620
        #print "fit1_withSol1350 =", fit1_withSol1350
        #print "========================"
        #print "We have to have: fit2_1350 < fit2_withSol1620"
        #print "fit2_1350 = ", fit2_1350
        #print "fit2_withSol1620 = ", fit2_withSol1620
        #print "========================"
        #print "We have to have: fit2_withSol1620 < fit1_1620"
        #print "fit2_withSol1620 = ", fit2_withSol1620
        #print "fit1_1620 = ", fit1_1620
        #print "========================"
        #raw_input()

###=================================================== 
    #def verifyByDefinedFct(self):
        #coefs = []
        #for I in range(len(self.listOfFinalIndexes_priori)):
            ##print "self.listOfFinalIndexes_priori[I] =", self.listOfFinalIndexes_priori[I]
            ##print "sum(self.listOfFinalIndexes_priori[I]) = ", sum(self.listOfFinalIndexes_priori[I])
            ##print "1.0/(sum(self.listOfFinalIndexes_priori[I])+1.0) =", 1.0/(sum(self.listOfFinalIndexes_priori[I])+1.0)
            ##raw_input()
            #coef = 1.0 #/(sum(self.listOfFinalIndexes_priori[I])+1.0)
            #coefs.append(coef)
        #coefs_exact_np = np.asarray(coefs)    
        #b = np.dot(self.A_sparse_priori_tensorStructure, coefs_exact_np)
        #sol_LS =  np.linalg.lstsq(self.A_sparse_priori_tensorStructure,b)[0] #self.FinalTuckerCoeffs_sparse_priori = np.linalg.lstsq(self.A_sparse_priori_tensorStructure,b)[0]
        
        #diffCoefInf = self.normInf(coefs_exact_np, sol_LS)
        

        #diffCoefL2 = self.normL2(coefs_exact_np, sol_LS)
        
        #delta_sol_exact = self.verifyLeastSquare(self.A_sparse_priori_tensorStructure, coefs_exact_np, b)
        #delta_sol_LS = self.verifyLeastSquare(self.A_sparse_priori_tensorStructure, sol_LS, b)
        
        
        #print "delta_sol_LS = %s must be smaller than delta_sol_exact = %s" %(delta_sol_LS, delta_sol_exact)
        
        #print "=======Case defined function to verify least square================"
        #print "diffCoefInf = ", diffCoefInf
        #print "diffCoefL2 = ", diffCoefL2
        #print "=================================="
        
        #raw_input()

            
####===================================================   
#### Fin de la partie de vérification sur le creusage a priori avec le structire tensoriel 
###===================================================             
            
