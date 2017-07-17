# -*- coding: utf-8 -*-
import numpy as np
import warnings
from scipy.special.orthogonal import p_roots, t_roots
import numpy as np
from scipy.fftpack import ifft
import os
import random
from random import uniform

### INPUTS:
### Axis_k --> un entier qui represente le paramètre qui est finement discretizé.
### listOfBorders ---> c'est un dictionnaire 
    # clé: c'est la cahine correspondant à un numéro d'axe  ex '0', '1' etc
    # valeur est une liste soit une liste simple de bornes [min, max], soit une liste de listes (decoupa l'axe)
    # ex: If we subdivide a parameter, its interval will be decomposed, for example [a,b] = [[a,c],[c,d], [d,b]]
### listOfNbPointsOnEachIntervalForFinalDiscretization_Axis_k ---> c'est une liste de listes d'entiers
      #
     # Define the total number of point used in the discretization for each interval/subinterval.
### listOfTuckerGridNodes_Axis_k ---> Initial interval : [lower border, upper border] ---> discretize ----> finally: discretized points (by "CC", by "T") included in listOfTuckerGridNodes_Axis_k
###=====================================
class domainD:
    ###=====================================
    ### This class is used to define 1 Tucker grid corresponding to the axis "k" (Axis_k)
    ### This grid is finely discretized in Axis_k and coarsely in other axes
    ### There are in total of "d" Tucker grids but this class is describred only a discretization of one grid
    ##======================================
    def __init__(self,  Axis_k, listOfBorders, listOfNbPointsOnEachIntervalForFinalDiscretization_Axis_k, listOfTuckerGridNodes_Axis_k) :
        """
          
          Argumentes:
             "Axis_k"  [entier] : represente le paramètre qui est finement discretizé.
             "listOfBorders" [dictionnaire] : 
                + clé [chaine de caractères] : désigne l'axe. ex :'0', '1'
                + valeur [liste de flottants] : soit une liste simple de bornes [min, max], soit une liste de listes (decoupa l'axe)]
             "listOfNbPointsOnEachIntervalForFinalDiscretization_Axis_k" [liste de listes d'entiers]: ce sont les nombres de points dans la discrétisation
             "listOfTuckerGridNodes_Axis_k" [un dictionnaire]: 
             + clé : c'est la chain de caractère correspondant à un numéro d'axe  ex '0', '1' etc
             + valeur : les coordonnées de points dans la discrétisation (pour une grille de Tucker) pour tous les axes
          Retourn:
            rien (constructeur)
        """ 
        self.listOfBorders = listOfBorders        
        self.dimension = len(listOfBorders)     
        self.listOfLowerBorders = []
       
        for i in range(self.dimension) :
            self.listOfLowerBorders.append(listOfBorders[str(i)][0])
            
        self.listOfUpperBorders = []
        for i in range(self.dimension) :
            self.listOfUpperBorders.append(listOfBorders[str(i)][-1])
            
        self.lx = []
        for i in range(self.dimension) :
            self.lx.append(self.listOfUpperBorders[i] - self.listOfLowerBorders[i])
            
        self.weights = []
        wi = []
        for i in range(0, self.dimension):
            self.weights.append(wi)

        
        self.listOfTuckerGridNodes_Axis_k = listOfTuckerGridNodes_Axis_k
        
        self.listOfNbPointsOnEachIntervalForFinalDiscretization_Axis_k = listOfNbPointsOnEachIntervalForFinalDiscretization_Axis_k    
        self.listOfNbPointsOnEachIntervalForFinalDiscretization_Axis_kTotal = listOfNbPointsOnEachIntervalForFinalDiscretization_Axis_k[:]
        #print "self.listOfNbPointsOnEachIntervalForFinalDiscretization_Axis_kTotal", self.listOfNbPointsOnEachIntervalForFinalDiscretization_Axis_kTotal
        #print "self.listOfNbPointsOnEachIntervalForFinalDiscretization_Axis_k[Axis_k]", self.listOfNbPointsOnEachIntervalForFinalDiscretization_Axis_k[Axis_k]
        #raw_input()
        ### listOfNbPoints = [n0, n1, ...[n_k1, ...,n_kp], n_(k+1), ..., n_d]
        nbrPTotal_Axis_k = int(self.sumArr(self.listOfNbPointsOnEachIntervalForFinalDiscretization_Axis_k[Axis_k]))
        self.listOfNbPointsOnEachIntervalForFinalDiscretization_Axis_kTotal[Axis_k] = nbrPTotal_Axis_k
        
        self.Axis_k = Axis_k   
        self.NbPointsTotal = self.produitArr(self.listOfNbPointsOnEachIntervalForFinalDiscretization_Axis_kTotal)
##=====================================
    def produitArr(self, u):
        """
        focntion qui fait le produit des elements d'un tableau
        
        Argument:
            "u" : [liste ou array de flottants]
            
        retourne:
            un flottant qui est le produit des elements de le tableau "u"
        """
        n = len(u)
        s = 1.0
        for i in range(0,n):
            s = s*u[i]
        return s
##===================================== 
    def sumArr(self, u):
        """
        focntion qui fait la somme des elements d'un tableau
        
        Argument:
            "u" : [liste ou array de flottants]
            
        retourne:
            un flottant qui est la somme des elements de le tableau "u"
        """
        n = len(u)
        s = 0
        for i in range(n):
            s = s + u[i]
        return s

##=====================================    
    def getCoefAffine(self, xLOld, xUOld, xLNew, xUNew) :
        """
        focntion qui fait l'inverse de la transformation affine de l'intervalle [xLOld, xUOld] à la nouvelle intervalle [xLNew, xUNew]
        + xOld dans [a, b] ---> xNew dans [c, d] => xNew = (d-c)/(b-a)*xOld + (bc-ad)/(b-a)
        + coefficients de transformation :   (d-c)/(b-a) et   (bc-ad)/(b-a)
        + C-a-d:
     
            "xOld" est dans  [a, b] ---> "xNew" est dans [c, d] => xNew = (d-c)/(b-a)*xOld + (bc-ad)/(b-a)
              => xOld = (b-a)/d-c) *xNew + (ad-bc)/(d-c)
            f(xOld) sur [a, b] <==> f((b-a)/d-c) *xNew + (ad-bc)/(d-c)) sur [c,d]

        
        Argument:
            "xLOld" [un flottant] : min de l'intervalle initiale
            "xUOld" [un flottant] : max de l'intervalle initiale
            
            "xLNew" [un flottant] : min de la nouvelle intervalle 
            "xUNew" [un flottant] : max de la nouvelle intervalle 
            
        retourne:
            une liste de 2 flottants correspondant aux coefficients de transformation affine de  [xLNew, xUNew] à  [xLOld, xUOld]
        """ 
        return [(xUOld - xLOld)/(xUNew - xLNew), (xLOld*xUNew -xUOld* xLNew)/(xUNew - xLNew)]
##=====================================     
    def transformeAffine(self, x, coef) :
        """
        focntion qui calcule le coordonnée correspondant d'un point x par une transformation affine y = ax + b 
        
        Argument:
            "x" [un flottant] : coordonnées d'un point
            coef [une liste de deux flottants] : "a" et "b" dans la transformation affine y = ax + b 
            
        retourne:
            un flottant : le nouveau coordonnée de x après avoir transformé
        """
      
        return coef[0]*x + coef[1]
##===================================== 
    def getDxiAfterTransf(self, xi, Axis_i):
        """
        focntion qui calcule le nouveau coordonnée d'un point x par une transformation affine de l'intervalle [-1,1] à une nouvelle intervalle
        (+ [-1,1] intervalle par défaut de numpy dans discrétisation par des points de CC,
         + la nouvelle intervalle est les bornes dans notre discrétisation) 
        
        Argument:
            "xi" [un flottant ou un array] : coordonnées d'un poit ou plusieur points  
            Axis_i [un entier] : l'axe sur laquelle on travaille
            
        retourne:
            rien (le nouveau coordonnées sous la forme [un flottant ou un array] est mis dans la discrétisation de l'axe "Axis_i"
        """  
      
        coefX = self.getCoefAffine(self.listOfLowerBorders[Axis_i], self.listOfUpperBorders[Axis_i], -1.0, 1.0)
        self.listOfTuckerGridNodes_Axis_k[Axis_i] = self.transformeAffine(xi, coefX)
    
            
    def discretizeTrapezoidal(self, Axis_i): #no on [-1,1]
        if self.listOfNbPointsOnEachIntervalForFinalDiscretization_Axis_kTotal[Axis_i] !=  2 :
            print self.listOfNbPointsOnEachIntervalForFinalDiscretization_Axis_kTotal[Axis_i]
            raw_input()
            warnings.warn("error",'nombre of points must be equal to 2')
        else :
            self.listOfTuckerGridNodes_Axis_k[Axis_i] = np.array([self.listOfLowerBorders[Axis_i], self.listOfUpperBorders[Axis_i]])
            self.weights[Axis_i] = [self.lx[Axis_i]/2, self.lx[Axis_i]/2]
           
##===================================== 
    ### Trapezoidal method for two points c & d where a <=c < d <=b 
##=====================================    
    def discretizeTrapezoidal_bis(self, Axis_i): #
        Points = self.listOfTuckerGridNodes_Axis_k[Axis_i]
  
        if len(Points) != 2 or ( Points[0] < self.listOfLowerBorders[Axis_i]) or ( Points[1] > self.listOfUpperBorders[Axis_i]): #self.listOfNbPointsOnEachIntervalForFinalDiscretization_Axis_kTotal[Axis_i] !=  2 :
            warnings.warn("error",'nombre of points must be equal to 2 or out of bound')
        else :
    
            xL_i = self.listOfLowerBorders[Axis_i]
            xU_i = self.listOfUpperBorders[Axis_i]
      
            w0 = -(self.lx[Axis_i]/2)*(xL_i  +xU_i - 2* Points[1])/( Points[1] -  Points[0])
            w1 = (self.lx[Axis_i]/2)*(xL_i  +xU_i - 2* Points[0])/( Points[1] -  Points[0])
   
            self.weights[Axis_i] = [w0, w1] #[1.0, 1.0]
            
 ##==========Ajoute le 27/01/2015=========================== 
 ### f(x)/ x0 = a < x1 < ...< xi < ... < x(n-1 )= b
 ### f(x) = [f(x0)+f(x1)]*(x1-x0)/2 + [f(x1)+f(x2)]*(x2-x1)/2 + ... + [f(x_k)+f(x_(k-1))]*(x_k-x_(k-1))/2 + ...[f(x_(n-2))+f(x_(n))]*(xn-x_(n-1))/2
 ### = f(x0)*(x1-x0)/2 + f(x1) *(x2-x0)/2 + ...+ f(xk)*(x_(k+1)- x_(k-1))/2 + f(x_(n-1))*(x_(n-1)-x_(n-2))/2
 ##=========================================================
    def discretizeTrapezoidal_ter(self, Axis_i): #
       
        Points = self.listOfTuckerGridNodes_Axis_k[Axis_i]
      
        if  ( Points[0] != self.listOfLowerBorders[Axis_i]) or ( Points[-1] != self.listOfUpperBorders[Axis_i]): #self.listOfNbPointsOnEachIntervalForFinalDiscretization_Axis_kTotal[Axis_i] !=  2 :
            warnings.warn("error",'first point and final point must be extrema of direction')
        else :
            n = len(Points)
            w = np.ones(n)
            for i in range(n):
                if (i==0) :
                    w[i] = (Points [i+1] - Points[i])/2
                elif (i== n-1):
                    w[i] = (Points [i] - Points[i-1])/2
                else :
                    w[i] = (Points [i+1] - Points[i-1])/2            
   
            self.weights[Axis_i] = w
                 
      
##=====================================            
    def discretizeSimpson(self, Axis_i):
        #if (Axis_i >= self.dimension) or (Axis_i <0):
            #warnings.warn("error",'index of axis is out of bound')
        if (self.listOfNbPointsOnEachIntervalForFinalDiscretization_Axis_kTotal[Axis_i] - 1 )%2 !=  0 :
            warnings.warn("error",'nombre of points must be odd')
        else :
            x = np.zeros(self.listOfNbPointsOnEachIntervalForFinalDiscretization_Axis_kTotal[Axis_i])
            w = np.ones((self.listOfNbPointsOnEachIntervalForFinalDiscretization_Axis_kTotal[Axis_i]))
        
            hx = self.lx[Axis_i]/(self.listOfNbPointsOnEachIntervalForFinalDiscretization_Axis_kTotal[Axis_i] - 1)
            for i in range(self.listOfNbPointsOnEachIntervalForFinalDiscretization_Axis_kTotal[Axis_i]):
                x[i] = self.listOfLowerBorders[Axis_i] + i*hx
                self.listOfTuckerGridNodes_Axis_k[Axis_i] = x
        
            for i in range(1,self.listOfNbPointsOnEachIntervalForFinalDiscretization_Axis_kTotal[Axis_i] - 1):
                if i%2 == 0 :
                    w[i] = 2
                else :
                    w[i] = 4
            self.weights[Axis_i] = w*hx/3
            

##=====================================    
    
    def discretizeLegendre(self, Axis_i):
        [xi, self.weights[Axis_i]] = p_roots(self.listOfNbPointsOnEachIntervalForFinalDiscretization_Axis_kTotal[Axis_i])

        self.getDxiAfterTransf(xi, Axis_i) # Update self.listOfTuckerGridNodes_Axis_k[Axis_i]

        self.weights[Axis_i] = (self.lx[Axis_i]/2)*self.weights[Axis_i]        
##=====================================
        
    def discretizeChebyshev(self, Axis_i):
        [xi, self.weights[Axis_i]] = t_roots(self.listOfNbPointsOnEachIntervalForFinalDiscretization_Axis_kTotal[Axis_i])
        self.getDxiAfterTransf(xi, Axis_i) # Update self.listOfTuckerGridNodes_Axis_k[Axis_i]
        self.listOfcoefW[Axis_i] = self.lx[Axis_i]/2
    #self.listOfNbPointsOnEachIntervalForFinalDiscretization_Axis_kTotal[Axis_i]= self.listOfNbPointsOnEachIntervalForFinalDiscretization_Axis_kTotal[Axis_i]
 
    def clencurt(self,n1):

        """
        source : http://www.scientificpython.net/1/post/2012/04/clenshaw-curtis-quadrature.html 
        focntion qui calcule les coordonées des points de Clenshaw Curtis et ses poids dans l'intervalle [-1,1]

        Argument:
            "n1" [un entier] : le nombre de points de  Clenshaw Curtis dans [-1,1]
            
        retourne:
            "x" [un array] : les coordonnées des points CC
            "w" [un array] : les poids associé
        """  
        if n1 == 1:
            x = 0
            w = 2
        else:
            n = n1 - 1
            C = np.zeros((n1,2))
            k = 2*(1+np.arange(np.floor(n/2)))
            C[::2,0] = 2/np.hstack((1, 1-k*k))
            C[1,1] = -n
            V = np.vstack((C,np.flipud(C[1:n,:])))
            F = np.real(ifft(V, n=None, axis=0))
            x = F[0:n1,1]
            w = np.hstack((F[0,0],2*F[1:n,0],F[n,0]))

        return x,w
##=====================================     

    def discretizeClenshawCurtis(self, Axis_i):
        """
        
        focntion qui calcule les coordonées des points de Clenshaw Curtis et ses poids associés pour un axe SANS coupage
        (si sur cet axe, on utilise les points CC)
        Argument:
            Axis_i [un entier] : l'axe sur laquelle on travaille
            
            
        retourne:
            rien (les points de CC et ses poids sont directement mis à jour dans la discrétisation de l'axe " Axis_i")
        """
      
        [xi, self.weights[Axis_i]] =  self.clencurt(self.listOfNbPointsOnEachIntervalForFinalDiscretization_Axis_kTotal[Axis_i])
        self.getDxiAfterTransf(xi, Axis_i) # Mettre à jours self.listOfTuckerGridNodes_Axis_k[Axis_i]
        self.weights[Axis_i] = (self.lx[Axis_i]/2)*self.weights[Axis_i]
       
##=====================================    
    def discretizeClenshawCurtisBis(self, Axis_i) :#listOfExtremesForCouplage, listOfNbrPsForCouplage):
        """
        
        focntion qui calcule les coordonées des points de Clenshaw Curtis et ses poids associés pour un axe AVEC le coupage
        (si sur cet axe, on utilise les points CC)
        Argument:
            Axis_i [un entier] : l'axe sur laquelle on travaille
            
            
        retourne:
            rien (les points de CC et ses poids sont directement mis à jour dans la discrétisation de l'axe " Axis_i")
        """ 
        x_CC = []
        w_CC = []
        
        #print "========================="
        #print "self.listOfBorders[str(Axis_i)]", self.listOfBorders[str(Axis_i)]
        #print "self.listOfNbPointsOnEachIntervalForFinalDiscretization_Axis_k", self.listOfNbPointsOnEachIntervalForFinalDiscretization_Axis_k
        #print "========================="
        #raw_input()
  
        if (len(self.listOfBorders[str(Axis_i)]) -1) !=  len(self.listOfNbPointsOnEachIntervalForFinalDiscretization_Axis_k[Axis_i]) :
            warnings.warn("error",'nb of extremes must be equal to nb of elements in list of nb of points/sub-interval ')
        else :
            for i in range(len(self.listOfNbPointsOnEachIntervalForFinalDiscretization_Axis_k[Axis_i])):
                [xi, wi] = self.clencurt(self.listOfNbPointsOnEachIntervalForFinalDiscretization_Axis_k[Axis_i][i])
       
                coefi = self.getCoefAffine(self.listOfBorders[str(Axis_i)][i], self.listOfBorders[str(Axis_i)][i+1], -1, 1) 
                xi = self.transformeAffine(xi, coefi) 
                wi = ((self.listOfBorders[str(Axis_i)][i+1] - self.listOfBorders[str(Axis_i)][i])/2.0)*wi
                ### Faire la concaténation entre les segments du coupage     
                x_CC = x_CC + list(xi)
                w_CC = w_CC + list(wi)
                #print "Axis_i =", Axis_i
                #print "i =", i
                #print "xi =", xi
                #print "wi = ", wi
                #print "x_CC = ", x_CC
                #print "w_CC =", w_CC
                #raw_input()
        self.listOfTuckerGridNodes_Axis_k[Axis_i] = np.asarray(x_CC)
        
        self.weights[Axis_i] = w_CC
        
        #print "======================="
        #print "self.listOfTuckerGridNodes_Axis_k[Axis_i] =", self.listOfTuckerGridNodes_Axis_k[Axis_i]
        #print "======================="
        #raw_input()
##=====================================        
    def discretizeLegendreBis(self, Axis_i) :
            x_L = []
            w_L = []
            if (len(self.listOfBorders[Axis_i]) -1) !=  len(self.listOfNbPointsOnEachIntervalForFinalDiscretization_Axis_k[Axis_i]) :
                warnings.warn("error",'nb of extremes must be equal to nb of elements in list of nb of points/sub-interval ')
            else :
                for i in range(len(self.listOfNbPointsOnEachIntervalForFinalDiscretization_Axis_k[Axis_i])):
                    [xi, wi] = p_roots(self.listOfNbPointsOnEachIntervalForFinalDiscretization_Axis_k[Axis_i][i])
                    coefi = self.getCoefAffine(self.listOfBorders[Axis_i][i], self.listOfBorders[Axis_i][i+1], -1, 1) 
                    xi = self.transformeAffine(xi, coefi) 
                    wi = ((self.listOfBorders[Axis_i][i+1] - self.listOfBorders[Axis_i][i])/2.0)*wi
                         
                    x_L = x_L + list(xi)
                    w_L = w_L + list(wi)
            self.listOfTuckerGridNodes_Axis_k[Axis_i] = np.asarray(x_L)
            self.weights[Axis_i] = w_L
##=====================================       
    def discretizeDomainEntire(self, listOfMethodsIntegrals):
        """
        
        focntion qui réalise les discrétisations (points + ses poids) selon le choix de méthode par axe
        Argument:                   
            "listOfMethodsIntegrals" [list of caractère] :  pour décrire les méthodes d'intégration pour chaque axe d'UNE grille de Tucker
        retourne:
            rien (une grille de Tucker est crée avec la "d" discrétisations sur "d" axes).
            Dans ce travail, un seul axe est discrétisé par "CC" ou "CCbis" (si coupage), 
            les autres par "T" (bis, ter) )
        """ 
       
        for Axis_i in range(0, self.dimension) :
            method = listOfMethodsIntegrals[Axis_i]
            
            if method == "T" :
                self.discretizeTrapezoidal(Axis_i)
                
            elif method == "Tbis" :
                self.discretizeTrapezoidal_bis(Axis_i)
                
            elif method == "Tter" :
                self.discretizeTrapezoidal_ter(Axis_i)
                
            elif method == "CC" :
                self.discretizeClenshawCurtis(Axis_i)
                
            elif method == "CCbis" :
                self.discretizeClenshawCurtisBis(Axis_i)   
                
            elif method == "Ch" :
                self.discretizeChebyshev(Axis_i)
                
            elif method == "L" :
                self.discretizeLegendre(Axis_i)
                
            elif method == "Lbis":
                self.discretizeLegendreBis(Axis_i)
      
                                
            elif method == "SS" :
                self.discretizeSimpson(Axis_i)
                
            else :
                warnings.warn("error",'out of method existe')
                
        #print "Après discrétisation numérique :"
        #print "self.listOfTuckerGridNodes_Axis_k", self.listOfTuckerGridNodes_Axis_k
        #print "self.weights ", self.weights
        #raw_input()

##=====================================                
    def transfIndexIntToArray(self, I, d, listOfNbPoints): 
        """
        
        focntion qui transforme un indice entier à une liste I[entier] ---> [i_0, ...,i_(d-1)]
         Argument:
         I [un entier] : indice de points sur une grille en "d" dimensions
         d [un entier] : la dimension = le nombre de paramètres CRN
         listOfNbPoints [une liste de entiers] : la base pour transformer l'indice
        
            
        retourne:
        "i" : une liste de "d" entier qui désigne la transformation de I dans la base "listOfNbPoints"
            
        """ 
        I = int(I)
        #i = np.zeros(d)
        i = []
        ii = [0]
        for k in range(0,d):
            i.append(ii)
            
        index = range(0,d)
        indexInverse = index[::-1] 
        
        for j in indexInverse:
            if j > 0 :
                i[j] = int(I/self.produitArr(listOfNbPoints[:j] ))
                I = I - i[j]* self.produitArr(listOfNbPoints[:j] )
                I = int(I)
        
            else :
                i[j] = int(I) 
        
        return i
##=====================================        
       
    def getXi(self, i):
        """

         Argument:
         i [un entier] : indice de l'axe         
            
        retourne:
        x [une liste des flottants] : les coordonnées des points dans la discrétisation de l'axe "i"
            
        """ 
        x = []
        
        for j in range(0, self.dimension):
            x.append(self.listOfTuckerGridNodes_Axis_k[j][i[j]])           
        return x

##=====================================    
    def getMatrixOfCrossSectionMkI(self, crossSection):
        """
         Argument:
         crossSection :  un dictionnaire avec la clé est le nombre de section qui est lue de dkzip      
            
        retourne:
        MkI [une matrice] : les valeurs de section efficace "crossSection" sur la grille de Tucker définie dans cette classe
                          : Cette matrice est utilisée dans la décomposition de Karhunen-Loève discrète
                          (taille : nbr de points dans l'axe principal * 2^{d-1}  )  
        """ 
   
        self.nbOfGroupedPointsY = self.NbPointsTotal/self.listOfNbPointsOnEachIntervalForFinalDiscretization_Axis_kTotal[self.Axis_k]
        listOfNbPointsNew = self.listOfNbPointsOnEachIntervalForFinalDiscretization_Axis_kTotal[:] ##Faire une copie, n'utilise pas "=" car il y aura probleme
        listOfNbPointsNew.remove(listOfNbPointsNew[self.Axis_k])
     
        M0 = np.zeros( self.NbPointsTotal)
        MkI = M0.reshape(self.listOfNbPointsOnEachIntervalForFinalDiscretization_Axis_kTotal[self.Axis_k], self.nbOfGroupedPointsY)
           
        for ik in range(0, int(self.listOfNbPointsOnEachIntervalForFinalDiscretization_Axis_kTotal[self.Axis_k])):
            for iI in range(0, int(self.nbOfGroupedPointsY)):
                i = self.transfIndexIntToArray(iI, self.dimension - 1, listOfNbPointsNew)
                i.insert(self.Axis_k,ik)
                x = self.getXi(i)
                MkI[ik,iI] = crossSection.evaluate(x)
     
        return MkI
##=====================================        
##=====================================    
      
    def getProdWi(self,i):
        """
         Argument:
         i [array of entiers ] :  les entiers sont des indices de coordonnées d'un point dans l'axe correspondant
        retourne:
        prodWi [un flottant] : le poid d'intégration pour le point indecé par "i" dans la méthode d'intégration dans multidimensions
        """ 
        prodWi = 1.0        
        for j in range(0, self.dimension):
            prodWi = prodWi*self.weights[j][i[j]]
        return prodWi
##=====================================    
        
    def getWeights(self):
        """
        fonction qui calcule les poids d'intégration associés aux matrice des valeurs de la section efficace sur la grille de Tucker
         Argument:
        rien
        retourne : 
        "WkI"  [une matrice des flottants]  : matrice des poids d'intégration (taille : nbr de points dans l'axe principal * 2^{d-1}  )        
        """ 
        
      
        self.minWeight_Axis_k = min(self.weights[self.Axis_k])
        
        self.minWeight_groupedAxis = 1E12
        
        self.nbOfGroupedPointsY = int(self.NbPointsTotal/self.listOfNbPointsOnEachIntervalForFinalDiscretization_Axis_kTotal[self.Axis_k])
        listOfNbPointsNew = self.listOfNbPointsOnEachIntervalForFinalDiscretization_Axis_kTotal[:] ##Faire une copie, n'utilise pas "=" car il y aura probleme
        listOfNbPointsNew.remove(listOfNbPointsNew[self.Axis_k])
        
        W = np.zeros(self.NbPointsTotal) ##self.listOfNbPointsOnEachIntervalForFinalDiscretization_Axis_kTotal[self.Axis_k]*self.nbOfGroupedPointsY)
        WkI = W.reshape(self.listOfNbPointsOnEachIntervalForFinalDiscretization_Axis_kTotal[self.Axis_k], self.nbOfGroupedPointsY)
        for ik in range(0, self.listOfNbPointsOnEachIntervalForFinalDiscretization_Axis_kTotal[self.Axis_k]):
            for iI in range(0, self.nbOfGroupedPointsY):
                i = self.transfIndexIntToArray(iI, self.dimension - 1, listOfNbPointsNew)
                i.insert(self.Axis_k,ik)
                WkI[ik,iI] = self.getProdWi(i)
                Weight_iI = WkI[ik,iI]/self.weights[self.Axis_k][ik] ### 
                self.minWeight_groupedAxis = min(self.minWeight_groupedAxis, Weight_iI)
                
        #print "self.weights[self.Axis_k]", self.weights[self.Axis_k]
        #print " self.minWeight_Axis_k",  self.minWeight_Axis_k
        #print "self.minWeight_groupedAxis = ", self.minWeight_groupedAxis
        #raw_input()
        
        return WkI
###=====================================     
