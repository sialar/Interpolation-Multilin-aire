# -*- coding: utf-8 -*-
import numpy as np
import warnings
from scipy.special.orthogonal import p_roots, t_roots
import numpy as np
from scipy.fftpack import ifft
import os
import random
from random import uniform
###=====================================
### This class is used to define 1 Tucker grid corresponding to the axis "k" (Axis_k)
### This grid is finely discretized in Axis_k and coarsely in other axes
### There are in total of "d" Tucker grids.
##======================================
### INPUTS:
### Axis_k --> fine discretization
### listOfBorders ---> min, max values of each parameter
### If we subdivide a parameter, its interval will be decomposed, for example [a,b] = [[a,c],[c,d], [d,b]]
### listOfNbPointsOnEachIntervalForFinalDiscretization_Axis_k ---> Define the total number of point used in the discretization for each interval/subinterval.
###listOfTuckerGridNodes_Axis_k --->
###=====================================
class domainD:
    def __init__(self,  Axis_k, listOfBorders, listOfNbPointsOnEachIntervalForFinalDiscretization_Axis_k, listOfTuckerGridNodes_Axis_k) :
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
        n = len(u)
        s = 1.0
        for i in range(0,n):
            s = s*u[i]
        return s
##===================================== 
    def sumArr(self, u):
        n = len(u)
        s = 0
        for i in range(n):
            s = s + u[i]
        return s
##=====================================
   ## xOld \in [a, b] ---> xNew \in [c, d] => xNew = (d-c)/(b-a)*xOld + (bc-ad)/(b-a)
   ##               => xOld = (b-a)/d-c) *xNew + (ad-bc)/(d-c)
   ## f(xOld)/[a, b] ---> f((b-a)/d-c) *xNew + (ad-bc)/(d-c))/[c,d]
##=====================================    
    def getCoefAffine(self, xLOld, xUOld, xLNew, xUNew) :
        return [(xUOld - xLOld)/(xUNew - xLNew), (xLOld*xUNew -xUOld* xLNew)/(xUNew - xLNew)]
    
    def transformeAffine(self, x, coef) :
        return coef[0]*x + coef[1]
##===================================== 
# Source : http://en.wikipedia.org/wiki/Gaussian_quadrature
# https://github.com/scipy/scipy/blob/v0.13.0/scipy/special/orthogonal.py#L713 
##=====================================
    
    def getDxiAfterTransf(self, xi, Axis_i):
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
    ### Methode Trapezoidal dans le cas deux points ne sonts pas extrêmes. a <=c < d <=b
    ### S = Trapezoidal (a->b) = -(a+b-2d)/(d-c)*f(c) + (a+b-2c)/(d-c)*f(d)
    ### Existe meuilleur méthode ??? Points milieux ?
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
 ### Objectif : utilise des points physique de Burnup
 ### f(x)/ x0 = a < x1 < ...< xi < ... < x(n-1 )= b
 ### f(x) = [f(x0)+f(x1)]*(x1-x0)/2 + [f(x1)+f(x2)]*(x2-x1)/2 + ... + [f(x_k)+f(x_(k-1))]*(x_k-x_(k-1))/2 + ...[f(x_(n-2))+f(x_(n))]*(xn-x_(n-1))/2
 ### = f(x0)*(x1-x0)/2 + f(x1) *(x2-x0)/2 + ...+ f(xk)*(x_(k+1)- x_(k-1))/2 + f(x_(n-1))*(x_(n-1)-x_(n-2))/2
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

        self.getDxiAfterTransf(xi, Axis_i) # Mettre à jours self.listOfTuckerGridNodes_Axis_k[Axis_i]

        self.weights[Axis_i] = (self.lx[Axis_i]/2)*self.weights[Axis_i]        
##=====================================
        
    def discretizeChebyshev(self, Axis_i):
        [xi, self.weights[Axis_i]] = t_roots(self.listOfNbPointsOnEachIntervalForFinalDiscretization_Axis_kTotal[Axis_i])
        self.getDxiAfterTransf(xi, Axis_i) # Mettre à jours self.listOfTuckerGridNodes_Axis_k[Axis_i]
        self.listOfcoefW[Axis_i] = self.lx[Axis_i]/2
    #self.listOfNbPointsOnEachIntervalForFinalDiscretization_Axis_kTotal[Axis_i]= self.listOfNbPointsOnEachIntervalForFinalDiscretization_Axis_kTotal[Axis_i]
    
##=====================================
##source : http://www.scientificpython.net/1/post/2012/04/clenshaw-curtis-quadrature.html    
    def clencurt(self,n1):
## Computes the Clenshaw Curtis nodes and weights """
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
        [xi, self.weights[Axis_i]] =  self.clencurt(self.listOfNbPointsOnEachIntervalForFinalDiscretization_Axis_kTotal[Axis_i])
        self.getDxiAfterTransf(xi, Axis_i) # Mettre à jours self.listOfTuckerGridNodes_Axis_k[Axis_i]
        self.weights[Axis_i] = (self.lx[Axis_i]/2)*self.weights[Axis_i]
        #self.listOfcoefW[Axis_i] = self.lx[Axis_i]/2    
        #self.listOfNbPointsOnEachIntervalForFinalDiscretization_Axis_kTotal[Axis_i]= self.listOfNbPointsOnEachIntervalForFinalDiscretization_Axis_kTotal[Axis_i] 
##=====================================    
    def discretizeClenshawCurtisBis(self, Axis_i) :#listOfExtremesForCouplage, listOfNbrPsForCouplage):
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
    def discretizeLegendreBis(self, Axis_i) :#listOfExtremesForCouplage, listOfNbrPsForCouplage):
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
    def transfIndexIntToArray(self, I, d, listOfNbPoints):   #### i = array([i_0, ...,i_(d-1)])
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
    # i = (i_0, ...,i_k-1,i_k+1,...i_d-1) concatene avec i_k pour avoir (i_0, ...,i_k-1,i_k,i_k+1,...i_d-1)    
    def concatenateIndexes(self, k, ik, i):
        return i.insert(k,ik)
##=====================================        
    def getXi(self, i):
        x = []
        #x = np.array(np.zeros(self.dim))
       
        for j in range(0, self.dimension):
            x.append(self.listOfTuckerGridNodes_Axis_k[j][i[j]])
            #x[j] = self.listOfTuckerGridNodes_Axis_k[j][i[j]]
        return x

##=====================================            
    def getCrossSectionOnSampleTest(self, crossSection, listOfDx,  listOfNbPoints):
    #correspondant avec la transformation d'indice dessus    
        nX = self.produitArr(listOfNbPoints[:self.dim] ) ### 
        fx = np.zeros(nX)      
        listOfxiOnNewDx = []
        for I in range(0, nX) :
            i = self.transfIndexIntToArray(I, self.dim ,listOfNbPoints)
            x = []
        #x = np.array(np.zeros(self.dim))       
            for j in range(0, self.dim):
                x.append(self.listOfTuckerGridNodes_Axis_k[j][i[j]])
           
            fx[I] = crossSection.evaluate(x)
        return fx
##=====================================    
    def getMatrixOfCrossSectionMkI(self, crossSection):
        #self.getNbPointsTotal()
    #correspondant avec la transformation d'indice dessus
        self.nbOfGroupedPointsY = self.NbPointsTotal/self.listOfNbPointsOnEachIntervalForFinalDiscretization_Axis_kTotal[self.Axis_k]
        listOfNbPointsNew = self.listOfNbPointsOnEachIntervalForFinalDiscretization_Axis_kTotal[:] ##Faire une copie, n'utilise pas "=" car il y aura probleme
        listOfNbPointsNew.remove(listOfNbPointsNew[self.Axis_k])
     
        M0 = np.zeros( self.NbPointsTotal)##self.listOfNbPointsOnEachIntervalForFinalDiscretization_Axis_kTotal[self.Axis_k]*self.nbOfGroupedPointsY)
        MkI = M0.reshape(self.listOfNbPointsOnEachIntervalForFinalDiscretization_Axis_kTotal[self.Axis_k], self.nbOfGroupedPointsY)
           
        for ik in range(0, int(self.listOfNbPointsOnEachIntervalForFinalDiscretization_Axis_kTotal[self.Axis_k])):
            for iI in range(0, int(self.nbOfGroupedPointsY)):
                i = self.transfIndexIntToArray(iI, self.dimension - 1, listOfNbPointsNew)
                i.insert(self.Axis_k,ik)
                x = self.getXi(i)
                MkI[ik,iI] = crossSection.evaluate(x)
        
        
        
       
        self.minOfCrossSectionValue = abs(MkI).min()
        self.meanOfCrossSectionValue = np.mean(MkI)
        self.meanOfSquaredCrossSectionValue = np.mean(np.power(MkI, 2))
        #print "MkI =", MkI
        #print "self.Axis_k = ",self.Axis_k 
        #print "self.minOfCrossSectionValue", self.minOfCrossSectionValue
        #raw_input()
        return MkI
##=====================================        
    def getMatrixOfFunctionMkI(self, f, args):
      
        self.nbOfGroupedPointsY = int(self.NbPointsTotal/self.listOfNbPointsOnEachIntervalForFinalDiscretization_Axis_kTotal[self.Axis_k])
        listOfNbPointsNew = self.listOfNbPointsOnEachIntervalForFinalDiscretization_Axis_kTotal[:] ##Faire une copie, n'utilise pas "=" car il y aura probleme
        listOfNbPointsNew.remove(listOfNbPointsNew[self.Axis_k])
        
        M0 = np.zeros(self.NbPointsTotal)
        MkI = M0.reshape(self.listOfNbPointsOnEachIntervalForFinalDiscretization_Axis_kTotal[self.Axis_k], self.nbOfGroupedPointsY)
        self.fxk_0 = [] ## Pour stocker les valeurs de f en tous les points de x_k et 
            ## les autres directions sont prise en premier point d'extrême dont  
            ## indice (0,0,...,x_k,0, 0, 0) => iI = 0
        print "==============="
        print "self.listOfNbPointsOnEachIntervalForFinalDiscretization_Axis_kTotal[self.Axis_k]", self.listOfNbPointsOnEachIntervalForFinalDiscretization_Axis_kTotal[self.Axis_k]
        print "==============="
        raw_input()
        for ik in range(0, int(self.listOfNbPointsOnEachIntervalForFinalDiscretization_Axis_kTotal[self.Axis_k])):
            print self.nbOfGroupedPointsY
            raw_input()
            for iI in range(0, int(self.nbOfGroupedPointsY)):
                print iI
                raw_input()
                i = self.transfIndexIntToArray(iI, self.dimension - 1, listOfNbPointsNew)
                
                i.insert(self.Axis_k,ik)             
                x = self.getXi(i)    
                MkI[ik,iI] =f(x,args)
        if iI == 0 :
            self.fxk_0.append(MkI[ik,iI])

        return MkI
##=====================================    
 
        
    def getProdWi(self,i): ## i = (i_0, ....,i_(k-1),i_(k+1),...i_(d-1)) sans i_k
        prodWi = 1.0        
        for j in range(0, self.dimension):
            prodWi = prodWi*self.weights[j][i[j]]
        return prodWi
##=====================================    
        
    def getWeights(self):
      
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
