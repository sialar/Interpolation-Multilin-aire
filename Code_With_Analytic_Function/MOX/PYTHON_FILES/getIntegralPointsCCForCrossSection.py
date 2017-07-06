# -*- coding: utf-8 -*-

import numpy as np
from numpy import linalg as LA
import scipy
from scipy import interpolate
from scipy.special.orthogonal import p_roots, t_roots

from scipy.fftpack import ifft

"""
  Ce script est pour pré-déterminer les points de Clenshaw-Curtis sur une interevalle donnée. 
  Ces points sont ensuite mis dans l'un de deux fichiers : mdb_irr.py ou mdb_crn.py pour lancer GAB
"""

"""    
  Pour dklib_5parameters sur le domain étendu 
  Burnup : 0, 80000
  fuel_temperature : 20, 1800 
  coolant_temperature : 20, 343
  coolant_density : 0.4, 1 
    b10_concentration :  0 2200
  Xe : 0.0001, '%', 200 '%'

  Nominal Values = {
  'fuel_temperature' : 550.0, 
  'coolant_density': 0.71692097187, 
  'b10_concentration': 3.97352914661e-06, 
  'xenon_level': 1 = 100 (%) 
  'coolant_temperature': 304.600006104}
"""

def getCoefAffine(xLOld, xUOld, xLNew, xUNew) :
    return [(xUOld - xLOld)/(xUNew - xLNew), (xLOld*xUNew -xUOld* xLNew)/(xUNew - xLNew)]

def transformeAffine(x, coef) :
    return coef[0]*x + coef[1]   
        
def getDxiAfterTransf(xi, Axis_i, x):
    coefX = getCoefAffine(x[Axis_i][0], x[Axis_i][1], -1.0, 1.0)
    domainDiscret_i = transformeAffine(xi, coefX)
    return domainDiscret_i
    
def discretizeLegendre(Axis_i,nbrP_i,x):
    
    [xi, weights_xi] = p_roots(nbrP_i)
    return getDxiAfterTransf(xi, Axis_i, x) 

def clencurt(n1):
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
    
def discretizeClenshawCurtis(Axis_i, nbrP_i, x):
    [xi, weights_xi] =  clencurt(nbrP_i)
    return getDxiAfterTransf(xi, Axis_i, x) # Mettre à jours self.domainDiscret[Axis_i]


if __name__ == "__main__":
    #methodIntegralNumerique = sys.argv[1]
    
    dimens = 1 
    
    x = np.zeros(dimens*2).reshape(dimens,2)   
    
    ### Rho :  nominal : 0.71692097187
    
    ###############################
    # A MODIFER (bornes de l'interevalle des CRN)
    ###############################
    
    x[0][0] = 150
    x[0][1] = 10000
    
    ###############################
    # FIN MODIFER
    ###############################
    
     
        
    method = "CC" #"CC"
    print "================\n"
    print " Domain =  ", x
    print "================\n" 
    
    nomFichier = "listOfPointsLegendesForDKlibs.txt"
    nbPOfAxis_k = int(raw_input("Donnez nbPOfAxis_k =  ")) 
     
    axisName = ["BURNUP", "fuel_temperature","coolant_density","b10", "xenon" ]
    with open(nomFichier, "a") as fichier :
        fichier.write("\n")
        fichier.write("******INFORMATION GENERAL*********\n\n")
        fichier.write("Information general of domain :\n")
        fichier.write("Dimensions : %s\n\n" %dimens)
        fichier.write("Points : %s\n\n" %method)
        #fichier.write("inforAxis(crossSection) %s" %inforAxis(crossSection))
        fichier.write("***************************************\n\n")
  
        for Axis_k in range(dimens):
            fichier.write("=========\n") 
            fichier.write("\t %s. Name of axis = %s\n"  %(Axis_k,axisName[Axis_k]))   
            fichier.write("Borders : [%s,%s]\n" %(x[Axis_k][0], x[Axis_k][1]) )
            #fichier.write("=========\n") 
            if method == "L" :
                domainDiscret_i = discretizeLegendre(Axis_k, nbPOfAxis_k, x)
            elif method == "CC" :
                domainDiscret_i = discretizeClenshawCurtis(Axis_k, nbPOfAxis_k, x)
            print "Axis_k = ", Axis_k
            print "domainDiscret_i ", domainDiscret_i
            fichier.write("==================\n")
            fichier.write( "Index - Value of " + axisName[Axis_k] + "\n")
            for i in range(len(domainDiscret_i)):
                
                fichier.write("%s.\t %s \n" %(i, domainDiscret_i[i]) )
    fichier.close() 

        
