# -*- coding: utf-8 -*-
import copy
class discretizationChoice:
    def __init__(self, listOfDomainBorders,listOfNbPointsOnEachIntervalForFinalDiscretization,\
            listOfSubdivisionDirection, listOfValuesForSubdivision, listOfNumberOfPointsForSubdivision):
        self.dimension = len(listOfDomainBorders)
        self.listOfDomainBorders = listOfDomainBorders
        self.listOfNbPointsOnEachIntervalForFinalDiscretization = listOfNbPointsOnEachIntervalForFinalDiscretization
        self.listOfSubdivisionDirection = listOfSubdivisionDirection
       
        self.listOfValuesForSubdivision = listOfValuesForSubdivision
        self.listOfNumberOfPointsForSubdivision = listOfNumberOfPointsForSubdivision
        self.listOfIndexesForSubdivision = {}
        self.listOfMethodsIntegrals = {}
        ### listOfNombresOfPointsForSubdivision on peut en déduit de dklib de couplage + listOfValuesForSubdivision
        
    def getDiscretization(self):
       #===========15/2/2016=============== 
       ## Tbis --> T
        methodIntegralDefaut = "T" #"Tbis"
        methodIntegralDefautForSubdivision= "CCbis"
        methodIntegralForStudiedDirection = "CC"#"CCbis"
        
        
        listOfMethodsIntegrals = {}  # list of list, sous liste a dimens composantes
        
        for Axis_k in range(self.dimension):
            methodIntegral_k = []
            for i in range(self.dimension):
                methodIntegral_k.append(methodIntegralDefaut)
           
                      
            methodIntegral_k[Axis_k] = methodIntegralForStudiedDirection             
       
            
            if str(Axis_k) in self.listOfSubdivisionDirection and self.listOfSubdivisionDirection[str(Axis_k)] == 1:
                methodIntegral_k[Axis_k] = methodIntegralDefautForSubdivision 
                
            self.listOfMethodsIntegrals[str(Axis_k)] = methodIntegral_k 
            #===========15/2/2016===============
            #====== Forcer 2 points d'intégrations dans la méthode "T", sinon, ca sera le nbr de points que l'on a calculé par dklib===
            for axis in range(self.dimension):
                
                if self.listOfMethodsIntegrals[str(Axis_k)][axis] == "T":          
                    self.listOfNbPointsOnEachIntervalForFinalDiscretization[str(Axis_k)][axis] = 2

        
        listOfIndexesPointsForSubdivision = {}
        
        for Axis_k in range(self.dimension):
            if  str(Axis_k) in self.listOfSubdivisionDirection and self.listOfSubdivisionDirection[str(Axis_k)] == 1:
                  for j in range(len(self.listOfValuesForSubdivision[str(Axis_k)])) :
                      val = self.listOfValuesForSubdivision[str(Axis_k)][j]
              
                      self.listOfDomainBorders[str(Axis_k)].insert(-1, val)### Position est compté 1, 2, ..
                 
                      currentNbOfPoint_Axis_k =  self.listOfNbPointsOnEachIntervalForFinalDiscretization[str(Axis_k)][Axis_k][-1]
                      
                      newNbOfPoint_Axis_k = currentNbOfPoint_Axis_k + 1 ### the commun point is considered as two difference points
                     
                      nbOfPointsInLastInterval_Axis_k =  newNbOfPoint_Axis_k \
                                                        -  self.listOfNumberOfPointsForSubdivision[str(Axis_k)][j]
                      ### Update the new discretization by coupling
                
                      self.listOfNbPointsOnEachIntervalForFinalDiscretization[str(Axis_k)][Axis_k].remove(self.listOfNbPointsOnEachIntervalForFinalDiscretization[str(Axis_k)][Axis_k][-1])                          
                      self.listOfNbPointsOnEachIntervalForFinalDiscretization[str(Axis_k)][Axis_k].append(self.listOfNumberOfPointsForSubdivision[str(Axis_k)][j])
                      self.listOfNbPointsOnEachIntervalForFinalDiscretization[str(Axis_k)][Axis_k].append(nbOfPointsInLastInterval_Axis_k)
             
                  
                  self.listOfIndexesForSubdivision[str(Axis_k)] =  self.listOfNbPointsOnEachIntervalForFinalDiscretization[str(Axis_k)][Axis_k][:] ### [:] for copy
                  self.listOfIndexesForSubdivision[str(Axis_k)].insert(0,0)                   
                      
                  for p in range(1,len(self.listOfIndexesForSubdivision[str(Axis_k)])):                 
                      self.listOfIndexesForSubdivision[str(Axis_k)][p] = self.listOfIndexesForSubdivision[str(Axis_k)][p-1] \
                                                                + self.listOfIndexesForSubdivision[str(Axis_k)][p]
       