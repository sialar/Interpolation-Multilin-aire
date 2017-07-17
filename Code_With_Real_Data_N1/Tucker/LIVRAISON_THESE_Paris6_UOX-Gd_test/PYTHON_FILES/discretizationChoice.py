# -*- coding: utf-8 -*-
import copy

class discretizationChoice:
##=====================================
### This class is used to define respectively "d" Tucker grids, each Tucker grid has "d" discretizations for "d" axes, in which:
### + only a fine discrtization with the method CC ou CCbis (if having subdivision). Here CC = Clenshaw- Curtis
### + "d-1" coarse discretizations with 2 points, defined by method "T". Here T = Trapezoidal
##=====================================
    def __init__(self, listOfDomainBorders,listOfNbPointsOnEachIntervalForFinalDiscretization,\
            listOfSubdivisionDirection, listOfValuesForSubdivision, listOfNumberOfPointsForSubdivision):
        """
          
          Argumentes:
              "listOfDomainBorders "  [un dictionnaire]
                  + clé [une chain de caractère d'un entier] : ex '0', '1' ---> désigne l'axe
                  + Valeur [list of 2 flottant] : e.x [0,80000.0] ---> désigne min, max d'un paramètre
                  
              
              "listOfNbPointsOnEachIntervalForFinalDiscretization" [un dictionnaire]:
                  + clé [une chain de caractère d'un entier] : ex '0', '1' ---> désigne l'axe
                  + Valeur [une liste de "d" entiers - il y a une poisition est mis dans une liste ou liste de liste (si coupage) pour désigner l'axe principale dans la KL] : 
                        e.x [[3,9,9], 2, 2, 2, 2] ou [2,[5], 2,2,2]
               
              "listOfSubdivisionDirection" [un dictionnaire]:
                  + clé [une chain de caractère d'un entier] : ex '0', '1' ---> désigne l'axe sur lequel on a le coupage
                  + Valeur [un entier, valeur = 1] : il y a un coupage
              
              "listOfValuesForSubdivision"  [un dictionnaire]:
                  + clé [une chain de caractère d'un entier] : ex '0', '1' ---> désigne l'axe sur lequel on a le coupage
                  + Valeur [une liste de flottants] : les valeurs entre min, max qui sont les bornes de segments dans le coupage
               
               
              "listOfNumberOfPointsForSubdivision" [un dictionnaire]:
                  + clé [une chain de caractère d'un entier] : ex '0', '1' ---> désigne l'axe sur lequel on a le coupage
                  +  Valeur [une liste d'entiers] : e.x [3,9] (2 valeurs si on coupe en 3 segments, le nombre de points dans le 
                   dernier segment est en déduit de le nombre total - les nombre de points indiqué dans cette liste.

             
          Retourn:
            rien (constructeur)
       """      
              
              
              
        self.dimension = len(listOfDomainBorders)
        self.listOfDomainBorders = listOfDomainBorders
        self.listOfNbPointsOnEachIntervalForFinalDiscretization = listOfNbPointsOnEachIntervalForFinalDiscretization
        self.listOfSubdivisionDirection = listOfSubdivisionDirection
       
        self.listOfValuesForSubdivision = listOfValuesForSubdivision
        self.listOfNumberOfPointsForSubdivision = listOfNumberOfPointsForSubdivision
        """
            - listOfIndexesForSubdivision[un dictionnaire]: 
                  + Clé [une chain de caractère d'un entier] : ex '0', '1' ---> désigne l'axe sur lequel on a le coupage
                  + Valeur [une liste des entiers] : les indices des points correspondant aux extremes des segments dans le coupage
                  e.x : [0, 3, 12, 23] (4 valeurs -> 3 segments -> 2 points communs -> total : 23 -2 =21 points -> nbr de points par segments :
                  [3-0 , 12-3, 23-12-2] = [3,9,9] (à revoir)
        
            - listOfIntegralMethods [un dictionnaire]:
                  + Clé [une chaine de caractère d'un entier] : ex '0', '1' ---> désigne le nom de l'axe 
                  + Valeur [une liste de "d" composantes, chaque composante est une chaine de caractères] : e.x ['T', 'CC', 'T', 'T,'T']
            
        """
        self.listOfIndexesForSubdivision = {}
        self.listOfIntegralMethods = {}
        
        
    def getDiscretization(self):
        """
       Fonction qui définit quelle méthode d'intégration est utilisée pour chaque axe dans une grille de Tucker
       Répéter "d" fois pour définir "d" grilles
       Par défaut : l'axe étudié : "CC" ou "CCbis" (si coupage)
                    l'axe non étudié : "T" (Méthode des trapèzes avec 2 points extremes)
        """
       #===========15/2/2016=============== 
       ## Tbis --> T
        integralMethodDefaut = "T" #"Tbis"
        integralMethodDefautForSubdivision= "CCbis"
        integralMethodForStudiedDirection = "CC"
        
        
        listOfIntegralMethods = {}  


        for Axis_k in range(self.dimension):
            integralMethod_k = []
            for i in range(self.dimension):
                integralMethod_k.append(integralMethodDefaut)
           
                      
            integralMethod_k[Axis_k] = integralMethodForStudiedDirection             
       
            
            if str(Axis_k) in self.listOfSubdivisionDirection and self.listOfSubdivisionDirection[str(Axis_k)] == 1:
                integralMethod_k[Axis_k] = integralMethodDefautForSubdivision 
                
            self.listOfIntegralMethods[str(Axis_k)] = integralMethod_k 
            
            #====== Force to have only 2 integral points in the "T" method, otherwise, this will be the number of points that were calculated by dklib===
            for axis in range(self.dimension):
                
                if self.listOfIntegralMethods[str(Axis_k)][axis] == "T":          
                    self.listOfNbPointsOnEachIntervalForFinalDiscretization[str(Axis_k)][axis] = 2
      
        listOfIndexesPointsForSubdivision = {}
        
        for Axis_k in range(self.dimension):
            if  str(Axis_k) in self.listOfSubdivisionDirection and self.listOfSubdivisionDirection[str(Axis_k)] == 1:
                  for j in range(len(self.listOfValuesForSubdivision[str(Axis_k)])) :
                      val = self.listOfValuesForSubdivision[str(Axis_k)][j]
              
                      self.listOfDomainBorders[str(Axis_k)].insert(-1, val)### Position is counted by 1, 2, ..
                 
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
       