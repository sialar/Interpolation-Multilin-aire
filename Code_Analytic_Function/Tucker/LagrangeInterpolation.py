# -*- coding: utf-8 -*-

##===================================================  
### For a fct in the class of LagrangePolynomial
### fct is a list of sub-functions
### + Either: fct = [f]
### + Or fct = [f1, f2, ...,fn]
##===================================================   
def getInterpolation(fct,x):
    
    n = len(fct)
    ##===================================================
    if n == 1:
        return fct[0].getInterpolation(x)
        
        
    if x < min(fct[0].axisValues): # or x > min(fct[n-1].axisValues):
        return fct[0].getInterpolation(x)
        
        
    if x > max(fct[n-1].axisValues):
        return fct[n-1].getInterpolation(x)        
       
    ##===================================================      
    if n > 1 and (min(fct[0].axisValues) <= x and x <= max(fct[n-1].axisValues)):
        j = 0
        while j <= n-1:
             
            if min(fct[j].axisValues) <= x and x <= max(fct[j].axisValues):
                
                return fct[j].getInterpolation(x)
                break
            else :
                j = j +1
##===================================================                  
def getInterpolationArr(fct,x):
    listOfValues = []
    for xi in x:
        listOfValues.append(getInterpolation(fct,xi))
    return listOfValues
##===================================================    