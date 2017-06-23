# -*- coding: utf-8 -*-
import fctMultiVariates as fct
import numpy as np
import matplotlib.pyplot as plt
import copy
#import Gnuplot, Gnuplot.funcutils
if __name__ == "__main__":
    #argument1 = sys.argv[1]

    argsF = [0.5, 0.25, 0.5]
    csName = "fct3d"
    dimension = len(argsF)
    f = fct.function(argsF, csName, dimension)
    nbP = 200
    testValues_arr = np.linspace(0, 1.0,nbP)
    testValues = testValues_arr.tolist()
    #=========================
    #f.root_df_dx0()
    #=========================
    #print "x =", x
    #raw_input()
    #=========================
    concernedParameterName = "x1" #"x0"#"x2"
    concernedParameterValues = testValues
    concernedParameterIndex = 1#0 #2
    #=========================

    #x0 = 0.1
    #x1 = 0.2
    #x2 = 0.8

    list_x = [[0.0, 0.0, 0.0], [0.1, 0.1, 0.1],  [0.5, 0.1, 0.9], [0.5, 0.9, 0.1], [0.5, 0.5,0.5], [0.9, 0.9, 0.9], [1.0, 1.0, 1.0]] #, [1.5, 1.5, 1.5]]
    count = 0
    for x in list_x:
        x_origine = copy.copy(x)
        count = count + 1
        listValuesFctByConcernedParameter = []
        for x_parameter in concernedParameterValues:
            x[concernedParameterIndex] = x_parameter

            #print "value of feedback parameters = ", x
            #print "x_origine", x_origine
            #raw_input()
            valueFct = f.evaluate(x)
            listValuesFctByConcernedParameter.append(valueFct)

        fctValues_arr = np.asarray(listValuesFctByConcernedParameter)

        minFct = np.min(fctValues_arr)
        maxFct = np.max(fctValues_arr)

        nomFichierValPro = "MinMaxFctInFct" + str(concernedParameterName)+"_extendDomaine.txt"
        with open(nomFichierValPro, "a") as fichierVal :
            fichierVal.write("\n")
            fichierVal.write("*****************************************\n")
            fichierVal.write("======Fct = %s\n" %(f.name))
            fichierVal.write("x initial = %s \n" %(x))
            fichierVal.write("Varier en  %s \n" %(concernedParameterName))
            fichierVal.write("Min = %s \n" %(minFct) )
            fichierVal.write("Max = %s \n" %(maxFct) )
            fichierVal.write("*******************Fin de test **********************\n")
            fichierVal.write("\n")
            fichierVal.close()

        plt.plot(concernedParameterValues, fctValues_arr)
        plt.xlabel(concernedParameterName)
        plt.ylabel('f('+concernedParameterName +')')
        plt.title('Fct en ' + concernedParameterName + ", variant en x =" +str(x_origine))
        plt.grid(True)
        plt.savefig(f.name + "_" + concernedParameterName + "_" + str(count) + ".pdf")
        plt.show()
        raw_input()
