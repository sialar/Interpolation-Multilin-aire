# -*- coding: utf-8 -*-
# Ok
import os
import dkmodel
import dktools
import dkbase
import array
import Gnuplot
from numpy import *
import Gnuplot, Gnuplot.funcutils
import Symbols
###========01/03/2016====================
### This method is used to provide all cross-section types associated with a given cross-section name
### Exemple : input : cross-section name = macro_totale/ output : macro_totale0 (fast group), macro_totale1 (thermal group)
### Used with NEW VERSION :  DKZIP
### Data provided by GAB V2.3.1
###=======================================
from Tools import DKLibAdapter
##=====================================
def getCrossSections(dklibname, crossSectionName):
    print "This is the name of dklibs :", dklibname
 
    # Open the dklib in reading 
    #dklib = dktools.DKLib(dklibname)
    
    ###========01/03/2016====================
    ##=======NEW VERSION WITH DKZIP=======  
    dklib = DKLibAdapter.DKLibAdapter(dklibname)
    ##======================================= 
    
    print dklib

    # Browse DKLib
    pattern = dklib.getPatternNames()[0]

    print " pattern = " , pattern
    rodbank = dklib.getRodbankNames(pattern)[0]

    print " rodbank = " , rodbank

    #homogeneisation = dklib.getHomogenisationNames(pattern, rodbank)[1]
    homogeneisation = dklib.getHomogenisationNames(pattern, rodbank)[0]

    print " homogeneisation = " , homogeneisation

    material = dklib.getMaterialNames(pattern, rodbank, homogeneisation)[0]

    print " material = " , material


    path = dklib.getPath(pattern, rodbank, homogeneisation, material)

    print " path = \n", path

    # End of this path
    # Loading into the memory of the macroscopic cross sections and stored in the "domain" container 
    domain = dkmodel.Domain_Create("")
   
    #data = NeutronicLibraryManager.DKLibMaterialData(sections = {"macro": Symbols.DKLIB_MACRO_REACTIONS[:] + ("diffusion", 'absorption')},
                                                     #miscellaneousItems = ["b2", "diffusionCoefficient", "keff", "kinf"]
                                                     #)
    dklib.loadCrossSectionsForDomain(domain, path, dktools.LINEAR, ["macro"], [])
    dklib.loadMiscDataForDomain(domain, path, dktools.LINEAR)

    print "============= DOMAIN ============"
    print domain

    ## crossSection,diffusionCoefficient : have two fields corresponding to two energy lavels  
    # Recovery of the total section 
    if crossSectionName == "kinf" :
        kinf =  domain.getField("kinf")
        crossSection = kinf
        #print " kinf = ", kinf
        #raw_input()
    elif crossSectionName == "keff" :
        keff = domain.getField("keff")
        crossSection = keff
        #print " keff = ", keff
        #raw_input()
    elif crossSectionName == "b2":
        b2 = domain.getField("b2")
        crossSection = b2
        
    elif crossSectionName == "db2":
        b2 = domain.getField("b2")
        diffusionCoefficient= domain.getField("diffusionCoefficient")
        db2 =  b2 * diffusionCoefficient
        crossSection = db2
        #print "db2 = ", db2
        #raw_input()
    else:
        crossSection     = domain.getField(crossSectionName)
  
    return crossSection
###=====================================     
def getListOfNumberAgrs(dklibname, crossSectionName):
    
    crossSection = getCrossSections(dklibname, crossSectionName)
    listOfNbOfAgrs = []
   
    if isinstance(crossSection, dkmodel.TFEnergyXsLinear_Ptr): 
      listOfNbOfAgrs = [crossSection.getMesh().nCells()]
      #print "nb groupes = ", crossSection.getMesh().nCells()
    elif isinstance(crossSection, dkmodel.TFAnisotropyTFEnergyTFEnergyXsLinear_Ptr):
      anisotropyOrder = crossSection.getMesh().nCells()
      nbArrivalGroups = crossSection[0].getMesh().nCells()
      nbDepartGroups= crossSection[0][0].getMesh().nCells()
      listOfNbOfAgrs = [anisotropyOrder, nbArrivalGroups, nbDepartGroups]

    return listOfNbOfAgrs
###=====================================   
def getListOfArgs(dklibname, crossSectionName):
    listOfNbOfAgrs = getListOfNumberAgrs(dklibname, crossSectionName)
    listOfArgsOfCrossSection = []
   
    if len(listOfNbOfAgrs) == 1 :
        for energyLevel in range(listOfNbOfAgrs[0]):
            listOfArgsOfCrossSection.append([energyLevel])
            
    elif len(listOfNbOfAgrs) == 3 :            
            for anisotropy in range(listOfNbOfAgrs[0]):
                for arrivalGroup in range(listOfNbOfAgrs[1]):
                    for departGroup in range(listOfNbOfAgrs[2]):
                        ### case of macro_scattering_0_ArrivalGroup_DepartGroup
                        if anisotropy== 0 :
                            listOfArgsOfCrossSection.append([anisotropy, arrivalGroup, departGroup])          
            
    return listOfArgsOfCrossSection
    
###=====================================   
def getCrossSectionNameWithArg(crossSectionName, args):
    strArgs = ""
    for i in range(len(args)):
        strArgs = strArgs + str(args[i])
        
    finalCrossSectionName = crossSectionName + strArgs   
    return finalCrossSectionName
###=====================================  
def getfinalCrossSection(dklibname, crossSectionName, args):
    crossSection = getCrossSections(dklibname, crossSectionName)
    if crossSectionName=="kinf" or crossSectionName=="keff" or crossSectionName=="b2" or args == [] :
        return crossSection
   
    else :
        n = len(args)
        
        if crossSectionName == "macro_scattering" and n != 3 :          
            raise ValueError("The number of arguments for this cross section must be equal: 3 (anisotropy order, group of energy arrive, group of energy depart )")
        elif crossSectionName != "macro_scattering" and n != 1 :
            raise ValueError("The number of arguments this cross section must be equal: 1 (level of group of enerny)")
          
        if crossSectionName == "macro_scattering" :
            anisotropyOrder = int(args[0]) 
            arrivalGroup = int(args[1]) 
            departGroup = int(args[2])             
            finalCrossSection = crossSection[anisotropyOrder][arrivalGroup][departGroup] 
        else :
            levelGroupEnergy = int(args[0]) 
            finalCrossSection = crossSection[levelGroupEnergy]            
        return finalCrossSection
###=====================================      
def getFinalListOfCrossSections(dklibname, crossSectionName):
    listOfArgsOfCrossSection = getListOfArgs(dklibname, crossSectionName)
    listOfFinalCrossSections = {}
    
    if crossSectionName=="kinf" or crossSectionName=="keff" or crossSectionName=="b2" or listOfArgsOfCrossSection == []:
        cs_final = getCrossSections(dklibname, crossSectionName)
        listOfFinalCrossSections[crossSectionName] = cs_final
        return listOfFinalCrossSections
    else :
       
        #print "listOfArgsOfCrossSection =", listOfArgsOfCrossSection
        #raw_input()
        #listOfFinalCrossSections = []
        for indexArg in range(len(listOfArgsOfCrossSection)):
            print "indexArg =", indexArg
            args = listOfArgsOfCrossSection[indexArg]
            strArgs = ""
            for i in range(len(args)):
                strArgs = strAgrs + str(args[i])
            #print "agrs =", agrs
            #raw_input()
            finalCrossSectionName = crossSectionName + strArgs
            cs_final = getfinalCrossSection(dklibname, crossSectionName, args)
            listOfFinalCrossSections[finalCrossSectionName] = cs_final
        return listOfFinalCrossSections
    
###=====================================      
def getFeedbackMesh(dklibFileName):
    crossSectionName = "macro_totale"
    crossSectionTotal =  getCrossSections(dklibFileName, crossSectionName)
    crossSection =  crossSectionTotal[0] ### Prendre une section efficace, ici, niveau d'énergie 0
    
    feedbackMesh = crossSection.getMesh()
    ### Copier de getFeedbackMesh(dklibFileName, pattern) dans evalueAREVA
    bu_values = None
    density_values = None
    tc_values = None
    atb10_values = None
    cb_values = None
    ct_values = None
    xe_values = None
    for iaxis in range(feedbackMesh.nAxis()):
        param  = feedbackMesh.getAxisName(iaxis)
        fid    = Symbols.FEEDBACK_PARAMETERS.get(param)
        if fid is None:
            raise KeyError("Le parametre de CRN %s n'est pas correctement traite par COCAGNE actuellement." % param)
        values = list(feedbackMesh.getRectilinearAxisSteps(iaxis))
        print 'values ', values
        if fid == Symbols.BURNUP:
            bu_values = values
        elif fid == Symbols.COOLANT_DENSITY:
            density_values = values #[ conditionsNominales['coolant_density'] ]
        elif fid == Symbols.FUEL_TEMPERATURE:
            tc_values = values #[ conditionsNominales['fuel_temperature'] ]
        elif fid == Symbols.XENON_LEVEL:
            xe_values = values #[ conditionsNominales['xenon_level'] ]
        elif fid == Symbols.ATOMS_OF_BORON10:
            atb10_values = values #[ conditionsNominales['b10_concentration'] ] 
        elif fid == Symbols.BORON_CONCENTRATION:
            cb_values = values #[ conditionsNominales['boron_concentration']]
        elif fid == Symbols.COOLANT_TEMPERATURE: # bibliotheque pour accidents : comment traiter ?
            ct_values = values #[ conditionsNominales['coolant_temperature'] ]
        else:
            raise KeyError("Le parametre de CRN %s n'est pas traite actuellement par cet utilitaire de test." % param)

    if  not bu_values :
        raise SystemError("on n a pas trouve d'axe de burnup sur le maillage de CRN de %s" %dklibFileName)

    if atb10_values and cb_values:
        raise ValueError("on a a la fois ATOMS_OF_BORON10 et BORON_CONCENTRATION dans le maillage de CRN de %s" %dklibFileName)

    # on transforme le nombre d'atomes de bore 10 en CB (ppm)
    if atb10_values:
        cb_values = atb10_values
        # coolant_density_nominal = conditionsNominales['coolant_density']
        # # recuperation valeur nominale de densité
        # cb_values = []
        # alpha = Symbols.NAvogadro * Symbols.PB10 / Symbols.MBNAT() / 1.E+24 * 1.E-6
        # for atb10 in atb10_values:
        #     cb_values.append(atb10 / alpha / coolant_density_nominal )
        #print 'cb_values ', cb_values
    else:
        coolant_density_nominal = None
   
    # parameter orders : for 6 parameters (case if)/ for 5 parameters (case else)
    if ct_values == None:
        return bu_values,  tc_values, density_values, cb_values, xe_values , ct_values
    else:
        return bu_values, tc_values, ct_values, density_values,  cb_values,  xe_values

            
###=====================================          
def prepareDomain(dklibName, material, data):
    domain = dkmodel.Domain_Create("%s" % material)
    dklib = dktools.DKLib(dklibName)
    Geometry._storeDKLibInfosInMaterial(domain, dklibName, material, interpolationMethod=dktools.LINEAR)
    loadFromDKLib.loadMaterialDataFromDKLib(dklib, domain, data)
    dklib.close()
    return domain     

