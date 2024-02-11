from os import path
from json import load
from numpy import arange, mean, corrcoef, sqrt
import matplotlib.pyplot as plt

clusterDic = {}
with open("firstPass.dat") as dicFile:
    #load dictionary from LuminosityFunctions.py into code
    clusterDic = load(dicFile)
    
#hard-coded column name, start of slice, end of slice tuples for table III properties
DATA_FILE_COLUMNS = [
        ("radialVelocity", 13, 20), #this should have no correlation, good to check
        ("centralConcentration", 49, 55),
        ("coreObsRadiusArcMin", 59, 65), #in arcmin, need cluster distance
        ("luminosityDensity", 80, 87),
        ("velocityDispersion", 36, 42),
    ]

GENERATE_PLOTS = [ #code name, table header tuples
        ("radialVelocity", "Pulsar Count vs Heliocentric Radial Velocity", "Radial Velocity (km s^-1)"),
        ("centralConcentration", "Pulsar Count vs King-model Central Concentration", "Central Concentration"),
        ("luminosityDensity", "Pulsar Count vs Central Luminosity Density", "Luminosity Density, log(L_☉ pc^-3)"),
        ("coreRadiusPc", "Pulsar Count vs Core Radius", "Core Radius (pc)"),
        ("velocityDispersion", "Pulsar Count vs Central Velocity Dispersion", "Central Velocity Dispersion, log(km/s)"),
        ("metallicity", "Pulsar Count vs Cluster Metallicity", "Metallicity (Fe/H)"),
        ("absMag", "Pulsar Count vs Absolute Visual Magnitude", "Absolute Magnitude")
    ]

MIN_OBSERVATIONS = 0 #adjust to eliminate low-data clusters
    
with open("clusterData.dat") as dataFile:
    #here we go. First, load into a list:
    dataFileMainList = dataFile.readlines()
    #Restrict to table III (velocities and structural parameters) and table I (positional data)
    dataFileTableIII = dataFileMainList[433:590]
    dataFileTableII = dataFileMainList[252:410]
    dataFileTableI = dataFileMainList[72:230]
    #Extract cluster names and validate cluster names in generated dictionary
    idListTableIII = [row[:12].strip() for row in dataFileTableIII]
    idListTableII = [row[:12].strip() for row in dataFileTableII]
    idListTableI = [row[:12].strip() for row in dataFileTableI]
    for clusterName in clusterDic.keys():
        if clusterName not in idListTableIII:
            raise ValueError(f"{clusterName} is invalid. Perhaps it has a different name?")
        #Retrieve cluster distance from table I
        rowNumI = idListTableI.index(clusterName)
        clusterDic[clusterName]["distanceKPc"] = float(dataFileTableI[rowNumI][67:74]) 
        rowNumMetal = idListTableII.index(clusterName)
        clusterDic[clusterName]["metallicity"] = float(dataFileTableII[rowNumMetal][12:19])
        clusterDic[clusterName]["absMag"] = float(dataFileTableII[rowNumMetal][47:54])
        #Go through desired properties defined in DATA_FILE_COLUMNS and add to clusterDic
        rowNum = idListTableIII.index(clusterName)
        for prop, start, end in DATA_FILE_COLUMNS:
            try:
                clusterDic[clusterName].update({
                        prop : float(dataFileTableIII[rowNum][start:end])
                    })
            except ValueError:
                clusterDic[clusterName].update({
                        prop: 0
                    })
        #set a core-collapsed boolean, just in case
        clusterDic[clusterName]["coreCollapsed"] = "c" in dataFileTableIII[rowNum][56:59]

#We've got our dictionary! Still needs some massaging though. Let's calculate physical radius of cluster cores first.
for clusterName in clusterDic.keys():
    clusterData = clusterDic[clusterName]
    d = clusterData["distanceKPc"] * 1000 #parsecs
    r_o = clusterData["coreObsRadiusArcMin"] * 0.000290888 #radians
    clusterData.update({
            "coreRadiusPc" : d*r_o #small-angle approximation is our friend
        })

properties = {}
    
for valueName, title, yLable in GENERATE_PLOTS:
    xValues = []
    yValues = []
    for clusterName in clusterDic.keys():
        yValues.append(clusterDic[clusterName]["probableCount"])
        xValues.append(clusterDic[clusterName][valueName])
        
    properties[valueName] = list(zip(xValues, yValues))

#find corrrelation coefficient via Pearson product-moment
#use numpy.corrcoef to check implementation

coeffs = {}

for prop in properties.keys():
    valueTuples = properties[prop]
    propList, countList = ([x[0] for x in valueTuples], [x[1] for x in valueTuples])
    propMean, countMean = (mean(propList), mean(countList))
    num = 0
    sumPropSquared = 0
    sumCountSquared = 0
    for i in range(len(propList)):
        num += (propList[i] - propMean)*(countList[i] - countMean)
        sumPropSquared += (propList[i] - propMean)**2
        sumCountSquared += (countList[i] - countMean)**2
    coeff = num/sqrt(sumPropSquared * sumCountSquared)

    npcoeff = corrcoef(propList, y=countList)[0][1]
    coeffs[prop] = {
            "mine" : coeff,
            "numpy" : npcoeff
        }

with open('correlation_coefficients.txt', "+w") as file:
    file.write("Property; This implementation; Numpy-calculated\n")
    for prop in coeffs.keys():
        file.write(f"{prop.ljust(23)} |{str(round(coeffs[prop]['mine'], 3)).ljust(6)}|{str(round(coeffs[prop]['numpy'], 3)).ljust(5)}\n")