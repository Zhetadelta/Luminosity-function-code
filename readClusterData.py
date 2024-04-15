from os import path
from json import load
from numpy import arange, mean, corrcoef, sqrt, log10, pi, linspace
import matplotlib.pyplot as plt
from matplotlib import rc
from simulation import *

MAKE_CLUSTER_PLOTS_FLAG = False #change this to regen plots
SIMULATION_ROUNDS = 0 #number of rounds to run simulation
SIMULATION_MU = -1.1
SIMULATION_SIGMA = 0.9
MIN_OBSERVATIONS = 2 #adjust to eliminate low-data clusters


font = {
        "weight" : "bold",
        "size" : 14
    }

rc('font', **font)


clusterDic = {}
with open("firstPass.dat") as dicFile:
    #load dictionary from LuminosityFunctions.py into code
    rawClusterDic = load(dicFile)
    for cluster in rawClusterDic:
        if rawClusterDic[cluster]["obsCount"] >= MIN_OBSERVATIONS:
            clusterDic[cluster] = rawClusterDic[cluster]
    
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
        ("absMag", "Pulsar Count vs Absolute Visual Magnitude", "Absolute Magnitude"),
        ("encounterRate", "Pulsar Count vs Encounter Rate", "Encounter Rate")
    ]

    
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
    
    #read encounter data in too
    encounterList = None
    with open("clusterEncounters.dat", "r") as encounterFile:
        encounterList = encounterFile.readlines()
        
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
                
        #attempt to estimate total cluster mass via color + luminosity
        clusterColor = float(dataFileTableII[rowNumMetal][60:67])
        avgTemp = 4600.*(1/(0.92*clusterColor+1.7)+1/(0.92*clusterColor+0.62)) #Ballesteros' formula
        #R = M^0.8; L = 4*pi*R^2*sigma*T^4; L = M^4
        avgMass = (avgTemp/5770)**(4/3.2) #algebra in notebook; in solar masses
        avgLum = avgMass**4
        clusterLum = 10**(0.4*(4.85-clusterDic[clusterName]["absMag"]))
        clusterDic[clusterName]["totalMass"] = clusterLum/avgLum*avgMass #VERY rough; solar masses 
        clusterDic[clusterName]["numPerMass"] = clusterDic[clusterName]["probableCount"]/clusterDic[clusterName]["totalMass"]
        #set a core-collapsed boolean, just in case
        clusterDic[clusterName]["coreCollapsed"] = "c" in dataFileTableIII[rowNum][56:59]
        
        #encounter rate stuff
        #if there's no correlation here, it won't be anywhere
        #convert cluster name
        encName = "".join(clusterName.lower().split(" "))+" "
        encIndex = None
        for i in range(0,len(encounterList)):
            if encName in encounterList[i]:
                encIndex = i
                break
        outliers = None
        with open("encRateOutliers.dat", "r") as outliersFile:
            outliers = [line.strip()+" " for line in outliersFile.readlines()]
        if encIndex is None:
            print(f"Cluster {clusterName} not in external dataset")
        #trying a thing
        #elif encName in outliers:
            #print(f"Cluster {clusterName} is excluded from encounter analysis")
        else:
            encRateString = encounterList[encIndex][58:66]
            encRateMan, encRateExp = float(encRateString.split("E")[0]), int(encRateString.split("E")[1])
            encRate = encRateMan*(10**encRateExp)
            clusterDic[clusterName]["encounterRate"] = encRate
        
        
        
            

#We've got our dictionary! Still needs some massaging though. Let's calculate physical radius of cluster cores first.
#sum all pulsar counts while we're here

totalCount = 0
for clusterName in clusterDic.keys():
    clusterData = clusterDic[clusterName]
    totalCount += clusterData["probableCount"]
    d = clusterData["distanceKPc"] * 1000 #parsecs
    r_o = clusterData["coreObsRadiusArcMin"] * 0.000290888 #radians
    clusterData.update({
            "coreRadiusPc" : d*r_o #small-angle approximation is our friend
        })
 



properties = {}

for valueName, title, yLable in GENERATE_PLOTS:
    plt.clf()
    xValues = []
    xValuesRatios = []
    yValues = []
    xErrorMin = []
    xErrorMax = []
    for clusterName in clusterDic.keys():
        if clusterDic[clusterName]["obsCount"] >= MIN_OBSERVATIONS:
            xValues.append(clusterDic[clusterName]["probableCount"])
            xValuesRatios.append(clusterDic[clusterName]["numPerMass"])
        
            #two clusters dont have encounter rates
            #handle that here
            removed = False
            try:
                yValues.append(clusterDic[clusterName][valueName])
            except KeyError:
                if valueName == "encounterRate":
                    xValues.pop() #keep lists aligned
                    xValuesRatios.pop()
                    removed = True
                else:
                    raise KeyError("oops")
            
            if not removed:
                xErrorMin.append(clusterDic[clusterName]["95min"])
                xErrorMax.append(clusterDic[clusterName]["95max"])
        
    properties[valueName] = list(zip(xValues, yValues))
        
    if MAKE_CLUSTER_PLOTS_FLAG: #only if flag is set
        xError = [xErrorMin, xErrorMax]
        #plt.scatter(xValues, yValues)
        plt.yscale("log")
        if valueName == "encounterRate":
            plt.xscale("log")
        plt.title(title)
        plt.ylabel("Most probable count of pulsars")
        plt.xlabel(yLable)
        plt.errorbar(yValues, xValues, yerr=xError, fmt='or', capsize=0) #swap the x and y axes the messy way
        plt.savefig(path.join(".","plots","properties",f"{valueName}.png"))
        plt.clf()
        
        #do it again with the ratio
        xError = [xErrorMin, xErrorMax]
        #plt.scatter(xValues, yValues)
        #plt.yscale("log")
        #plt.xscale("log")
        plt.title(title+" (per solar mass)")
        plt.ylabel("Pulsar count per solar mass")
        plt.xlabel(yLable)
        plt.errorbar(yValues, xValuesRatios, fmt='or', capsize=0) #swap the x and y axes the messy way
        plt.savefig(path.join(".","plots","properties",f"{valueName}perSolarMass.png"))
        plt.clf()
        
        #do it again with the ratio (log edition)
        xError = [xErrorMin, xErrorMax]
        #plt.scatter(xValues, yValues)
        plt.yscale("log")
        #plt.xscale("log")
        plt.title(title+" (per solar mass)")
        plt.ylabel("Pulsar count per solar mass")
        plt.xlabel(yLable)
        plt.errorbar(yValues, xValuesRatios, fmt='or', capsize=0) #swap the x and y axes the messy way
        plt.savefig(path.join(".","plots","properties",f"{valueName}perSolarMassLog.png"))
        plt.clf()
        
    plt.clf()

#generate histogram of real data before simulation
binCount = 40

total = []
for clusterName, clusterProps in clusterDic.items():
    total.extend([log10(float(lum)) for lum in clusterProps["lumList"]])
    
plt.hist(total, bins=binCount)
plt.ylabel("N(L)")
plt.xlabel("log L")
plt.title("Observed GC PSR Population by Luminosity")
plt.yticks(ticks=[])
plt.show()

plt.clf()
total.sort()
total.reverse()
lumValues = []
lumRank = []
for i in range(len(total)):
    lumValues.append(total[i])
    lumRank.append(log10(i))

plt.scatter(lumValues, lumRank)
plt.ylabel("log N(L > L_0)")
plt.xlabel("log L_0")
plt.yticks(ticks=[])
plt.title("Observed GC PSR Cumulative Count")
plt.show()

#now a histogram of cluster populations
plt.clf()
counts = [int(clusterDic[cluster]["obsCount"]) for cluster in clusterDic]
maxcount = max(counts)
plt.hist([clusterDic[cluster]["obsCount"] for cluster in clusterDic], bins=list(range(maxcount+1)))
plt.show()

#print sum-of-minima and sum-of-maxima information for uncertainties
small = 0
likely = 0
big = 0
for cluster in clusterDic:
    small += clusterDic[cluster]["95min"]
    likely += clusterDic[cluster]["probableCount"]
    big += clusterDic[cluster]["95max"]
print(f"Minima: {small}, Likely: {likely}, Big: {big}")

simTotal = []
for i in range(SIMULATION_ROUNDS): #simulation stuff
    thisSim = []
    for clusterName, clusterProps in clusterDic.items():        
        thisSim.extend(generatePop(SIMULATION_MU, SIMULATION_SIGMA, clusterProps["obsCount"], clusterProps["minLum"]))
    simTotal.extend(thisSim)

if SIMULATION_ROUNDS > 0:
    binCount = 40
    
    plt.hist(simTotal, bins=binCount)
    plt.ylabel("N(L)")
    plt.xlabel("log L")
    plt.title("Simulated GC PSR Population by Luminosity")
    plt.yticks(ticks=[])
    plt.show()

    plt.clf()
    simTotal.sort()
    simTotal.reverse()
    lumValues = []
    lumRank = []
    for i in range(len(simTotal)):
        lumValues.append(simTotal[i])
        lumRank.append(log10(i))

    plt.scatter(lumValues, lumRank)
    plt.ylabel("log N(L > L_0)")
    plt.xlabel("log L_0")
    plt.yticks(ticks=[])
    plt.title("Simulated GC PSR Cumulative Count")
    plt.show()

    
#find corrrelation coefficient via Pearson product-moment
#use numpy.corrcoef to check implementation

coeffs = {}
datapoints = len(clusterDic)

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
    
    coeffRange = linspace(-1, 1, 2000)
    coeffPDF = []
    for i in coeffRange:
        coeffPDF.append( 
                (1-i**2)**((datapoints-1)/2)/((1-i*coeff)**(datapoints-1.5))*(1+(1/(datapoints-0.5))*((1+i*coeff)/(8)))
            )
    plt.clf()
    plt.plot(coeffRange, coeffPDF)
    plt.yticks(ticks=[])
    plt.xticks(ticks=[-1, 0, coeff, 1], labels=["-1", "0", "ρ", "1"])
    plt.subplots_adjust(bottom=0.4)
    print(f"{prop}")
    #plt.show()
    


with open('correlation_coefficients.txt', "+w") as file:
    file.write("Property; This implementation; Numpy-calculated\n")
    for prop in coeffs.keys():
        file.write(f"{prop.ljust(23)} |{str(round(coeffs[prop]['mine'], 3)).ljust(6)}|{str(round(coeffs[prop]['numpy'], 3)).ljust(5)}\n")
