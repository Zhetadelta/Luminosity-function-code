#  infer number of pulsars in a globular cluster given a luminosity function, number of observed pulsars and minimum luminosity
from calendar import c
import math, csv
import numpy as np
from json import dump, load
from os import path
from scipy.special import erfc
import matplotlib.pyplot as plt

# set up the general model parameters of interest
a = -1.1 # mean of lognormal
b = 0.9  # sigma of lognormal
minCount = 0 #minimum number of pulsars in cluster to be eligible

# locate data file and read clusters into dictionary
DATAFILENAME = "Data03272024.csv"

with open("firstpass.dat", "r") as firstpass:
    clusters = load(firstpass)

# with open(DATAFILENAME, newline='') as datafile:
#     fileReader = csv.reader(datafile)
#     for row in fileReader: 
#         clusterName = row[1]
#         if clusterName == "47 Tucanae": #Force 47 Tuc to use NGC name for ReadClusterData.py processing
#             clusterName = "NGC 104"
#         lumValue = row[9]
#         if row[0] != "Name" and lumValue: #exclude header row, only count pulsars with luminosities
#             if clusterName not in clusters:
#                 clusters[clusterName] = { #add new cluster to dictionary
#                         "minLum" : lumValue,
#                         "count" : 1
#                     }
#             else:
#                 clusters[clusterName].update({ #update existing cluster
#                         "minLum" : min(clusters[clusterName]["minLum"], lumValue),
#                         "count" : clusters[clusterName]["count"] + 1
#                     })
#                 try:
#                     clusters[clusterName]["lumList"].append(lumValue)
#                 except KeyError:
#                     clusters[clusterName]["lumList"] = [lumValue]


###
### FUNCTION DEFINITIONS
###
# Function to calculate fraction of the luminosity function above lmin
def fraction_above(lmin):
    return 0.5 * erfc((math.log10(lmin) - a) / b / math.sqrt(2))

# Function to calculate the likelihood using the Binomial formula
def likelihood(n):
    f = fraction_above(l)
    numerator = math.factorial(n)
    denominator = math.factorial(o) * math.factorial(n - o)
    result = (numerator / denominator) * (f ** o) * ((1 - f) ** (n - o))
    return result

###
### PROCESSING
###

clusterInfo = {}

for clusterName in clusters.keys(): #loop through every cluster in dataset
     
    if clusters[clusterName]["count"] >= minCount:
        #read in cluster-specific values
        o = int(clusters[clusterName]["count"])   # number of observed pulsars in the cluster
        l = float(clusters[clusterName]["minLum"]) # minimum luminosity probed (mJy kpc^2)

        # find the most probable value and decide on range of pulsars to plot over
        nhat = int(o / fraction_above(l))
        n_values = np.arange(o, 2*nhat)  # Adjust the range as needed

        # Calculate the likelihood for each value of n, convert to normalized cumulative sum
        results = [likelihood(n) for n in n_values]
        results = np.cumsum(results)
        normResults = results / results[-1]

        # Plot the cumulative probability distribution of N
        plt.clf() #clear figure
        plt.plot(n_values, normResults)
        plt.title(clusterName)
        plt.xlabel(f"Estimate of true number of pulsars \nBest estimate: {nhat} @ {round(likelihood(nhat)*100,2)}%")
        plt.ylabel('Cumulative likelihood of lower than estimate')
        plt.subplots_adjust(bottom=0.25)
        #Add text summarizing assumptions and results
        plt.gcf().text(0.98,0.03,f"{o} pulsars \nabove {l} mJy kpc^2\n({round(fraction_above(l)*100,2)}% probed)", horizontalalignment='right')
        plt.gcf().text(0.02,0.05,f"Assuming log-normal of\nmean {a} and sigma {b}")
        plt.savefig(path.join(".","plots","cumulative",f"{clusterName}.png"))
        
        #Strip lowest and highest 2.5% to get 95% interval
        pairedPoints = list(zip(n_values, normResults))
        intervalMin = None
        intervalMax = None
        for i in range(len(pairedPoints)):
            if pairedPoints[i][1] > .025 and intervalMin is None:
                intervalMin = pairedPoints[i][0]
            elif pairedPoints[i][1] > .975 and intervalMax is None:
                intervalMax = pairedPoints[i-1][0]
                break


        #Add cluster to dictionary
        clusterInfo.update({
                clusterName : {
                        "obsCount": clusters[clusterName]["count"],
                        "probableCount" : nhat,
                        "95min" : int(intervalMin),
                        "95max" : int(intervalMax),
                        "minLum" : l,
                        "lumList" : clusters[clusterName]["lumList"]
                }
            })
        
with open("firstPass.dat", "w+") as out:
    dump(clusterInfo, out, indent = 4)