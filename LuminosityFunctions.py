#  infer number of pulsars in a globular cluster given a luminosity function, number of observed pulsars and minimum luminosity
import math
import csv
import numpy as np
from scipy.special import erfc
import matplotlib.pyplot as plt

# set up the general model parameters of interest
a = -1.1 # mean of lognormal
b = 0.9  # sigma of lognormal
minCount = 0 #minimum number of pulsars in cluster to be eligible

# locate data file and read clusters into dictionary
DATAFILENAME = "Data1102023.csv"
clusters = {}

with open(DATAFILENAME, newline='') as datafile:
    fileReader = csv.reader(datafile)
    for row in fileReader: 
        clusterName = row[1]
        lumValue = row[9]
        if row[0] != "Name" and lumValue: #exclude header row, only count pulsars with luminosities
            if clusterName not in clusters:
                clusters[clusterName] = { #add new cluster to dictionary
                        "minLum" : lumValue,
                        "count" : 1
                    }
            else:
                clusters[clusterName].update({ #update existing cluster
                        "minLum" : min(clusters[clusterName]["minLum"], lumValue),
                        "count" : clusters[clusterName]["count"] + 1
                    })


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

for clusterName in clusters.keys(): #loop through every cluster in dataset
    if clusters[clusterName]["count"] >= minCount:
        #read in cluster-specific values
        o = int(clusters[clusterName]["count"])   # number of observed pulsars in the cluster
        l = float(clusters[clusterName]["minLum"]) # minimum luminosity probed (mJy kpc^2)

        # find the most probable value and decide on range of pulsars to plot over
        nhat = int(o / fraction_above(l))
        n_values = np.arange(o, 2*nhat)  # Adjust the range as needed

        # Calculate the likelihood for each value of n
        results = [likelihood(n) for n in n_values]

        # Plot the probability distribution of N
        plt.plot(n_values, results)
        plt.xlabel(f"Estimate of true number of pulsars in {clusterName}")
        plt.ylabel('Model likelihood')
        plt.show()

        # Summarize the assumptions and results of this analysis
        print (f"In {clusterName}, we observe",o,"pulsars above",l,"mJy kpc^2")
        print ("Assuming a log-normal LF with mean",a,"and sigma",b)
        print ("In this case, the fraction of the LF probed is",fraction_above(l))
        print ("Best estimate for this cluster is",nhat,"pulsars")
        print ("Maximum likelihood in this case is",likelihood(nhat))