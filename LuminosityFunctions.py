#  infer number of pulsars in a globular cluster given a luminosity function, number of observed pulsars and minimum luminosity
from calendar import c
import math, csv
import numpy as np
from json import dump, load
from os import path
from scipy.special import erfc
import matplotlib.pyplot as plt
from matplotlib import rc

# set up the general model parameters of interest
a = -1.1 # mean of lognormal
b = 0.9  # sigma of lognormal
minCount = 0 #minimum number of pulsars in cluster to be eligible

font = {
        "family" : 'normal',
        "weight" : "bold",
        "size" : 14
    }

rc('font', **font)

if __name__ == "__main__":

    # locate data file and read clusters into dictionary
    DATAFILENAME = "Data03272024.csv"
    clusters = {}

    with open(DATAFILENAME, newline='') as datafile:
        fileReader = csv.reader(datafile)
        for row in fileReader: 
            clusterName = row[1]
            if clusterName == "47 Tucanae": #Force 47 Tuc to use NGC name for ReadClusterData.py processing
                clusterName = "NGC 104"
            lumValue = row[9]
            if row[0] != "Name" and lumValue: #exclude header row, only count pulsars with luminosities
                if clusterName not in clusters:
                    clusters[clusterName] = { #add new cluster to dictionary
                            "minLum" : lumValue,
                            "count" : 1,
                            "lumList" : [lumValue]
                        }
                else:
                    clusters[clusterName].update({ #update existing cluster
                            "minLum" : min(clusters[clusterName]["minLum"], lumValue),
                            "count" : clusters[clusterName]["count"] + 1
                        })
                    try:
                        clusters[clusterName]["lumList"].append(lumValue)
                    except AttributeError:
                        clusters[clusterName]["lumList"] = [lumValue]


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

            # Calculate the likelihood for each value of n
            results = [likelihood(n) for n in n_values]

            # Plot the probability distribution of N
            plt.clf() #clear figure
            plt.plot(n_values, results)
            plt.title(clusterName)
            plt.xlabel(f"Best estimate: {nhat} @ {round(likelihood(nhat)*100,2)}%\n{o} pulsars above {l} mJy kpc^2 ({round(fraction_above(l)*100,2)}% probed)")
            plt.yticks(ticks=[])
            plt.ylabel('Model likelihood')
            plt.subplots_adjust(bottom=0.2, left=0.05, right=0.95)
            #Add text summarizing assumptions and results
            #plt.gcf().text(0.5,0.03,f"{o} pulsars \nabove {l} mJy kpc^2\n({round(fraction_above(l)*100,2)}% probed)", horizontalalignment='center')
            #plt.gcf().text(0.02,0.05,f"Assuming log-normal of\nmean {a} and sigma {b}")
            plt.savefig(path.join(".","plots",f"{clusterName}.png"))
            #Add cluster to dictionary
            clusterInfo.update({
                    clusterName : {
                            "obsCount": clusters[clusterName]["count"],
                            "probableCount" : nhat,
                            "minLum" : clusters[clusterName]["minLum"],
                            "lumList" : clusters[clusterName]["lumList"],
                            "count" : len(clusters[clusterName]["lumList"])
                    }
                })
        
    with open("firstPass.dat", "w+") as out:
        dump(clusterInfo, out, indent = 4)