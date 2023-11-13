#  infer number of pulsars in a globular cluster given a luminosity function, number of observed pulsars and minimum luminosity
import math
import numpy as np
from scipy.special import erfc

# set up the parameters of interest
o = 29   # number of observed pulsars in the cluster
l = 0.6 # minimum luminosity probed (mJy kpc^2)
a = -1.1 # mean of lognormal
b = 0.9  # sigma of lognormal

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

# find the most probable value and decide on range of pulsars to plot over
nhat = int(o / fraction_above(l))
n_values = np.arange(o, 2*nhat)  # Adjust the range as needed

# Calculate the likelihood for each value of n
results = [likelihood(n) for n in n_values]

# Plot the probability distribution of N
import matplotlib.pyplot as plt

plt.plot(n_values, results)
plt.xlabel('Estimate of true number of pulsars in cluster')
plt.ylabel('Model likelihood')
plt.show()

# Summarize the assumptions and results of this analysis
print ("In this cluster, we observe",o,"pulsars above",l,"mJy kpc^2")
print ("Assuming a log-normal LF with mean",a,"and sigma",b)
print ("In this case, the fraction of the LF probed is",fraction_above(l))
print ("Best estimate for this cluster is",nhat,"pulsars")
print ("Maximum likelihood in this case is",likelihood(nhat))