from numpy.random import default_rng
from json import load

def generatePop(mu, sigma, count, mini):
    """
    mu and sigma are properties of the lognormal luminosity function that we assume all globular cluster pulsars follow.
    count is the number of pulsars above luminosity mini.
    returns list of luminosities
    """
    rng = default_rng()
    out = []
    while count > 0:
        newPsr = rng.normal(loc=mu, scale=sigma)
        if 10**newPsr > mini:
            count -= 1
            out.append(newPsr)
    return out

def countAbove(pop, mini):
    count = 0
    for lum in pop:
        if lum >= mini:
            count = count + 1 
    return count