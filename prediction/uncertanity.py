from calibrate import misfit
import numpy as np
from numpy.core.fromnumeric import mean
from model_fitting.cm_solver_functions import *
from data.interpolate_data import *
from numpy.random import default_rng

from prediction.cycle import const_flow
from model_fitting.cm_solver_functions import *
from data.interpolate_data import *
from scipy.optimize import curve_fit


 
def findGridBorders(Tstm, p, stepPortion = 0.1, cutoff = 10, initPars = False):
    # Inputs:
    #   p: best fit of parameters to make a grid around
    #   stepPortion: float - what decimal to increment each parameter (e.g. 0.1 = +- 10% of orginal value each time)
    #   cutoff: float - what ratio smaller than initial Prob should be considered done (higher will go longer)
    #   initPars: bool - to pass to misfit() as init_is_pars (are initial values considered as parameters)
    # Outputs:
    #   pMaxDelta: array of length p - amount you need to change each variable indivuidally to give negligible prob  

    pInc = [stepPortion*param for param in p] # set step sizes

    # setup cost / pdf for initial guess
    cost = misfit(Tstm, p, initPars)
    initialheight = np.exp(cost*-0.5)
    pMaxDelta = [0]*len(p)

    for i in range(len(p)):
        pNew = p
        ratio = 1
        while ratio < cutoff:
            pNew[i] += pInc[i]
            newCost = misfit(Tstm, pNew, initPars)

            ratio = initialheight / np.exp(newCost*-0.5) 
        pMaxDelta[i] = pNew[i] - p[i]  

    return pMaxDelta


def gridSample(p, pMaxDelta, steps, isInitPars, Tstm, numSamples):
    # Inputs:
    #   p: list - best guess at parameters
    #   pMaxDelta: list same length as p - most each parameter can be +- by
    #   steps: int (3+) - number of points to generate in each dimension - so total points to sample from will be steps**len(p)
    #   numSamples: int - how many models to generate

    bestMisfit = misfit(Tstm, p, isInitPars)
    # this factor is so that the 'best' model will always pass the random rejection process - make everything more likely to pass
    #   I assume so that it doesn't spend eons rejecting stuff
    factor = np.exp(bestMisfit*-0.5)

    pmin = p - pMaxDelta
    pmax = p + pMaxDelta

    coords = [0]*len(p)

    for i in range(len(p)):
        coords[i] = np.linspace(pmin[i], pmax[i], steps)
    
    maxIndex = [len(vector) for vector in coords]

    params = [0]*len(p)
    models = []

    while len(models) < numSamples:

        # generate a random model in the grid of options
        for i in range(len(p)): 
            # select a random value in the vector
            j = np.random.randint(0,len(coords[i]))
            # make a list of randomly chosen paramters to check
            params[i] = coords[i][j]

        # find the cost and randomly reject based on their misfit (higher misfit - more likely to be rejected)
        cost = misfit(Tstm, params, init_is_pars=isInitPars)
        r = np.random.rand()*factor
        if r < np.exp(cost*-0.5):
            models.append(params)


    return models


def normSample(p, Pcov, Tcov, samples):
    # assumes that pressure and temp parameters are unrelated
    # a better but harder to implement and more expensive form would be to use the Ppars to recompute the T everything and then sample from that
    # I don't think the parameters accross domains should affect each other so eh

    r = default_rng()

    pP = p[0:3]
    pT = p[3:6]
    Ppars = r.multivariate_normal(pP, Pcov, samples)
    Tpars = r.multivariate_normal(pT, Tcov, samples)
    s = [Ppars, Tpars]
    models = [np.append(s[0][i], s[1][i]) for i in range(samples)] 

    return models


        







