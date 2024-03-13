#############################################################################
# Create a graph similar to https://www.nature.com/articles/s41592-018-0293-7
# Where we choose a length l and calculate the probability that we find a mutation
# at distance i+l given that there was a mutation at i
#############################################################################

from lmfit import Parameters, Minimizer, minimize
import matplotlib.pyplot as plt
import math
import numpy as np
import pandas as pd
import random


def Power(a, b):
    """compute power"""
    return a**b


def const_r1(x, fBar, phiC, w):
    """calculate r1 assuming constant fragment size"""
    return np.where(x < fBar, w*phiC*x, w*phiC*fBar)


def calcP2(thetaS, r1, r2, ds, a):
    """
    calcP2 using expression computed using Mathematica CForm
    """
    v = (2*(r2*thetaS + ds*r1*(1 + r1 + r2 + a*thetaS))* \
        (r2*Power(thetaS,2) + Power(ds,2)*(1 + r1 + r2 + a*thetaS)* \
        (2*Power(r1,2) + r2 + 3*r1*r2 + Power(r2,2) + a*(r1 + 2*r2)*thetaS) - \
        ds*thetaS*(2*r2 + Power(r1 + r2,2) + a*(r1 + 3*r2)*thetaS)))/ \
        (Power(r1 + r2,2)*(1 + 2*r1 + r2 + 2*a*thetaS)* \
        (-(thetaS*(r1 - r2 + a*thetaS)) + ds*(2*r1 + a*thetaS)* \
        (1 + r1 + r2 + a*thetaS)))
    return v


def fcn2min(params, xvalues, yvalues, r1_func):
    """function 2 min"""
    thetaS = params['thetaS']
    phiS = params['phiS']
    f = params['f']
    w = params['w']
    r1 = r1_func(xvalues, f, phiS, w)
    r2 = phiS * w * f - r1
    ds = params['ds']
    a = params['a']
    p2 = calcP2(thetaS, r1, r2, ds, a) / ds
    return p2 - yvalues


def fit_model(xvalues, yvalues, d_sample, r1_func):
    """fitting correlation profile using lmfit"""
    params1 = Parameters()
    params1.add('ds', value=d_sample, vary=False)
    params1.add('thetaS', value=0.00001, min=0, max=d_sample)
    params1.add('f', value=1000, min=3, max=300000)
    params1.add('phiS', value=0.00005, min=0, max=1)
    params1.add('w', value=2.0/3.0, vary=False)
    params1.add('a', value=4.0/3.0, vary=False)
    params1.add('thetaP', expr='(ds*(1 + phiS*w*f + a*thetaS)-thetaS)/ \
                                ((1 - a*ds)*(phiS*w*f + a*thetaS)-(a*ds))')
    params1.add('phiP', expr='phiS*thetaP/thetaS')
    params1.add('c', expr='w*phiS*f/(1+w*phiS*f+thetaS*a)')
    params1.add('dp', expr='thetaP/(1+a*thetaP)')
    params1.add('dc', expr='thetaS/(1+a*thetaS)')
    result = minimize(fcn2min, params1, args=(xvalues, yvalues, r1_func),
                      method="least_squares", max_nfev=int(1e6))
    return result


# Read the data.
data_file = '/home/ada/Desktop/PinkBerry_scripts_paper/psb_scripts/recombination_rate/mcorr_PSB_data_2024.csv'

xvalues, yvalues = [], []

with open(data_file, 'r') as f:
    for line in f:
        x, y = line.split(',')
        xvalues.append(int(x.strip()))
        yvalues.append(float(y.strip()))

xvalues = np.array(xvalues)
yvalues = np.array(yvalues)

pl_results = []
l_list = []
for l in range(1,600):
    ds = 0.5011   # calculated
    thetaS = 1e-02
    f = 300
    phiS = 7e-03
    w = 2.0/3.0
    a = 4.0/3.0
    thetaP = (ds*(1 + phiS*w*f + a*thetaS)-thetaS)/ ((1 - a*ds)*(phiS*w*f + a*thetaS)-(a*ds))
    phiP = phiS*thetaP/thetaS
    c = w*phiS*f/(1+w*phiS*f+thetaS*a)
    dp = thetaP/(1+a*thetaP)
    dc = thetaS/(1+a*thetaS)
    c_0 = (1+2*thetaS*a) / (1 + 2*thetaS*a + phiS*w*(f+l))
    c_1 = (2*phiS*w*l) / (1+2*thetaS*a + phiS*w*(f+l))
    c_2 = (phiS*w*(f-l)) / (1 + 2*thetaS*a + phiS*w*(f+l))
    d_2theta = (2*thetaS) / (1 + 2*thetaS*a)
    qp_part1 = 2 * math.pow((thetaP / (1+thetaP*a)), 2)
    qp_part2 = (1 + thetaP*a + phiP*w*l) / (1 + 2*thetaP*a + 2*phiP*w*l)
    qp = qp_part1 * qp_part2 
    pl = (c_0 * d_2theta*ds + c_1*ds*dp + c_2*qp) / ds

    pl_results.append(pl)
    l_list.append(l)

plt.scatter(xvalues, yvalues, alpha=0.5)
plt.plot(l_list, pl_results)
plt.show()
