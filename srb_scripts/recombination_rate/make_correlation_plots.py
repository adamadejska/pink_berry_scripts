#############################################################################
# Create a graph similar to https://www.nature.com/articles/s41592-018-0293-7
# Where we choose a length l and calculate the probability that we find a mutation
# at distance i+l given that there was a mutation at i
#############################################################################

from lmfit import Parameters, Minimizer, minimize
import matplotlib.pyplot as plt
from matplotlib import gridspec
import math
import numpy as np
import pandas as pd


def Power(a, b):
    """compute power"""
    return a**b


def const_r1(x, fBar, phiC, w):
    """calculate r1 assuming constant fragment size"""
    return np.where(x < fBar, w*phiC*x, w*phiC*fBar)


def exp_r1(x, fBar, phiC, w):
    """calculate r1 assuming exponetional decay of fragment size"""
    return w*phiC*fBar*(1.0 - np.exp(-x/fBar))


def geom_r1(x, fBar, phiC, w):
    """calculate r1 assuming geom distribution"""
    prob = 1.0/fBar
    return w*phiC*fBar*(1.0 - np.power(1-prob, x))


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
    params1.add('ds', value=d_sample, vary=False)   # vary=False will prevent the value from changing in the fit
    params1.add('thetaS', value=0.0009, min=0.0009, max=d_sample)
    params1.add('f', value=1000, min=3, max=4000)
    params1.add('phiS', value=0.00009, min=0.00009, max=1)
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


def plot_fit(xvalues, yvalues, fitres, title):
    """Fit all row data and do plotting for the full-recombination model"""
    fig = plt.figure(tight_layout=False)

    figsize = 4
    fig.set_figheight(figsize)
    fig.set_figwidth(figsize)
    gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1], hspace=0)
    ax1 = plt.subplot(gs[0, 0])
    ax1.scatter(xvalues, yvalues, s=20, facecolors='none', edgecolors='k')
    predictions = yvalues + fitres.residual
    ax1.plot(xvalues, predictions, 'r')
    ax1.set_ylabel(r'$P$')
    if np.min(yvalues) != np.max(yvalues):
        ax1.set_ylim([np.min(yvalues)*0.9, np.max(yvalues)*1.1])
    ax1.locator_params(axis='x', nbins=5)
    ax1.locator_params(axis='y', nbins=5)
    plt.setp(ax1.get_xticklabels(), visible=True)
    ax1.set_xlabel(r'$l$')
    plt.title(title)

    plt.show()


# Read the data.
#data_file = '/home/ada/Desktop/PinkBerry_scripts_paper/psb_scripts/recombination_rate/data/mcorr_SRB_data_2024_ql_ds_cov_3.csv'
data_file = '/home/ada/Desktop/PinkBerry_scripts_paper/psb_scripts/recombination_rate/data/srb_clades/mcorr_SRB_data_2024_ql_ds_cov_3_F_clade.csv'
xvalues, yvalues = [], []

with open(data_file, 'r') as f:
    for line in f:
        x, y = line.split(',')
        xvalues.append(int(x.strip()))
        yvalues.append(float(y.strip()))

xvalues = np.array(xvalues)
yvalues = np.array(yvalues)

d_sample = 0.25960
title = 'P(l) vs distance for SRB (coverage >= 3)\n mixing layer'
r1_func = const_r1
result = fit_model(xvalues, yvalues, d_sample, r1_func)
params = result.params.valuesdict()
thetaS = result.params["thetaS"]
phiS = result.params["phiS"]
f = result.params["f"]


print("thetaS (init)", thetaS.init_value)
print("f (init)", f.init_value)
print("phiS (init)", phiS.init_value)
print("ds", "thetaS", "f", "phiS", "thetaP", "phiP", "c", "dp", "dc")
print(params["ds"], params["thetaS"], params["f"], params["phiS"],params["thetaP"], params["phiP"], params["c"], params["dp"], params["dc"])

plot_fit(xvalues, yvalues, result, title)

