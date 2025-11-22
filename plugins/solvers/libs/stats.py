""" 
Custom statistics functions for GeoDT
"""
#initializations
import numpy as np

def exponential_trunc(nsam,bval=1.0,Mmax=5.0,Mwin=1.0,prob=0.1): # random samples from a truncated exponential distribution
    # zero-centered Mcap
    lamba = np.log(10)*bval
    Mcap = (-1.0/lamba)*np.log(prob/(np.exp(lamba*Mwin)-1.0+prob))
    # calculated Mmin
    Mmin = Mmax-Mcap
    # sample sets with resamples out-of-range values
    s0 = np.random.exponential(1.0/(np.log(10)*bval), nsam) + Mmin
    iters = 0
    while 1:
        iters += 1
        r_pl = s0 > Mmax
        if np.sum(r_pl) > 0:
            r_dr = np.random.exponential(1.0/(np.log(10)*bval), nsam) + Mmin
            s0 = s0*(1-r_pl) + r_dr*(r_pl)
        else:
            break
        if iters > 100:
            break
    return s0

def contact_trunc(nsam,weight=0.15,
                  bd_nom=0.002,stddev=0.5*0.002,
                  exp_B=0.5,exp_C=0.5/np.pi,
                  bd_min=0.0001,bd_max=3.0*0.002): # get random samples from a 'contact' distribution
    # binomial samples
    n1 = np.random.binomial(nsam,weight,(1))
    n2 = nsam - n1
    # exponential samples
    s1 = np.random.exponential(exp_B, n1)*exp_C*bd_nom + bd_min
    iters = 0
    while 1:
        iters += 1
        r_pl = s1 > bd_max
        if np.sum(r_pl) > 0:
            r_dr = np.random.exponential(exp_B, n1)*exp_C*bd_nom + bd_min
            s1 = s1*(1-r_pl) + r_dr*(r_pl)
        else:
            break
        if iters > 100:
            break
    # normal samples
    s2 = np.random.normal(0,1,n2)*stddev + bd_nom
    iters = 0
    while 1:
        iters += 1
        r_pl = (s2 > bd_max) + (s2 < bd_min)
        if np.sum(r_pl) > 0:
            r_dr = np.random.normal(0,1,n2)*stddev + bd_nom
            s2 = s2*(1-r_pl) + r_dr*(r_pl)
        else:
            break
        if iters > 100:
            break
    s0 = np.concatenate((s1,s2),axis=0)
    return s0

def lognorm_trunc(nsam,logmu=0.0,logdev=1.0,
                  loglo=-2.0,loghi=2.0): # get random samples from a log normal distribution
    # normal samples
    s0 = np.random.normal(0,1,nsam)*logdev + logmu
    iters = 0
    while 1:
        iters += 1
        r_pl = (s0 > loghi) + (s0 < loglo)
        if np.sum(r_pl) > 0:
            r_dr = np.random.normal(0,1,nsam)*logdev + logmu
            s0 = s0*(1-r_pl) + r_dr*(r_pl)
        else:
            break
        if iters > 100:
            break
    return 10.0**s0

def norm_trunc(nsam,mu=0.0,dev=1.0,
                  lo=-2.0,hi=2.0): # get random samples from a log normal distribution
    s0 = np.random.normal(0,1,nsam)*dev + mu
    iters = 0
    while 1:
        iters += 1
        r_pl = (s0 > hi) + (s0 < lo)
        if np.sum(r_pl) > 0:
            r_dr = np.random.normal(0,1,nsam)*dev + mu
            s0 = s0*(1-r_pl) + r_dr*(r_pl)
        else:
            break
        if iters > 100:
            r_dr = np.random.uniform(0,1,nsam)*(hi-lo) + lo
            s0 = s0*(1-r_pl) + r_dr*(r_pl)
            break
    return s0

def uniform(nsam,lo,hi):
    return np.random.uniform(0,1,nsam)*(float(hi)-float(lo))+float(lo)

def loguniform(nsam,lo,hi):
    return 10.0**(np.random.uniform(0,1,nsam)*(np.log10(float(hi))-np.log10(float(lo)))+np.log10(float(lo)))