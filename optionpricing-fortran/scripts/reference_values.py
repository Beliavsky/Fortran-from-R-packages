#!/usr/bin/env python3
"""Print independent Black-Scholes and Curran ECV reference values."""
import math
from statistics import NormalDist

N = NormalDist()

def pdf(x: float) -> float:
    return math.exp(-0.5*x*x)/math.sqrt(2.0*math.pi)

def bs_call(t: float, k: float, r: float, sigma: float, s0: float):
    d1 = (math.log(s0/k)+(r+0.5*sigma*sigma)*t)/(sigma*math.sqrt(t))
    d2 = d1-sigma*math.sqrt(t)
    return (s0*N.cdf(d1)-k*math.exp(-r*t)*N.cdf(d2),
            N.cdf(d1), pdf(d1)/(s0*sigma*math.sqrt(t)))

def eval_ecv(t=1.0, d=12, kstrike=100.0, r=0.05, sigma=0.2, s0=100.0):
    dt=t/d
    varx=d*(d+1)*(2*d+1)/6
    mus=math.log(s0)+(r-0.5*sigma*sigma)*dt*(d+1)/2
    sigmas=sigma/d*math.sqrt(dt*varx)
    zcut=(math.log(kstrike)-mus)/sigmas
    avec=[]
    running=0.0
    for j in range(d):
        running+=(d-j)/math.sqrt(varx)
        avec.append(sigma*math.sqrt(dt)*running)
    sum1=sum(math.exp(r*(i+1)*dt)*N.cdf(-zcut+avec[i]) for i in range(d))
    hk=s0/d*sum(math.exp(a*zcut+r*(i+1)*dt-0.5*a*a)
                for i,a in enumerate(avec))
    hdk=s0/d*sum(a*math.exp(a*zcut+r*(i+1)*dt-0.5*a*a)
                 for i,a in enumerate(avec))
    price=s0/d*sum1-kstrike*N.cdf(-zcut)
    delta=sum1/d+(hk-kstrike)*pdf(zcut)/(s0*sigmas)
    gamma=(2*hk*sigmas-hdk-(hk-kstrike)*(sigmas-zcut))*pdf(zcut)/(s0*sigmas)**2
    disc=math.exp(-r*t)
    return disc*price,disc*delta,disc*gamma

print("BS call:", *bs_call(0.25,100.0,0.05,0.2,100.0))
print("Curran ECV:", *eval_ecv())
