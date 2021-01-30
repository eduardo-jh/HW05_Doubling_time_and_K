#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
BE523 Biosystems Analysis & Design
HW5 - Problem 1. Calculation of exponential K
https://mathinsight.org/doubling_time_half_life_discrete

Created on Tue Jan 26 22:37:49 2021
@author: eduardo
"""
from scipy import stats
import matplotlib.pyplot as plt
import numpy as np

# Bacteria density data
B = [22.1, 23.4, 26.1, 27.5, 30.5, 34.4, 36.6]  # Exercise 4C
# B = [0.015, 0.021, 0.031, 0.040, 0.055, 0.075, 0.106]  # Exercise 4B

steps = len(B)  # Adjust the lengths of the vectors with the number of steps
dB = np.zeros(steps)
Bgraph = np.zeros(steps)
Bexp = np.zeros(steps)
dt = 5
t = np.linspace(0, (steps-1)*dt, steps)  # time vector

for i in range(1, steps):
    dB[i] = B[i] - B[i-1]  # compute the increment between time steps
    Bgraph[i] = B[i-1]

print(B, Bgraph)

# Perform a linear regression with the data
slope, intercept, r_value, p_value, std_err = stats.linregress(Bgraph, dB)

# Figure 1, plotting dB vs B
plt.figure(1)
plt.plot(Bgraph, dB, 'bx', Bgraph, slope*Bgraph + intercept, 'r-')  # plot and linear eq.
plt.legend(['data', 'linear regression $R^2$=%.2f' % r_value**2], loc='best')
plt.xlabel('B')
plt.ylabel('dB')
# plt.savefig('p1_growth_%dsteps_linear.png' % steps, dpi=300, bbox_inches='tight')

print('slope =', slope, 'intercept =', intercept)

# Create the parameters for an exponential equation
tdouble = np.log(2)/np.log(1+slope)*dt  # time to double population
mu = np.log(2)/tdouble  # constant for exponential eq. with base on natural log
print('tdouble =', tdouble, 'mu=', mu)

# Generate an exponential equation (exact solution)
Bexp = B[0] * np.exp(mu*t)

# Generate predictions with the model: B(t+1) = B[0]*(r+1)^t (numerical solution)
Bmodel = B[0]*(1+slope)**(t/dt)

# Figure 2, plotting B vs t
plt.figure(2)
plt.plot(t, B, 'bx', t, Bmodel, 'r-', t, Bexp, 'k+')  # Plot data vs exponential growth eq.
plt.legend(['data',
            'equation B=%.2f*%.3f^t' % (B[0], slope+1),
            'exponential B=%.2f * exp(%.3f*t)' % (B[0], mu)],
           loc='best')
plt.xlabel('time (min)')
plt.ylabel('B')
plt.savefig('p1_growth_%dsteps.png' % steps, dpi=300, bbox_inches='tight')
