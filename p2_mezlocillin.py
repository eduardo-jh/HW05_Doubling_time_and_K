#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
BE523 Biosystems Analysis & Design
HW5 - Problem 2. Exercise 3 Mezlocillin model
https://mathinsight.org/penicillin_clearance_model_exercises

Created on Sat Jan 23 21:29:34 2021
@author: eduardo
"""
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm  # allows linear regression without intercept

# The mezlocillin concentration that will be used to find the model
M1 = np.array([71, 56, 45, 33, 25])
M5 = np.array([490, 390, 295, 232, 182])
assert len(M1) == len(M5), "The arrays should have same length"
steps = len(M5)
dt = 5  # time interval, min
t = np.linspace(0, (steps-1)*dt, steps)  # actual time vector
dM1 = np.zeros(steps)
dM5 = np.zeros(steps)

# Get the decrement of mezlocillin in every time interval
for i in range(0, steps-1):
    dM1[i] = M1[i+1] - M1[i]
    dM5[i] = M5[i+1] - M5[i]

print('M1=', M1, '\ndM1=', dM1)
print('M5=', M5, '\ndM5=', dM5)

# Plot mezlocillin concentration change (dM5) versus mezlocillin concentration (M5)
# Perform a linear regression with M5 and dM5, the slope of the line will give us
# the constant rate of mezlocillin clearance in every time step
# If the model is correct, the points should lie on a line that goes through
# the origin, no intercept is desired
# For 1 g
model1g = sm.OLS(dM1, M1)  # No intercept by default
results1g = model1g.fit()
slope_1g = results1g.params[0]
print("slope 1g=", slope_1g)

# For 5 g
model5g = sm.OLS(dM5, M5)  # No intercept by default
results5g = model5g.fit()
slope_5g = results5g.params[0]
print("slope 5g=", slope_5g)

# Figure 1, plotting dM5 vs M5
plt.figure(1)
plt.plot(M1, dM1, 'kx', M1, slope_1g*M1, 'r--', M5, dM5, 'k+', M5, slope_5g*M5, 'b--')
plt.legend(['data 1g injection',
            'linear reg. 1g $R^2$=%.2f' % results1g.rsquared,
            'data 5g injection',
            'linear reg. 5g $R^2$=%.2f' % results5g.rsquared],
           loc='best')
plt.xlabel('Mezlocillin concentration')
plt.ylabel('Change dM')
plt.savefig('p2_mezlocillin_linear.png', dpi=300, bbox_inches='tight')

# Generate an exponential equation (exact solution) for the 1g injection
tdouble_1g = (np.log(0.5)/np.log(1+slope_1g))*dt
mu_1g = np.log(0.5)/tdouble_1g
M1exp = M1[0] * np.exp(mu_1g*t)

# ... and for the 5g injection
tdouble_5g = (np.log(0.5)/np.log(1+slope_5g))*dt
mu_5g = np.log(0.5)/tdouble_5g
M5exp = M5[0] * np.exp(mu_5g*t)
print('tdouble_1g =', tdouble_1g, 'mu_1g =', mu_1g, 'b_1g=', slope_1g+1)
print('tdouble_5g =', tdouble_5g, 'mu_5g =', mu_5g, 'b_5g=', slope_5g+1)

# Make predictions with the model M(t+1) = M[0]*(slope+1)^t
M1model = M1[0]*(slope_1g+1)**(t/dt)
M5model = M5[0]*(slope_5g+1)**(t/dt)

# Figure 2, plotting M5 (from data and model) vs t
plt.figure(2)
plt.plot(t, M1, 'kx', t, M1model, 'r-', t, M1exp, 'g^', t, M5, 'k+', t, M5model, 'b-', t, M5exp, 'c^')
plt.legend(['data 1g injection',
            'eq $M_{1}$=%d$\cdot$%.3f$^t$' % (M1[0], slope_1g+1),
            'exponential 1g',
            'data 5g injection',
            'eq $M_{5}$=%d$\cdot$%.3f$^t$' % (M5[0], slope_5g+1),
            'exponential 5g'],
           loc='best')
plt.xlabel('Time (min)')
plt.ylabel('Mezlocillin conc. ($\mu g/ml$)')
plt.savefig('p2_mezlocillin_model.png', dpi=300, bbox_inches='tight')