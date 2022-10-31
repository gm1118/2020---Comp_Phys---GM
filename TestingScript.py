#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 16 10:47:32 2020

@author: georgemercieca
"""

import numpy as np
import matplotlib.pyplot as plt
import CodeModule as c
#%%
#Results for Q1
# a)
x=0.25
a,b = c.floaterror(x)  #a is the nearest upper float abd b the nearest lower float to 0.25
print(a,b)
#%%
#b)
a1,a2 = c.floaterror(a, result='float') #a1 and a2 are the two fractional rounding ranges near to a, a1 the higher and a2 the lower
b1,b2 = c.floaterror(b, result='float') #b1 and b2 are the likewise for b. The function also prints out the floats nearest to them.

#%%
#Results for Q2
#b)
A=np.array([[3,1,0,0,0],[3,9,4,0,0],[0,8,20,10,0],[0,0,-22,31,-25],[0,0,0,-35,61]])
L,U = c.CroutLU(A, Seperate=True)
d=c.DetTri(U)    #Due to Dolittle choice, only Upper matrix counts towards determinant as lower matrix has determinant of 1
print(d)        
#%%
#d)
b = [2,5,-4,8,9]
x = c.LUSub(L,U,b)  #Taking previously calculated U and L matrices
print(x)
#%%
#e
Inv = c.LUInv(A)   #Complete function for determining matrix inverse
print(Inv)
#%%
#Results for Q3
#c)
ydata = np.array([0.10,0.30,0.47,0.66,0.60,0.54,0.30,0.15,-0.32,-0.54,-0.60,-0.47,-0.08])
xdata = np.array([-0.75,-0.5,-0.35,-0.1,0.05,0.1,0.23,0.29,0.48,0.6,0.92,1.05,1.5])
xt = np.arange(-0.75,1.5,0.0001)  #High enough density of points for smooth interpolation

fix,ax1 = plt.subplots()
ax1.plot(xdata, ydata, 'kx', label = 'Original Data')
ax1.plot(xt, c.InterpolatePoly(xt, ydata, xdata), 'g', label = 'Lagrange')
ax1.plot(xt, c.CubicSplineArr(xt,ydata,xdata),'b-', label='Spline')

ax1.set(xlabel='X axis', ylabel='Y axis', title = 'Lagrange Polynomial & Cubic Spline Interpolation', ylim=(-1,1))
ax1.grid()
ax1.legend(loc='upper right', ncol=2)
plt.show
#%%
#Results for Q4
N=500
T=10
c.convolutehg(500,10)
#%%
#Results for Q5
#c)
#Setting base variables
T = 10              #Time Period 
V0 = 1              #Initial Vout, set to 1 so it is incorporated in results/axis
t0 = 0              #Initial time value, normally should be 0
ht = 0.01           #Time step
R = 2               #Resistance
C = 3               #Capacitance

c.VoutShow(T,ht,V0,t0,R,C, Vin_no = 1, RK = True, AB = True, An = True)
#%%
#d) 1st Graph
ht = 0.005
c.VoutShow(T,ht,V0,t0,R,C, Vin_no = 1, RK = True, AB = True, An = True)
#%%
#d) 2nd Graph
ht = 0.02
c.VoutShow(T,ht,V0,t0,R,C, Vin_no = 1, RK = True, AB = True, An = True)
#%%
#e) 1st Graph
ht = 0.01
Tperiod = 2*R*C
c.VoutShow(T,ht,V0,t0,R,C, Vin_no = 1, RK = True, AB = False, An = False)
#%%
#e) 2nd Graph
Tperiod = 2*R*C
c.VoutShow(T,ht,V0,t0,R,C, Vin_no = 1, RK = True, AB = False, An = False)
