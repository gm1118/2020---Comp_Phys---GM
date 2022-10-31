
import numpy as np
import matplotlib.pyplot as plt
import ProjectCodeModule as P
import random
import time


#%%
"""
Question 2c Results:

First decide on values for the constants I decided to create a constant 
value called 'C' which is equal to m*omega/hbar. It helps simplify the 
integral and the reasons for its use are explained in more detail in the 
report. Ultimately it should not affect the result of the integral 
(although it will vary within the accepted error, even for the trapezium data),
so I have set it to 2 to help the figures have nicer axes and numbers, 
If you wish to test this, set C to any valid value it can take (i.e. C>0) and
the 1D integral will stay the same
"""
C = 2


#The user-defined 1D integrand used for the following results, note this is the pdf, not the wavefunction!
def Psi0(x):
    return ((C/np.pi)**0.5)*np.exp(-C*x**2)


def psi1(x):
    return (2*C*x**2)*((C/np.pi)**0.5)*np.exp(-C*x**2)

def a(x):
    return -4*x/7 + (1/1.4) + 2*1.4/7


#%%
# Code for creating Figure 1
C = 2
xt = np.arange(0,2/(C**0.5),0.0001)  #High enough density of points for smooth line

fix,ax1 = plt.subplots()
ax1.plot(xt, Psi0(xt), 'k-', label = '|Psi0(x)|^2')
ax1.set(xlabel='x', ylabel='Probability', title = 'Pdf of Psi0(x) within Integrating Period', xlim = (0,1.5))
ax1.grid()
ax1.legend(loc='upper right', ncol=2)
plt.show

#%%
# Code for creating Figure 2
fix,ax1 = plt.subplots()
ax1.plot(xt, psi1(xt), 'kx', label = 'Original Data')

ax1.set(xlabel='X axis', ylabel='Y axis', title = 'Lagrange Polynomial & Cubic Spline Interpolation')
ax1.grid()
ax1.legend(loc='upper right', ncol=2)
plt.show



#%%
# Code for Question 2 Results Trapezoidal
start = time.time()
I, N = P.ExtTrap(Psi0, 0, 2/(C**0.5), e = 10**-6, Open = 'None', n=all, Nreturn = True)
end = time.time()
runtime = end - start
n = 2**N
print('Estimate of',I,'Calculated using ',n,'function evaluations, taking',runtime,'seconds')

#%%
# Code for Question 2 Results Simpson's
start = time.time()
I, N = P.ExtSimp(Psi0, 0, 2/(C**0.5), e = 10**-6, Open = 'None', n=all, Nreturn = True)
end = time.time()
runtime = end - start
n = 2**(N+1)
print('Estimate of',I,'Calculated using ',n,'function evaluations, taking',runtime,'seconds')


#%%
#3a results
Shots = 100000
start = time.time()
I, err = P.MonteCarloIntegr(Psi0, 0, 2/(C**0.5), Shots, P.SampleRandom, GiveErr = True)
end = time.time()
runtime = end - start
print('Estimate of ',I,'Calculated using',Shots,'random sample points with an error of',err,', taking',runtime,'seconds')
#%%
A = 4/7
def SamplePdf(x0,x1,Shots, pdf = False):
    """
    Function for obtaining a non-uniform sample between x0 and x1
    with a weighting determined by the Transformation method being
    applied using a pdf of the form 'y = Ax+b' where A is within a 
    specific range to ensure normalisation and b can be defined from
    A.
    
    A has to be chosen, so the function definition has been placed here 
    rather than in P so that a global variable A can be used to adjust A.
    """
    b = ((x1-x0)**-1)-A*(x0+x1)/2
    xlist = []
    #Create data points
    for i in range(Shots):
        y = random.random()
        x = -b/A + (((b**2)+2*A*((x0**2)*A/2+b*x0+y))**0.5)/A      #transform y samples to x data points in the x0,x1 range
        xlist.append(x)
    if pdf == True:
        def Pdf(x):
            return (A*x+b)
        return xlist, Pdf
    else:
        return xlist
#%%
#3b results
Shots = 10000000
start = time.time()
I, err = P.MonteCarloIntegr(Psi0, 0, 2/(C**0.5), Shots, SamplePdf, GiveErr = True)
end = time.time()
runtime = end - start
print('Estimate of ',I,'Calculated using',Shots,'random sample points with an error of',err,', taking',runtime,'seconds')
#%%
#4a results
Box_Edges = [0,2/(C**0.5),0,2/(C**0.5),0,2/(C**0.5)]
Uniform_Samples = [P.SampleRandom,P.SampleRandom,P.SampleRandom]
def Psi0_3D(x,y,z):
    fx = Psi0(x)
    fy = Psi0(y)
    fz = Psi0(z)
    return fx*fy*fz

Uniform_Samples2 = [SamplePdf,SamplePdf,SamplePdf]

def Psi1_3D(x,y,z):
    return ((C**5/np.pi**3)**0.5)*(x**2+y**2)*np.exp(-C*(x**2+y**2+z**2))



#All of the Results used in 4:
#%%
start = time.time()
I, err = P.ThreeDMonteCarlo(Psi0_3D, Box_Edges, 10000000, Uniform_Samples, GiveErr = True)
end = time.time()
runtime = end - start
print('Estimate of ',I,'Calculated using',Shots,'random sample 3D points with an error of',err,', taking',runtime,'seconds')
#%%
start = time.time()
I, N = P.ThreeDNewtonCoates(Psi0_3D, Box_Edges, e=0.001, n=all, Nreturn = True)
end = time.time()
runtime = end - start
n = 2**(N*3)
print('Estimate of ',I,'Calculated using',n,'function evaluations taking',runtime,'seconds')

#%%
start = time.time()
I, err = P.ThreeDMonteCarlo(Psi0_3D, Box_Edges, 10000000, Uniform_Samples2, GiveErr = True)
end = time.time()
runtime = end - start
print('Estimate of ',I,'Calculated using',Shots,'random sample 3D points with an error of',err,', taking',runtime,'seconds')
#%%
start = time.time()
I, err = P.ThreeDMonteCarlo(Psi1_3D, Box_Edges, 100000, Uniform_Samples, GiveErr = True)
end = time.time()
runtime = end - start
print('Estimate of ',I,'Calculated using',Shots,'random sample 3D points with an error of',err,', taking',runtime,'seconds')
#%%
start = time.time()
I, N = P.ThreeDNewtonCoates(Psi1_3D, Box_Edges, e=0.001, n=all, Nreturn = True)
end = time.time()
runtime = end - start
n = 2**(3*N)
print('Estimate of ',I,'Calculated using',n,'function evaluations taking',runtime,'seconds')



