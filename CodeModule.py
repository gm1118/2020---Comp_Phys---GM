
import numpy as np
import matplotlib.pyplot as plt
import numpy.fft as fft

def floaterror(x, result='float'):
    """
    Function for finding the two nearest reals to a given float, can 
    also give the fractional rounding range where reals are rounded 
    to the given float using the kwarg 'result'
    """
    y = 0
    if bit==64:
        m_range=52
        e_range=10
    elif bit==32:
        m_range=23
        e_range=7
    if x==0:            #Consider exceptional case of x=0 as log of 0 is undefined
        m=0
        e=-(2**e_range)+2
    else:
        a = np.log2(np.abs(x))  #find log base 2 of float
        e = int(a)      #log of mantissa <1 as therefore exponent is a nearby integer
        if e>((2**e_range)+1):
            e=(2**e_range)+1
        elif e<(-(2**e_range)+2):
            e=-(2**e_range)+3        
        if (a-e) == 0:      #but must account for smaller difference for the lower float due to the smaller exponent in that area 
            y = 1           #No need to change anything if mantissa is 1 even if exponent is negative
        elif a<0:           #the mantissa increases the log, which is fine if the exponent is positive
            e -= 1          #as the int() function rounds down absolutely this needs to be compensated when negative
        m = x/(2**e)        #the mantissa
    u = m+2**-m_range        #the next possible mantissa and thus the one of the nearest upper float
    l = m-2**(-m_range-y)    #the previous possible mantissa and thus the one of the nearest lower float
    uf = u*(2**e)       #the nearest upper float
    lf = l*(2**e)       #the nearest lower float
    if x==0:            #Additional considerations for x=0 as fractions are impossible with relation to 0
        ru=uf/2
        fractru='N/A'
        rl=-lf/2
        fractrl='N/A'
        fract='N/A'
        r=ru+rl
    else:
        ru = (uf-x)/2     #The range of the upper float where real numbers are rounded to the float
        fractru = ru/x     #Converted to a fraction
        rl = (x-lf)/2      #Likewise for lower
        fractrl = rl/x
        fract = fractru+fractrl  
        r = ru+rl

    print('Its nearest upper float is',uf,' and its nearest lower float is',lf)
    if result=='fract':
        print('This means the total fractional rounding range of the orginal float is ',fract,' of the original float or ±',r/2,)
        return fractru, fractrl    #Return ranges for either side to show whether they're even or not
    elif result=='diff':
        return ru,rl     #return the range where reals are rounded to the float but not as a fraction
    else:
        return uf,lf     #return the two nearest floats, the lower first, then the upper


def CroutLU(M, Seperate=False):
    """
    Function for applying the Crout Method using the Doolittle choice to an NxN Matrix.
    The resulting LU matrix overwrites the original matrix, although
    """
    if isinstance(M,np.ndarray) == True:
        n=len(M)
        if n != len(M[0]):
            raise Exception('Function only works for NxN Matrix')
        else:
            pass                          #Arrary works better for representing matrices
    elif isinstance(M,list) == True:      #List can be converted to matrix if necessary
        n=len(M)
        if isinstance(M[0],list) != True:
            raise Exception('Function only works for NxN Matrix')
        elif n != len(M[0]):
            raise Exception('Function only works for NxN Matrix')
        else:
            m=np.zeros((n,n))
            for r in range(n):
                for c in range(n):
                    m[r][c] = M[r][c]
            M=m
    else:                               #Otherwise matrix cannot be represented
        raise Exception('Function only works for Arrays or Lists')
    L = np.zeros((n,n))
    U = np.zeros((n,n))
    for j in range(n):       #Going through each column
        L[j][j] = 1          #Dolittle Choice
        for i in range(n):   #Going through each row in the column
            if j>(i-1):
                U[i][j] = M[i][j]
                for k in range(i):            #Go through the U values in the column
                    U[i][j] -= L[i][k]*U[k][j]
            else:
                S = M[i][j]                  #Then go through the L values
                for k in range (j): 
                    S -= L[i][k]*U[k][j]     #Creating sum value so that it can then be divided by the required U value
                L[i][j] = S/(U[j][j])
    if Seperate == True:
        return L,U                           #Two variables returned when the matrices are seperate
    else:
        for j in range(n):
            L[j][j] = 0
        N = L+U
        return N                             #One variable when a single combined matrix is wanted

def DetTri(M):
    """
    Function for determining the determinant of a triangular matrix/a combined LU Matrix
    It finds the product of the diagonal terms
    """
    if isinstance(M,np.ndarray) != True:
        raise Exception('Function only works for NxN Matrix')   #Make sure it's an array
    elif len(M) != len(M[0]):
        raise Exception('Function only works for NxN Matrix')   #Make sure it's a square matrix
    else:
        Det = 1
        for i in range (len(M)):    
            Det = Det*M[i][i]       #Get the product of all the diagonal terms
        return Det
    
def LUSub(L,U,b):
    """
    Function for using forward and backwards substitution on the LU matrices of an 
    original matrix to solve the matrix equation L.U.x=b. 
    """
    if isinstance(b,np.ndarray) == True:
        if len (b[0]) != 1:
            raise Exception('b must be a 1XN vector')    #Checks to make sure b can be solved
        elif len(b) != len(L):
            raise Exception('b must be a 1xN vector')
        else:
            pass
    elif isinstance(b,list) == True:
        if isinstance(b[0],list) == True:
            raise Exception('b must be a 1XN vector')   #Can also use a list to represent the vector
        elif len(b) != len(L):
            raise Exception('b must be a 1xN vector')
        else:
            pass
    else:
        raise Exception('b must be a 1xN vector')
    if isinstance(L,np.ndarray) == True:
        if len (L) != len(U):
            raise Exception('L and U must have the same dimensions')
        elif len (L[0]) != len(U[0]):
            raise Exception('L and U must have the same dimensions')  #Checking L and U match dimensions
        elif isinstance (U,np.ndarray) == True:
            pass
        else:
            raise Exception('U must be an array')
    else:
        raise Exception('U must be an array')
    y = np.zeros((len(L),1))                 #Creating template array for y to have values inserted in
    for i in range (len(L)):
        if isinstance(b[i],int) == True:
            b[i] = (b[i]+0.0)               #Making sure no problems arrive with float/int combinatinos
        else:
            pass
        S = b[i]                 #Creating variable for sum
        for j in range (i):
            S = S-L[i][j]*y[j]   
        y[i] = S/(L[i][i])       #Determining each y value
    x = np.zeros((len(L),1))
    for i in range(len(U)):
        k = (len(U)-1-i)
        S = y[k]
        for j in range (len(U)):
            if j<(k+1):
                pass
            else:
                S = S-U[k][j]*x[j]     #repeating with x values
        x[k] = S/(U[k][k])
    return x

def LUInv(A):
    """
    Function for calculating the inverse of a matrix using its LU decomposition matrices
    """
    L,U = CroutLU(A, Seperate = True)
    l = len(A)
    I = np.zeros((l,l))
    for i in range(l):
        b = np.zeros((l,1))
        b[i] = 1.0
        c = LUSub(L,U,b)          #Applying LUSub to each column of an identity matrix to determine inverse
        for j in range (l):
            I[j][i] = c[j]
    return I
        


def InterpolatePoly(x, yin, xin):
    """
    Function for a lagrange-polynomial interpolation on xy data performed for
    a specific x-point
    """
    np.seterr(divide='ignore', invalid='ignore')
    if len(xin) != len(yin):
        raise Exception("Size of data sets must be equal")
    l = len(xin)

    nom = 1
    Intersum = 0
    Inter = np.zeros(l)      #create array for terms in the sum
    
    for i in range(l):
        nom *= (x-xin[i])      #Create variable for nominator
        Inter[i] = yin[i]

    for i in range(l):
        prod = 1               #Variable for the product of the product terms
        for j in range(l):
            if i >= j:
                pass
            else:
                p = (xin[i]-xin[j])
                prod /= p
                Inter[j] /= -p   #ensure each 'p' can be reused in a future sum where possible to reduce calculation time
        Intersum += prod*Inter[i]*nom/(x-xin[i])   #Add term to sum by multiplying product and y-value

    return Intersum


def CubicSplinef(yin, xin):
    """
    Function for a cubic spline interpolation on xy data 
    determining all of the double derivative values
    """
    if len(xin) != len(yin):
        raise Exception("Size of data sets must be equal")
    l = len(xin)
    n = l-1
    A = np.zeros((n-1,n-1))   #Create base array of zeros for the equation terms
    b = np.zeros((n-1,1))     #Create base array of zeros for equation results
    for j in range(n-1):
        i = j+1
        if j == 0:           #Special situation for first equation as must account for f''(0)=0
            fi = (yin[i+1]-yin[i])/(xin[i+1]-xin[i])
            b[j] = fi - (yin[i]-yin[j])/(xin[i]-xin[j])
            xi = (xin[i+1]-xin[i])
            A[j][j] = (xin[i+1]-xin[j])/3
            A[j][i] = xi/6
        elif j == (n-2):    #Similar situation for last equation
            b[j] = (yin[i+1]-yin[i])/(xin[i+1]-xin[i]) - fi
            xii = (xin[i+1]-xin[i])
            A[j][j-1] = xi/6
            A[j][j] = (xi+xii)/3
        else:
            fii = (yin[i+1]-yin[i])/(xin[i+1]-xin[i]) #Can reuse variables for later equations
            b[j] = (fii-fi)
            fi = fii                  #Update variable
            xii = (xin[i+1]-xin[i])
            A[j][j-1] = xi/6         #Note how because of the shortened size of this matrix (using natural spline)
            A[j][j] = (xi+xii)/3     #Its index numbers for the same valuable are different
            A[j][i] = (xii)/6        #So j here would be i in xin or yin
            xi = xii                 #Updating another variable
    L,U = CroutLU(A, Seperate=True)   #Reuse matrix solver from earlier
    f = LUSub(L,U,b)
    F = np.zeros((l,1))
    for i in range(n-1):         #Making sure values from natural spline aren't forgotten
        j = i+1                    
        F[j] = f[i]
    return F

def CubicSplinex(x,yin,xin,F):
    """
    Function for a cubic spline interpolation on xy data 
    determining individual x values
    """
    l = len(xin)
    d = np.zeros((l,1))
    for j in range(l):
        d[j] = np.abs(xin[j]-x)        #Determining nearest x values
    i = d.argmin()
    if (x-xin[i]) == 0:                #Ensuring that x is between the nearest two values
        return yin[i]
    elif (xin[i]-x)>0:
        i -= 1
    j = i+1
    X = (xin[j]-xin[i])              #Creating Cubic Spline equation
    A = (xin[j]-x)/X
    B = 1-A
    C = ((A**3-A)*(X**2))/6
    D = ((B**3-B)*(X**2))/6
    f = A*yin[i]+B*yin[j]+C*F[i]+D*F[j]   #Using y values and double derivatives calculated previously
    return f

def CubicSplineArr(xarr,yin,xin):
    """
    Function combining the previous two CubicSpline functions to provide an interpolation
    for a range of x values (the array)
    """
    F = CubicSplinef(yin,xin)         #Using first cubic spline function
    vals = np.zeros((len(xarr),1))    
    for i in range(len(xarr)):
        vals[i] = CubicSplinex(xarr[i],yin,xin,F)   #Using Second cubic spline Function
    return vals


def h(t):                       #Definition of H provided
    if 5 <= t <= 7:
        return 4
    else:
        return 0

def g(t):                     #Definition of G provided
    return ((2*np.pi)**(-0.5))*np.exp(-(t**2)/4)

def convolutehg(N,T, show=True):
    """
    Function for Convoluting two predetermined functions h(t) and g(t)
    Requires a Time period and amount of points used for the Discrete 
    Fourier Transform. Results can be shown on a plot using kwarg 'show'
    """
    dT = T/N
    tarray = np.zeros(N)
    harray = np.zeros(N)
    garray = np.zeros(N)
    for i in range(N):
        tarray[i] = i*dT
        harray[i] = h(i*dT)
        garray[i] = g(i*dT)

    gtildarray = (fft.fft(garray,n=(N)))
    htildarray = (fft.fft(harray,n=(N)))
    convtildarray = htildarray*gtildarray
    convarray = (fft.ifft(convtildarray, n=N))  
    if show == False:
        return convarray
    else:
        fix,ax2 = plt.subplots()
        ax2.plot(tarray,harray,'k-',label="h(t) with "+str(N)+" points")
        ax2.plot(tarray,garray,'g--', label='g(t) with '+str(N)+' points')
        ax2.plot(tarray,convarray,'b--',label="(g*h)(t) with "+str(N)+" points")
        ax2.set(title='Convolution of h(t) and g(t)', xlabel='t')
        ax2.legend(loc = 'upper left')
        plt.show()


def Vin1(t):        #Function to describe Vin values as they are in initial parts of 5
    if t<0:
        return V0   #Requires V0 to have been predetermined
    else:
        return 0
    
def Vin2(t):      #Function to describe Vin values as they are later in 5 with the square pulses
    if t<0:
        return V0
    for i in range(1,100):
        if Tperiod*(i-1)<=t<Tperiod*i/2:    #Tperiod needs to be predetermined
            return 0
        elif Tperiod*i/2<=t<Tperiod*i:
            return V0

            


def VoutSolver(T,h,V0,t0,R,C, Method = "AB", Vin_no = 1):
    """
    Function for solving Vout using different possible methods, default
    is the Adams Bashforth, but using kwarg 'Method', Runge-Kutta method 
    can be used here also.
    """
    def Vin(t):            #Say which Vin is being used here
        if Vin_no == 2:
            return Vin2(t)
        else:
            return Vin1(t)
    
    N = ((T-t0)/h)        #Determine number of discrete values calculated
    if (int(N)-N) == 0:
        N = int(N)
    else:
        N = int(N)
        print('The desired time step ',h,' cannot be used')
        h = ((T-t0)/N)
        print('The nearby time step ',h,' will be used instead')
        
    
    tn = np.zeros(N+1)
    Vout = np.zeros(N+1)
    for i in range(N+1):     #Create time values for each step
        tn[i] = t0+i*h
    Vout[0] = V0
    if Method == "AB":
        for i in range(0,3):    #Use Euler Method to calculate initial values
            Vout[i+1] = Vout[i]+h*(Vin(tn[i])-Vout[i])/(R*C)   #Would divide by RC but as those aren't given they can be added to units
        fd = h*(Vin(tn[0])-Vout[0])/(R*C)
        fc = h*(Vin(tn[1])-Vout[1])/(R*C)
        fb = h*(Vin(tn[2])-Vout[2])/(R*C)
        fa = h*(Vin(tn[3])-Vout[3])/(R*C)
        for i in range(3,N):   #Use 4th Order AB for rest
            Vout[i+1] = Vout[i]+(h/24)*(55*fa-59*fb+37*fc-9*fd)
            fd = fc            #Update fa,fb,fc,fd values for each n
            fc = fb
            fb = fa
            fa = h*(Vin(tn[i+1])-Vout[i+1])/(R*C)
        return Vout, tn
    elif Method != "RK":
        raise Exception("Only possible methods are Adams-Bashforth (AB) or Runge-Katta (RK)")
    else:
        for i in range(N):      #Runge Kutta Method, each step is effectively independent, no resued values
            fa = (Vin(tn[i])-Vout[i])/(R*C)
            fb = (Vin(tn[i]+h/2)-Vout[i]-h*fa/2)/(R*C)
            fc = (Vin(tn[i]+h/2)-Vout[i]-h*fb/2)/(R*C)
            fd = (Vin(tn[i]+h)-Vout[i]-h*fc)/(R*C)
            Vout[i+1] = Vout[i]+(h/6)*(fa+2*fb+2*fc+fd)  
        return Vout, tn
        
def VoutShow(T,h,V0,t0,R,C, Vin_no = 1, RK = True, AB = True, An = False):
    """
    Function for showing the results of solving the ODE for Vout on
    a graph with the possibility of comparing different methods
    """
    if AB == True:
        Vout1,tn = VoutSolver(T,h,V0,t0,R,C, Method = "AB", Vin_no = Vin_no)
    if RK == True:
        Vout2,tn = VoutSolver(T,h,V0,t0,R,C, Method = "RK", Vin_no = Vin_no)
    if An == True:
        Vanalytic = np.zeros(len(tn))
        for i in range(len(tn)):
            Vanalytic[i] = np.exp(-(t0+i*h)/(R*C))

    fix,ax3 = plt.subplots()
    if AB == True:
        ax3.plot(tn,Vout1,'k-',label="AB Method",alpha=1)
    if RK == True:
        ax3.plot(tn,Vout2,'g--', label='RK Method',alpha=1)
    if An == True:
        ax3.plot(tn,Vanalytic,'b-', label='Analytic',alpha=0.5)
    ax3.set(title='ODE Methods', xlabel='Time (s)', ylabel='Vout (V/V0)')
    ax3.legend(loc = 'upper right')
    plt.show()           
             
