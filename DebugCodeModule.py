import numpy as np
import math
import random

def ExtTrap(f, x0, x1, e = 0.01, Open = 'None', n=all, Nreturn = False):
    """
    Function for implementing the extended Trapeziodal rule.
    Takes the integrand function, start point, end point and 
    relative accuracy as arguments. Also has kwarg for deciding if
    it is open at its ends or closed and kwarg for set number of 
    iterations instead of accuracy requirement.
    """
    # Checks to make sure input variables are compatible
    #if type(f) != 'function':
     #   raise Exception("Integrand function must be a function")
    if isinstance(x0,int) == True or isinstance(x0,float) == True:
        pass
    else:
        raise Exception('Start point must be a number')
    if isinstance(x1,int) == True or isinstance(x1,float) == True:
        pass
    else:
        raise Exception('End point must be a number')
    if x0 == x1:
        raise Exception('Start and End point must be different numbers')
    if isinstance(e,int) == True or isinstance(e,float) == True:
        pass
    else:
        raise Exception('Relative Accuracy must be a number')
    if isinstance(n,int) == True and n>1:
        pass
    elif n==all:
        pass 
    else:
        raise Exception("Number of Iterations must be an integer greater than 1, set to 'all' if you wish to calculate for a set accuracy, not a set number of iterations")
    N = 1                  # number of function evaluations    
    h = (x1-x0)         # Step size
    
    if Open == 'None':
        # First Solving Closed Trapezium rule
        # First Evaluate integrand at the given end points
        
        f0 = f(x0)
        f1 = f(x1)
        #checks to make sure divergent function values aren't present
        if math.isnan(f0) == True:
            raise Exception("Function value at ",x0," is undefined, set kwarg Open to 'Start' or 'Both', then run function again")
        elif math.isnan(f1) == True:
            raise Exception("Function value at ",x1," is undefined, set kwarg Open to 'End' or 'Both', then run function again")
        I1 = h*0.5*(f0+f1)
        
        # Evaluate Second trapzium rule iteration
        N += 1
        h /= 2
        fi = f((x0+x1)/2)
        # Divergent Value check, this check is applied to all values of f() used, with suggested alternatives provided 
        if math.isnan(fi) == True:
            raise Exception('Function value at ',(x0+x1)/2,' is undefined, use ExtTrapAsymp() to perform integral over your desired period')
        I2 = (h*f1 + I1)/2
        diff = np.abs((I2-I1)/I1)
        # Optional variable for bypassing accuracy check and only completing set number of iterations instead
        if n == N:
            if Nreturn == True: #Optional Kwarg for returning N also
                return I2, N
            else:
                return I2
        # Create loop adding additional iterations until accuracy check is satisfied
        while diff > e:
            I1 = I2
            I2 = 0
            h /= 2
            # Create all the new values for this iteration
            for i in range(2**(N-1)):
                xi = ((x1-x0)*(1+2*i)/(2**N))+x0
                fi = f(xi)
                if math.isnan(fi) == True:
                    raise Exception('Function value at ',xi,' is undefined, use ExtTrapAsymp() to perform integral over your desired period')
                I2 += fi
            I2 *= h
            I2 += I1/2                  #New integral estimation
            diff = np.abs((I2-I1)/I1)   #New accuracy value
            N += 1
            # If iteration check is reached before accuracy check, result is returned
            if n == N:
                return I2
        # Once accuracy check is completed result is provided along with number of iterations if requested
        if Nreturn == True:
            return I2, N
        else:
            return I2
    elif Open != 'Start' or Open != 'End' or Open != 'Both':
        raise Exception("kwarg 'Open' can only have values 'None', 'Start', 'End', or 'Both' ",Open," is not a valid value")
    else:
        #Solve Trapezium rule with open end(s)
        #Midpoint Rule iteration to begin
        f0 = f((x0+x1)/2)
        if math.isnan(f0) == True:
            raise Exception('Function value at ',(x0+x1)/2,' is undefined, use ExtTrapAsymp() to perform integral over your desired period')
        I1 = h*f0
        N += 1
        h /= 2
        # Different calculations for different values of 'Open', this couldn't be done earlier as N=1 means midpoint rule cannot be replaced at any point
        if Open == 'End':
            #If only the end needs to be open, trapezium rule can be used at start
            f1 = f(x0)
            if math.isnan(f1) == True:
                raise Exception("Function value at ",x0," is undefined, set kwarg Open to 'End' or 'Both', then run function again")
        else:
            #Otherwise use of midpoint rule is necessary again here
            f1 = f((x1-x0)/4+x0)
            if math.isnan(f1) == True:
                raise Exception('Function value at ',(x1-x0)/4+x0,' is undefined, use ExtTrapAsymp() to perform integral over your desired period')
        # and vice versa at the end of the period
        if Open == 'Start':
            f2 = f(x1)
            if math.isnan(f2) == True:
                raise Exception("Function value at ",x1," is undefined, set kwarg 'Open' to 'End' or 'Both', then run function again")
        else:
            f2 = f((x1-x0)*3/4+x0)
            if math.isnan(f1) == True:
                raise Exception('Function value at ',(x1-x0)*3/4+x0,' is undefined, use ExtTrapAsymp() to perform integral over your desired period')
        if Open == 'Both':
            I2 = h*(f1+f2)
        elif Open == 'Start':
            I2 = h*(f1+(f0/2)+(f2/4))
            f0 += (f2/2)
            f2 = 0
        elif Open == 'End':
            I2 = h*((f1/2)+(f0/2)+f2)
            f0 += (f1/2)
            f1 = 0
        diff = np.abs((I2-I1)/I1)
        f0 /= 2  #f0 represents a core of trapezoid iterations that can be reused from N=3 onwards and will be updated accordingly
        F1 = 0   #F1 & F2 will represent the function values used at the 'edges' of the trapezoid iterations
        F2 = 0   #They aren't used for N=1 or N=2 so they are set to 0 here
        #Iteration Check
        if n == N:
            if Nreturn == True:
                return I2, N
            else:
                return I2
           #Accuracy Checks
        while diff > e:
            I1 = I2
            f0 /= 2          #Adjusting 'core' for new trapezoidal iteration
            f0 += (F1+F2)    #Adding previous values to 'core' so they don't need to be recalculated
            F1 = f1/4        #Creating new 'edge' trapeziodal values from previous midpoint calculations
            F2 = f2/4
            h /= 2
            
            #Calculating new trapezoidal values for this iteration
            for i in range(2**(N-1)):
                if i == 0 or i == (2**(N-1)-1):
                    pass
                xi = (x1-x0)*(1+2*i)/(2**N)+x0
                fi = f(xi)
                if math.isnan(fi) == True:
                    raise Exception('Function value at ',xi,' is undefined, use ExtTrapAsymp() to perform integral over your desired period')
                fi *= h
                f0 += fi     #adding new trapezoidal values to the 'core'
            if Open == 'End':
                xi = (x1-x0)/(2**N)+x0
                fi = f(xi)
                if math.isnan(fi) == True:
                    raise Exception('Function value at ',xi,' isundefined, use ExtTrapAsymp() to perform integral over your desired period')
                fi *= h     
                f0 += fi    #Adding missing trapezoidal value
            else:
                # Creating new midpoint iteration
                f1 = f((x1-x0)/(2**(N+1))+x0)
                if math.isnan(f1) == True:
                    raise Exception('Function value at ',((x1-x0)/(2**(N+1))+x0),' is undefined, use ExtTrapAsymp() to perform integral over your desired period')
                f1 *= h
            #New trapezoidal/midpoint value for end of period also
            if Open == 'Start':
                xi = (x1-x0)*((2**N)-1)/(2**N)+x0
                fi = f(xi)
                if math.isnan(fi) == True:
                    raise Exception('Function value at ',xi,' is undefined, use ExtTrapAsymp() to perform integral over your desired period')
                fi *= h
                f0 += fi
            else:
                f2 = f((x1-x0)*((2**(N+1))-1/(2**(N+1)))+x0)
                if math.isnan(f2) == True:
                    raise Exception('Function value at ',((x1-x0)*((2**(N+1))-1/(2**(N+1)))+x0),' is undefined, use ExtTrapAsymp() to perform integral over your desired period')
                f2 *= h
            #Summing together different component for new estimate
            I2 = f0+f1+f2+F1+F2
            diff = np.abs((I2-I1)/I1)
            N += 1      #Updating number of iterations
            #Iteration Check
            if n == N:
                if Nreturn == True:
                    return I2, N
                else:
                    return I2
        if Nreturn == True:
            return I2, N
        else:
            return I2
    
def ExtTrapAsymp(f, x0, x1, Asymp, e = 0.01, n=all, listreturn = False):
    """
    Function for applying the extended trapezoidal rule when asymptotes
    within the integration period are known. These points must be given 
    as a list for the arg 'Asymp' and must include the start and/or end of 
    the function period if the user wishes these parts to be open also. 
    Otherwise uses same arguments as ExtTrap() so they can be immediately 
    implemented into the smaller iterations.
    
    N.B. Because each section has to be calculated seperately, there will
    be long runtimes for large lists of 'Asymp'
    """
    #First make sure it is a list
    if isinstance(Asymp,list) == False:
        raise Exception("Var, 'Asymp', must be given as a list")   
    S = 0           
    E = 0            
    l = len(Asymp)
    
    #Check to make sure integral isn't inverse so that sorting of Asymptotes aligns with x0 & x1
    if x0 > x1:
        Asymp.sort(reverse=True)
        #Check to ensure all asymptotes are within integration period, these are repeated
        if Asymp[0] > x0:
            raise Exception('Asymptotes must be within integration period')
        elif Asymp[0] == x0:
            S = 1        # Using dummy variable to indicate start is an asymptote
        if Asymp[len(l)] < x1:
            raise Exception('Asymptotes must be within integration period')
        elif Asymp[len(Asymp)] == x1:
            E = 1        # Using a seperate dummy variable to indicate end is also an asymptote
    #Repeated for typical order of x0 and x1 as well
    else:
        Asymp.sort()
        if Asymp[0] < x0:
            raise Exception('Asymptotes must be within integration period')
        elif Asymp[0] == x0:
            S = 1
        if Asymp[len(l)] > x1:
            raise Exception('Asymptotes must be within integration period')
        elif Asymp[len(Asymp)] == x1:
            E = 1
    # Calculating First Section of Trapezium rule, different depending on values of S
    if S == 1:
        I = ExtTrap(f, x0, Asymp[1] , e=e, Open = 'Both', n=n, Nreturn = listreturn)
        Ilist = [I]
    else:
        I = ExtTrap(f, x0, Asymp[0], e=e, Open = 'End', n=n, Nreturn = listreturn)
        Ilist = [I]
    # Running through each section of integral between the asymptotes
    for i in range(l-S-E):
        j = i+S
        I = ExtTrap(f, Asymp[j], Asymp[j+1], e=e, Open = 'Both', Nreturn = listreturn)
        Ilist.append(I)   #List of integration sections so they can be identified seperately
    # `Calculating final section
    if E == 1:
        I = ExtTrap(f, Asymp[l-1], x1, e=e, Open = 'Both', n=n, Nreturn = listreturn)
        Ilist.append(I)
    else:
        I = ExtTrap(f, Asymp[l], x1, e=e, Open = 'Start', n=n, Nreturn = listreturn)
    #Return lists if kwarg listreturn is satisfied
    if listreturn == True:
        Nlist = []
        # Taking original list and creating two seperate I and N lists
        for i in range(len(Ilist)):
            Nlist.append(Ilist[i][1])
            J = Ilist[i][0]
            Ilist[i] = J
        return sum(Ilist), Ilist, Nlist   #Note there are 3 returns if the list return kwarg is satisfied!
    else:
        return sum(Ilist)  #Returing cumulative sum of the integrations over the period

def ExtTrapAdd(f, x0, x1, oldI, oldN, e=0.001, Open = 'None', n = all, Nreturn = False):
    """
    Function for calculating additional Trapzeium rule iterations for
    a given previous calculation, provided its number of iterations is
    given. Can calculate to a new accuracy or for a set number of additional
    iterations depending on the inputs. 
    """
    # Checks to make sure input variables are compatible
    #if type(f) != 'function':
     #   raise Exception("Integrand function must be a function")
    if isinstance(x0,int) == True or isinstance(x0,float) == True:
        pass
    else:
        raise Exception('Start point must be a number')
    if isinstance(x1,int) == True or isinstance(x1,float) == True:
        pass
    else:
        raise Exception('End point must be a number')
    if x0 == x1:
        raise Exception('Start and End point must be different numbers')
    if isinstance(e,int) == True or isinstance(e,float) == True:
        pass
    else:
        raise Exception('Relative Accuracy must be a number')
    if isinstance(oldI,int) == True or isinstance(oldI,float) == True:
        pass
    else:
        raise Exception('Previous Trapzium calculation must be a number')
    if isinstance(oldN,int) != True:
        raise Exception('Previous number of Trapezium iterations must be an integer')
    if isinstance(n,int) == True and n>0:
        pass
    elif n==all:
        pass
    else:
        raise Exception("Number of Iterations must be an integer greater than 0, set to 'all' if you wish to calculate for a set accuracy, not a set number of iterations")
    N = oldN
    h = (x1-x0)/(2**(N-1))
    diff = e+1
    I2 = oldI
    if Open == 'None':
        #reused loop code from ExtTrap()
        while diff > e:
            I1 = I2
            I2 = 0
            h /= 2
            for i in range(2**(N-1)):
                xi = (x1-x0)*(1+2*i)/(2**N)+x0
                fi = f(xi)
                if math.isnan(fi) == True:
                    raise Exception('Function value at ',xi,' is undefined, restart and use ExtTrapAsymp() to perform integral over your desired period')
                I2 += fi
            I2 *= h
            I2 += I1/2                  
            diff = np.abs((I2-I1)/I1)  
            N += 1
            if n == (N-oldN):
                return I2
        if Nreturn == True:
            return I2, (N-oldN)
        else:
            return I2
    elif Open != 'Start' or Open != 'End' or Open != 'Both':
        raise Exception("kwarg 'Open' can only have values 'None', 'Start', 'End', or 'Both' ",Open," is not a valid value")
    else:
        #Solve Trapezium rule with open end(s)
        #Set-up to calculate old values of f0,f1,f2,F1,F2
        if Open == 'Start':
            f2 = 0
            F2 = 0
        else:
            f2 = f((x1-x0)*((2**(N+1))-1/(2**(N+1)))+x0)
            if math.isnan(f2) == True:
                raise Exception('Function value at ',(x1-x0)*((2**(N+1))-1/(2**(N+1)))+x0,'is undefined, Check to make sure oldI and oldN are correct')
            f2 *= h
            F2 = f((x1-x0)*((2**(N))-1/(2**(N)))+x0)
            if math.isnan(F2) == True:
                raise Exception('Function value at ',(x1-x0)*((2**(N))-1/(2**(N)))+x0,'is undefined, Check to make sure oldI and oldN are correct')
            F2 *= h/4
        if Open == 'End':
            f1 = 0
            F1 = 0
        else:
            f1 = f((x1-x0)/(2**(N+1))+x0)
            if math.isnan(f1) == True:
                raise Exception('Function value at ',(x1-x0)/(2**(N+1))+x0,'is undefined, Check to make sure oldI and oldN are correct')
            f1 *= h
            F1 = f((x1-x0)/(2**(N))+x0)
            if math.isnan(F1) == True:
                raise Exception('Function value at ',(x1-x0)/(2**(N))+x0,'is undefined, Check to make sure oldI and oldN are correct')
            F1 *= h/4    
        f0 = I2-(f1+F1+f2+F2)  #Calculate old f0 from other constants and I2/oldI
        
        #Copy of code from ExtTrap() for calculating further 'open' iterations
        while diff > e:
            I1 = I2
            f0 /= 2          
            f0 += (F1+F2)    
            F1 = f1/4        
            F2 = f2/4
            h /= 2
            for i in range(2**(N-1)):
                if i == 0 or i == (2**(N-1)-1):
                    pass
                xi = (x1-x0)*(1+2*i)/(2**N)+x0
                fi = f(xi)
                if math.isnan(fi) == True:
                    raise Exception('Function value at ',xi,' is undefined, use ExtTrapAsymp() to perform integral over your desired period')
                fi *= h
                f0 += fi     
            if Open == 'End':
                xi = (x1-x0)/(2**N)+x0
                fi = f(xi)
                if math.isnan(fi) == True:
                    raise Exception('Function value at ',xi,' isundefined, use ExtTrapAsymp() to perform integral over your desired period')
                fi *= h     
                f0 += fi    
            else:
                f1 = f((x1-x0)/(2**(N+1))+x0)
                if math.isnan(f1) == True:
                    raise Exception('Function value at ',(x1-x0)/(2**(N+1))+x0,' is undefined, use ExtTrapAsymp() to perform integral over your desired period')
                f1 *= h
            if Open == 'Start':
                xi = (x1-x0)*((2**N)-1)/(2**N)+x0
                fi = f(xi)
                if math.isnan(fi) == True:
                    raise Exception('Function value at ',xi,' is undefined, use ExtTrapAsymp() to perform integral over your desired period')
                fi *= h
                f0 += fi
            else:
                f2 = f((x1-x0)*((2**(N+1))-1/(2**(N+1)))+x0)
                if math.isnan(f2) == True:
                    raise Exception('Function value at ',(x1-x0)*((2**(N+1))-1/(2**(N+1)))+x0,' is undefined, use ExtTrapAsymp() to perform integral over your desired period')
                f2 *= h
            I2 = f0+f1+f2+F1+F2
            diff = np.abs((I2-I1)/I1)
            N += 1      
            if n == (N-oldN):
                if Nreturn == True:
                    return I2, (N-oldN)
                else:
                    return I2
        if Nreturn == True:
            return I2, (N-oldN)
        else:
            return I2

    
    
    

def ExtSimp(f, x0, x1, e=0.01, Open = 'None', n = all, Traphelp = [0,0], Nreturn = False):
    """
    Function for implementing the extended Simpson's rule by utilising
    previously calculated functions for Trapezium rule. It can also calculate
    the equivalent Simpsons rule for a given Trapezium rule and its
    number of iterations, however to do this one must provide a list of 
    the form [I,N] for the kwarg Traphelp.
    """
    N = 0
    oldN = 0
    if isinstance(Traphelp,list) != True:
        raise Exception('Traphelp must be given as a list of the form [I,N] from a previous use of the trapezium rule, set Nreturn = True for ExtTrap() to obtain this')
    if Traphelp[1] > 0:
        Tn = Traphelp[0]
        oldN = Traphelp[1]
    else:
        if isinstance(n,int) == True and n>1:
            pass
        elif n==all:
            pass 
        else:
            raise Exception("Number of Iterations must be an integer greater than 1, set to 'all' if you wish to calculate for a set accuracy, not a set number of iterations")
        Tn = ExtTrap(f, x0, x1, Open = 'None', n=2, Nreturn = False)
        N = 1
    if Open == 'None':
        Tn1 = ExtTrapAdd(f, x0, x1, Tn, (N+oldN), Open = Open, n = 1, Nreturn = False)
        Sn2 = (4/3)*Tn1-Tn/3
        if n == N:
            if Nreturn == True: #Optional Kwarg for returning N also
                return Sn2, N
            else:
                return Sn2
        diff = e+1
        while diff>e:
            Sn1 = Sn2
            N += 1
            Tn = Tn1
            Tn1 = ExtTrapAdd(f, x0, x1, Tn, (N+oldN), Open = Open, n = 1, Nreturn = False)
            Sn2 = (4/3)*Tn1-Tn/3
            diff = np.abs((Sn2-Sn1)/Sn1)
            if n == N:
                print(diff)
                if Nreturn == True: #Optional Kwarg for returning N also
                    return Sn1, N
                else:
                    return Sn1
        if Nreturn == True: #Optional Kwarg for returning N also
            return Sn2, N
        else:
            return Sn2 
    elif Open != 'Start' or Open != 'End' or Open != 'Both':
        raise Exception("kwarg 'Open' can only have values 'None', 'Start', 'End', or 'Both' ",Open," is not a valid value")
    else:
        #Find values of different variables within Integrand, these correspond to the variable names used in ExtTrap()
        I2 = Tn
        fs = []
        h = (x1-x0)/((2**(N+oldN))-1)
        for i in range (2**(N+oldN)):
            if i == 0:
                if Open == 'End':
                    f1 = 0
                    F1 = 0
                    xi = (x1-x0)*(1+2*i)/(2**(N+oldN+1))+x0
                    fi = f(xi)
                    if math.isnan(fi) == True:
                        raise Exception('Function value at ',xi,' is undefined, check to make sure oldI and oldN are correct or use ExtSimpAsymp() to perform integral over your desired period')
                    fi *= h
                    fs.append(fi)
                else:
                    f1 = f((x1-x0)/(2**(N+oldN+1))+x0)
                    if math.isnan(f1) == True:
                        raise Exception('Function value at ',(x1-x0)/(2**(N+oldN+1))+x0,' is undefined, check to make sure oldI and oldN are correct or use ExtSimpAsymp() to perform integral over your desired period')
                    f1 *= h
                    F1 = f((x1-x0)/(2**(N+oldN))+x0)
                    if math.isnan(F1) == True:
                        raise Exception('Function value at ',(x1-x0)/(2**(N+oldN))+x0,' is undefined, check to make sure oldI and oldN are correct or use ExtSimpAsymp() to perform integral over your desired period')
                    F1 *= h/3
            elif i == (2**(N+oldN)-1):
                if Open == 'Start':
                    f2 = 0
                    F2 = 0
                    xi = (x1-x0)*(1+2*i)/(2**(N+oldN+1))+x0
                    fi = f(xi)
                    if math.isnan(fi) == True:
                        raise Exception('Function value at ',xi,' is undefined, check to make sure oldI and oldN are correct or use ExtSimpAsymp() to perform integral over your desired period')
                    fi *= h
                    fs.append(fi)
                else:
                    f2 = f((x1-x0)*(1+2*i)/(2**(N+oldN+1))+x0)
                    if math.isnan(f1) == True:
                        raise Exception('Function value at ',(x1-x0)*(1+2*i)//(2**(N+oldN+1))+x0,' is undefined, check to make sure oldI and oldN are correct or use ExtSimpAsymp() to perform integral over your desired period')
                    f2 *= h 
                    F2 = f((x1-x0)*(1+i)/(2**(N+oldN))+x0)
                    if math.isnan(F2) == True:
                        raise Exception('Function value at ',(x1-x0)*(1+i)//(2**(N+oldN))+x0,' is undefined, check to make sure oldI and oldN are correct or use ExtSimpAsymp() to perform integral over your desired period')
                    F2 *= h/3 
            else:
                xi = (x1-x0)*(1+2*i)/(2**(N+oldN+1))+x0
                fi = f(xi)
                if math.isnan(fi) == True:
                    raise Exception('Function value at ',xi,' is undefined, check to make sure oldI and oldN are correct or use ExtSimpAsymp() to perform integral over your desired period')
                fi *= h
                fs.append(fi)
        Tn = I2-(f1+f2)  #Calculate which part is purely trapezoidal
        f0 = Tn-(F1+F2)
        #Obtain value of next trapezoidal interation of the purely trapezoidal part
        Tn1 = ExtTrapAdd(f, x0, x1, Tn, (N+oldN), Open = 'None', n = 1, Nreturn = False)
        
        #And then determine equivalent Simpson iteration and add midpoint rule parts
        Sn2 = ((4/3)*Tn1-Tn/3)+(f1+f2) 
        #Apply Iteration Check
        if n == N:
            if Nreturn == True: 
                return Sn2, N
            else:
                return Sn2
        diff = e+1
        #Accuracy Check for the following iterations
        while diff>e:
            Sn1 = Sn2
            f0 /= 2
            f0 += (F1+F2)
            f0 += (sum(fs))/4
            F1 = f1/6
            F2 = f2/6
            h /= 2
            for i in range(2**(N+oldN)):
                if i == 0:
                    if Open == 'End':
                        xi = (x1-x0)*(1+2*i)/(2**(N+oldN+1))+x0
                        fi = f(xi)
                        if math.isnan(fi) == True:
                            raise Exception('Function value at ',xi,' is undefined, check to make sure oldI and oldN are correct or use ExtSimpAsymp() to perform integral over your desired period')
                        fi *= h
                        fs.append(fi)
                    else:
                        f1 = f((x1-x0)/(2**(N+oldN+1))+x0)
                        if math.isnan(f1) == True:
                            raise Exception('Function value at ',(x1-x0)/(2**(N+oldN+1))+x0,' is undefined, check to make sure oldI and oldN are correct or use ExtSimpAsymp() to perform integral over your desired period')
                        f1 *= h
                elif i == (2**(N+oldN)-1):
                    if Open == 'Start':
                        xi = (x1-x0)*(1+2*i)/(2**(N+oldN+1))+x0
                        fi = f(xi)
                        if math.isnan(fi) == True:
                            raise Exception('Function value at ',xi,' is undefined, check to make sure oldI and oldN are correct or use ExtSimpAsymp() to perform integral over your desired period')
                        fi *= h
                        fs.append(fi)
                    else:
                        f2 = f((x1-x0)*(1+2*i)/(2**(N+oldN+1))+x0)
                        if math.isnan(f1) == True:
                            raise Exception('Function value at ',(x1-x0)*(1+2*i)//(2**(N+oldN+1))+x0,' is undefined, check to make sure oldI and oldN are correct or use ExtSimpAsymp() to perform integral over your desired period')
                        f2 *= h 
                else:
                    xi = (x1-x0)*(1+2*i)/(2**(N+oldN+1))+x0
                    fi = f(xi)
                    if math.isnan(fi) == True:
                        raise Exception('Function value at ',xi,' is undefined, check to make sure oldI and oldN are correct or use ExtSimpAsymp() to perform integral over your desired period')
                    fi *= h
                    fs.append(fi)
            Sn2 = f0+f1+f2+F1+F2+sum(fs)
            N += 1
            diff = np.abs((Sn2-Sn1)/Sn1)
            if n == N:
                if Nreturn == True: 
                    return Sn2, N
                else:
                    return Sn2
        if Nreturn == True: 
            return Sn2, N
        else:
            return Sn2 

def ExtSimpAsymp(f, x0, x1, Asymp, e = 0.01, n=all, listreturn = False):
    """
    Function for applying the extended Simpson's rule when asymptotes
    within the integration period are known. These points must be given 
    as a list for the arg 'Asymp' and must include the start and/or end of 
    the function period if the user wishes these parts to be open also. 
    Unlike ExtSimp() this does not have an option to convert from the trapezoid
    equivalent to reduce calculation, if this is desired it can be done manually 
    using ExtSimp() for each section individually, which have been provided 
    from ExtTrapAsymp() and given by setting listreturn to 'True'
    """
    #First make sure it is a list
    if isinstance(Asymp,list) == False:
        raise Exception("Var, 'Asymp', must be given as a list")   
    S = 0           
    E = 0            
    l = len(Asymp)
    
    #Check to make sure integral isn't inverse so that sorting of Asymptotes aligns with x0 & x1
    if x0 > x1:
        Asymp.sort(reverse=True)
        #Check to ensure all asymptotes are within integration period, these are repeated
        if Asymp[0] > x0:
            raise Exception('Asymptotes must be within integration period')
        elif Asymp[0] == x0:
            S = 1        # Using dummy variable to indicate start is an asymptote
        if Asymp[len(l)] < x1:
            raise Exception('Asymptotes must be within integration period')
        elif Asymp[len(Asymp)] == x1:
            E = 1        # Using a seperate dummy variable to indicate end is also an asymptote
    #Repeated for typical order of x0 and x1 as well
    else:
        Asymp.sort()
        if Asymp[0] < x0:
            raise Exception('Asymptotes must be within integration period')
        elif Asymp[0] == x0:
            S = 1
        if Asymp[len(l)] > x1:
            raise Exception('Asymptotes must be within integration period')
        elif Asymp[len(Asymp)] == x1:
            E = 1
    # Calculating First Section of Trapezium rule, different depending on values of S
    if S == 1:
        I = ExtSimp(f, x0, Asymp[1] , e=e, Open = 'Both', n=n, Nreturn = listreturn)
        Ilist = [I]
    else:
        I = ExtSimp(f, x0, Asymp[0], e=e, Open = 'End', n=n, Nreturn = listreturn)
        Ilist = [I]
    # Running through each section of integral between the asymptotes
    for i in range(l-S-E):
        j = i+S
        I = ExtSimp(f, Asymp[j], Asymp[j+1], e=e, Open = 'Both', Nreturn = listreturn)
        Ilist.append(I)   #List of integration sections so they can be identified seperately
    # `Calculating final section
    if E == 1:
        I = ExtSimp(f, Asymp[l-1], x1, e=e, Open = 'Both', n=n, Nreturn = listreturn)
        Ilist.append(I)
    else:
        I = ExtSimp(f, Asymp[l], x1, e=e, Open = 'Start', n=n, Nreturn = listreturn)
    #Return lists if kwarg listreturn is satisfied
    if listreturn == True:
        Nlist = []
        # Taking original list and creating two seperate I and N lists
        for i in range(len(Ilist)):
            Nlist.append(Ilist[i][1])
            J = Ilist[i][0]
            Ilist[i] = J
        return sum(Ilist), Ilist, Nlist   #Note there are 3 returns if the list return kwarg is satisfied!
    else:
        return sum(Ilist)  #Returing cumulative sum of the integrations over the period

    
#pdf function for uniform sampling

def SampleRandom(x0,x1,Shots, pdf=False):
    """
    Function for obtaining a uniform sample of random numbers between
    x0 and x1 with the option to return their pdf also using
    the 'pdf' kwarg.
    """
    xlist = []        #list to contain all data points
    #Creating the data points
    for i in range(Shots):
        x = random.uniform(x0,x1)
        xlist.append(x)
    if pdf == True:        #Also returning a pdf function for uniform sampling if requested by kwarg
        def uniform(x):
            return 1/(x1-x0)
        return xlist, uniform
    else:                  #Default is just a list
        return xlist



def SamplePdf(x0,x1,Shots, pdf = False):
    """
    Function for obtaining a non-uniform sample between x0 and x1
    with a weighting determined by the Transformation method being
    applied using a pdf of the form 'y = Ax+b' where A is within a 
    specific range to ensure normalisation and b can be defined from
    A.
    
    n.b. A is a global variable here so it must be chosen before this sampler is used
    """
    #Check A is within the allowed range
    if A <= -2*((x1-x0)**-2):
        raise Exception('The minimum value of A for this function is',-2*((x1-x0)**-2),' for your given period')
    elif A >= 2*((x1-x0)**-2):
        raise Exception('The maximum value of A for this function is',2*((x1-x0)**-2),' for your given period')
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
    
def MonteCarloIntegr(f, x0, x1, Shots, Sampler, GiveErr = False):
    """
    Function for performing Monte Carlo Integration using 
    Importance Sampling on a function f. The sampler must
    be provided seperately and be of the same form as the samplers
    used here (args, kwargs, returns, etc.). Function also 
    has option to provide the given error associated with the result
    """
    # Input Variable Checks
   # if type(f) != 'function':
 #       raise Exception("Integrand function must be a function")
    if isinstance(x0,int) == True or isinstance(x0,float) == True:
        pass
    else:
        raise Exception('Start point must be a number')
    if isinstance(x1,int) == True or isinstance(x1,float) == True:
        pass
    else:
        raise Exception('End point must be a number')
    if x0 == x1:
        raise Exception('Start and End point must be different numbers')
    Psamples, pdf = Sampler(x0, x1, Shots, pdf = True)
    Qsamples = []
    for i in range(Shots):
        fi = f(Psamples[i])
        pi = pdf(Psamples[i])
        qi = fi/pi
        Qsamples.append(qi)
    ExpQ = (sum(Qsamples))/Shots
    if GiveErr == True:
        Varsamples = []
        for i in range(Shots):
            Vi = (Qsamples[i]-ExpQ)**2
            Varsamples.append(Vi)
        Var = (sum(Varsamples))/(Shots*(Shots-1))
        Err = Var**0.5
        return ExpQ, Err
    else:
        return ExpQ


def ThreeDMonteCarlo(f, Boundaries, Shots, Samplers, GiveErr = False):
    """
    Function for applying Monte Carlo Integration to a 3D function
    for a given cubic volume fully containing the function whose boundaries
    are given by a list in its aptly named arg in the format[x0,x1,y0..]. 
    Seperate Samplers can be used for each of the three dimensions however 
    they must be given as a list of functions, and must satisfy 
    MonteCarloIntegr()'s requirements for Samplers.
    
    N.B. Ordering Convention for the inputs is x, then y, then z and the function treats it as such
    """
    # Input Variable Checks
    #if type(f) != 'function':
     #   raise Exception("Integrand function must be a function")
    if len(Boundaries) != 6:
        raise Exception('There should only be 6 numbers here to represent the max and min in each dimension')
    for i in range(len(Boundaries)):    
        if isinstance(Boundaries[i],int) == True or isinstance(Boundaries[i],float) == True:
            pass
        else:
            raise Exception('Start and End points must be numbers')
        if i%2 == 1 and Boundaries[i] == Boundaries[i-1]:
            raise Exception('Start and End point musts be different numbers and correctly ordered')
    #Obtain the x, y and z values
    xsamples, xpdf = Samplers[0](Boundaries[0], Boundaries[1], Shots, pdf = True)
    ysamples, ypdf = Samplers[1](Boundaries[2], Boundaries[3], Shots, pdf = True)
    zsamples, zpdf = Samplers[2](Boundaries[4], Boundaries[5], Shots, pdf = True)
    Qsamples = []
    for i in range(Shots):
        fi = f(xsamples[i],ysamples[i],zsamples[i])
        pi = xpdf(xsamples[i])*ypdf(ysamples[i])*zpdf(zsamples[i])
        qi = fi/pi
        Qsamples.append(qi)
    ExpQ = (sum(Qsamples))/Shots
    if GiveErr == True:
        Varsamples = []
        for i in range(Shots):
            Vi = (Qsamples[i]-ExpQ)**2
            Varsamples.append(Vi)
        Var = (sum(Varsamples))/(Shots*(Shots-1))
        Err = Var**0.5
        return ExpQ, Err
    else:
        return ExpQ
    
def ThreeDNewtonCoates(f, Boundaries, e=0.01, n=all, Nreturn = False):
    """
    Function for applying the Newton-Coates Method, the extended Trapezium
    rule to a 3D function for a given cubic volume. 
    """
    # Input Variable Checks
   # if type(f) != 'function':
    #    raise Exception("Integrand function must be a function")
    if len(Boundaries) != 6:
        raise Exception('There should only be 6 numbers here to represent the max and min in each dimension')
    for i in range(len(Boundaries)):    
        if isinstance(Boundaries[i],int) == True or isinstance(Boundaries[i],float) == True:
            pass
        else:
            raise Exception('Start and End points must be numbers')
        if i%2 == 1 and Boundaries[i] == Boundaries[i-1]:
            raise Exception('Start and End point musts be different numbers and correctly ordered')
    if isinstance(e,int) == True or isinstance(e,float) == True:
        pass
    else:
        raise Exception('Relative Accuracy must be a number')
    if isinstance(n,int) == True and n>1:
        pass
    elif n==all:
        pass 
    else:
        raise Exception("Number of Iterations must be an integer greater than 1, set to 'all' if you wish to calculate for a set accuracy, not a set number of iterations")
    N = 1                     
    #Input values from Boundaries into variables
    x0 = Boundaries[0]
    x1 = Boundaries[1]
    y0 = Boundaries[2]
    y1 = Boundaries[3]
    z0 = Boundaries[4]
    z1 = Boundaries[5]
    hx = x1 - x0
    hy = y1 - y0
    hz = z1 - z0
    hprod = hx*hy*hz
    I1 = 0
    #Creating Corner Values
    for i in range(2):
        for j in range(2):
            for k in range(2):
                fi = f(Boundaries[i],Boundaries[j+2],Boundaries[k+4])
                if math.isnan(fi) == True:
                    raise Exception('Function value at ',(Boundaries[i],' ',Boundaries[j+2],' ',Boundaries[k+4]),' is undefined')
                I1 += fi
    I1 *= hprod/8
    N += 1
    # Evaluate Second trapzium rule iteration
    hprod /= 8    #Normally would be factor of 2, but as it needs to be applied for each dimension becomes 8
    I2 = 0
    for i in range(3):                       #Run through each x value
        xi = ((x1-x0)*i/2)+x0
        for j in range(3):                   #And its possible y values
            yi = ((y1-y0)*j/2)+y0
            for k in range(3):               #And their possible z values
                zi = ((z1-z0)*k/2)+z0
                if i%2 == 0 and j%2 == 0 and k%2 == 0:    #If 2 of i, j, or k is even it means this is a repeated value as fraction can be simplified
                    pass                                  #Therefore don't recalculate
                else:
                    fi = f(xi,yi,zi)
                    if math.isnan(fi) == True:
                        raise Exception('Function value at ',xi,' ',yi,' ',zi,' is undefined')
                    if i == 0 or i == 2:
                        fi /= 2
                    if j == 0 or j == 2:
                        fi /= 2
                    if k == 0 or k == 2:
                        fi /= 2
                    I2 += fi
    I2 *= hprod
    I2 += I1/8
    diff = np.abs((I2-I1)/I1)
    # Optional variable for bypassing accuracy check and only completing set number of iterations instead
    if n == N:
        if Nreturn == True: #Optional Kwarg for returning N also
            return I2, N
        else:
            return I2
    # Create loop adding additional iterations until accuracy check is satisfied
    while diff > e:
        I1 = I2
        I2 = 0
        hprod /= 8
        # Create all the new values for this iteration
        for i in range((2**N)+1):
            xi = ((x1-x0)*i/(2**N))+x0
            for j in range((2**N)+1):
                yi = ((y1-y0)*j/(2**N))+y0
                for k in range((2**N)+1):
                    zi = ((z1-z0)*k/(2**N))+z0
                    if i%2 == 0 and j%2 == 0 and k%2 == 0:
                        pass
                    else:
                        fi = f(xi,yi,zi)
                        if math.isnan(fi) == True:
                            raise Exception('Function value at ',xi,' ',yi,' ',zi,' is undefined')
                        if i == 0 or i == (2**N)-1:
                            fi /= 2
                        if j == 0 or j == (2**N)-1:
                            fi /= 2
                        if k == 0 or k == (2**N)-1:
                            fi /= 2
                        I2 += fi
        I2 *= hprod
        I2 += I1/8                  #New integral estimation
        diff = np.abs((I2-I1)/I1)   #New accuracy value
        N += 1
        # If iteration check is reached before accuracy check, result is returned
        if n == N:
            return I2, N
    # Once accuracy check is completed result is provided along with number of iterations if requested
    if Nreturn == True:
        return I2, N
    else:
        return I2
                

