#COMPARING NEWTON-RAPHSON AND SECANT METHOD ERROR VS ITERATION FOR COLEBROOK EQN

import math
import numpy as np
import matplotlib.pyplot as plt

#function to solve function f by newton-raphson method
#where deriv_f is derivative of f wrt friction factor x
#kd is ratio of relative roughness to diameter of pipe
#Re is reynolds number
#eps is to come out of loop when f value becomes smaller than epsilon
def Newton(f, deriv_f,kd,Re, x, eps):
    f_value = f(kd,x,Re)
    f0= np.zeros([1,0],dtype=float)     #initiate array that contains f0 for given Re newton 
    i= np.zeros([1,0],dtype=float)     #initiate array that contains i for given Re newton 
    iteration_counter = 0
    while abs(f_value) > eps and iteration_counter < 10000:
        try:
            x = x - float(f_value)/deriv_f(kd,x,Re) #finding xn from xn-1
            f0=np.append(f0,[[f_value]])         #Store soln in array            
            i=np.append(i,[[iteration_counter]])         #Store soln in array            
        except ZeroDivisionError:
            print ("Error! - derivative zero for x = ", x)
            sys.exit(1)     # Abort with error

        f_value = f(kd,x,Re)
        iteration_counter += 1
    f0=np.append(f0,[[f_value]])         #Store soln in array            
    i=np.append(i,[[iteration_counter]])         #Store soln in array            
       
    # Here, either a solution is found, or too many iterations
    if abs(f_value) > eps:
        iteration_counter = -1
    return f0, i

#function to solve function f by secant method
#where deriv_f is derivative of f wrt friction factor x
#kd is ratio of relative roughness to diameter of pipe
#Re is reynolds number
#eps is to come out of loop when f value becomes smaller than epsilon
def Secant(f,kd,Re,x0,x1, eps):
    f0_value = f(kd,x0,Re)
    f1_value = f(kd,x1,Re)
    f0= np.zeros([1,0],dtype=float)     #initiate array that contains f0 for given Re newton 
    i= np.zeros([1,0],dtype=float)     #initiate array that contains i for given Re newton 
    iteration_counter = 0   
    f0=np.append(f0,[[f0_value]])         #Store soln in array
    i=np.append(i,[[iteration_counter]])         #Store soln in array            
    iteration_counter = 1   
    while abs(f1_value) > eps and iteration_counter < 10000:
        try:
            p=x1       #p is temp var to store original value of xn-1
            x1 = x1 - float(f1_value)*(x1-x0)/((f1_value)-(f0_value)) #finding xn from xn-1 and xn-2
            x0=p #putting value of xn-1 in x0
            f0=np.append(f0,[[f1_value]])         #Store soln in array
            i=np.append(i,[[iteration_counter]])         #Store soln in array            
        except ZeroDivisionError:
            print ("Error! - denominator zero for x1 = ", x1)
            sys.exit(1)     # Abort with error

        f1_value = f(kd,x1,Re)
        f0_value = f(kd,x0,Re)
        iteration_counter += 1

    f0=np.append(f0,[[f1_value]])         #Store soln in array
    i=np.append(i,[[iteration_counter]])         #Store soln in array            
       # Here, either a solution is found, or too many iterations
    if abs(f1_value) > eps:
        iteration_counter = -1
    return f0, i


#function f Here Colebrook Eqn
def f(kd,x,Re):
    return 1/(x**(1/2)) + 0.86858896*  math.log(kd /3.7 + ( 2.51/(Re * x**(1/2))))

#first derivative of f wrt x
def deriv_f(kd,x,Re):
    return  -1*(1/(2 * x**(3/2)) + (2.1801583/(Re*kd/3.7*x**(3/2)+2.51*x)))

 
re=4000    
kdv=0.000005 #value of relative pipe roughness to diameter
f_nr, no_iterations_nr = Newton(f, deriv_f,kdv,re,64/re, eps=1.0e-8) 
x=64/re - float(f(kdv,64/re,re)/deriv_f(kdv,64/re,re)) # second initial value for secant method
f_sec, no_iterations_sec = Secant(f,kdv,re,64/re,x, eps=1.0e-8) 

#Plotting the error vs iteration to compare Newton Raphson and Secant Method      
plt.yscale("log")
plt.xticks(np.arange(0, 20, 1))
plt.scatter(no_iterations_nr,f_nr, s=5, c="b",label="Newton-Raphson")
plt.scatter(no_iterations_sec,f_sec, s=5, c="r",label="Secant")
plt.legend()
plt.title('Comparing the rate of convergence')
plt.xlabel('Value of i')
plt.ylabel('Value Of F(x_i)')
plt.show()
