#MOODY DIAGRAM USING NEWTON-RAPHSON AND SECANT METHOD

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
    iteration_counter = 0
    while abs(f_value) > eps and iteration_counter < 10000:
        try:
            x0=x
            x = x - float(f_value)/deriv_f(kd,x,Re) #finding xn from xn-1
        except ZeroDivisionError:
            print ("Error! - derivative zero for x = ", x)
            sys.exit(1)     # Abort with error

        f_value = f(kd,x,Re)
        iteration_counter += 1

    # Here, either a solution is found, or too many iterations
    if abs(f_value) > eps:
        iteration_counter = -1
    return x, iteration_counter

#function to solve function f by secant method
#where deriv_f is derivative of f wrt friction factor x
#kd is ratio of relative roughness to diameter of pipe
#Re is reynolds number
#eps is to come out of loop when f value becomes smaller than epsilon
def Secant(f,kd,Re,x0,x1, eps):
    f0_value = f(kd,x0,Re)
    f1_value = f(kd,x1,Re)
    iteration_counter = 0
    while abs(f1_value) > eps and iteration_counter < 10000:
        try:
            p=x1       #p is temp var to store original value of xn-1
            x1 = x1 - float(f1_value)*(x1-x0)/((f1_value)-(f0_value)) #finding xn from xn-1 and xn-2
            x0=p #putting value of xn-1 in x0
        except ZeroDivisionError:
            print ("Error! - derivative zero for x1 = ", x1)
            sys.exit(1)     # Abort with error

        f1_value = f(kd,x1,Re)
        f0_value = f(kd,x0,Re)
        iteration_counter += 1

    # Here, either a solution is found, or too many iterations
    if abs(f1_value) > eps:
        iteration_counter = -1
    return x1, iteration_counter


#function f Here Colebrook Eqn
def f(kd,x,Re):
    return 1/(x**(1/2)) + 0.86858896*  math.log(kd /3.7 + ( 2.51/(Re * x**(1/2))))

#first derivative of f wrt x
def deriv_f(kd,x,Re):
    return  -1*(1/(2 * x**(3/2)) + (2.1801583/(Re*kd/3.7*x**(3/2)+2.51*x)))

 
sol_newton= np.zeros([1,0],dtype=float)     #initiate array that contains soln. for given Re newton 
iteration_newton=np.zeros([1,0],dtype=float)     #initiate array that contains no. of iterations newton
reynold_newton=np.zeros([1,0])              #initiate array that contains Re values newton
sol_lin= np.zeros([1,0],dtype=float)     #initiate array that contains soln. for given Re linear
reynold_lin=np.zeros([1,0])              #initiate array that contains Re values linear
sol_secant= np.zeros([1,0],dtype=float)     #initiate array that contains soln. for given Re secant
iteration_secant=np.zeros([1,0],dtype=float)     #initiate array that contains no. of iterations secant
reynold_secant=np.zeros([1,0])              #initiate array that contains Re values secant

#Friction factor in case of laminar flow
#f=64/Re values being inserted for Re 240 to 2000
re=250    
while re< 2300:
  sol_lin=np.append(sol_lin,[[64/re]])
  reynold_lin=np.append(reynold_lin,[[re]]) 
  re+=200  

kdv=0.000005 #value of relative pipe roughness to diameter
#loop to solve for various kd values
while kdv <= 0.005:
   i = 5000
#loop to solve for various Re value 5000 to 10^8
   while i < 100000000:
      solution_nr, no_iterations_nr = Newton(f, deriv_f,kdv,i,64/i, eps=1.0e-6) 
      x=64/i - float(f(kdv,64/i,i)/deriv_f(kdv,64/i,i)) # second initial value for secant method
      solution_sec, no_iterations_sec = Secant(f,kdv,i,64/i,x, eps=1.0e-6) 
      #store soln and no of itertion corresponding to given Re and kd values  
      #here 64/re is taken as initial guess for friction factor  
      if no_iterations_nr > 0:    # Solution found
        sol_newton=np.append(sol_newton,[[solution_nr]])         #Store soln in array
        iteration_newton=np.append(iteration_newton,[[1+2*no_iterations_nr]])#Store iterations in array
        reynold_newton=np.append(reynold_newton,[[i]])        #Store corresp. Re values
      else:
        print ("Solution not found!")   
        continue
      if no_iterations_sec > 0:    # Solution found
        sol_secant=np.append(sol_secant,[[solution_sec]])         #Store soln in array
        iteration_secant=np.append(iteration_secant,[[1+2*no_iterations_sec]])#Store iterations in array
        reynold_secant=np.append(reynold_secant,[[i]])        #Store corresp. Re values
      else:
        print ("Solution not found!")   
        continue
      i=i+5000                                    #Interval gap of Re values
   kdv=kdv*10                                     #Scale at which kd value varies

#Plotting the moody diagram using Newton Raphson      
plt.xscale("log")
plt.yscale("log")
plt.scatter(reynold_lin,sol_lin, s=4, c="b")
plt.scatter(reynold_newton,sol_newton, s=4, c="b")
plt.title('Moody Diagram')
plt.xlabel('Re (Reynolds Number on log scale)')
plt.ylabel('f (Friction Factor)')
plt.show()
#Plotting the moody diagram using secant      
plt.xscale("log")
plt.yscale("log")
plt.scatter(reynold_lin,sol_lin, s=4, c="r")
plt.scatter(reynold_secant,sol_secant, s=4, c="r")
plt.title('Moody Diagram')
plt.xlabel('Re (Reynolds Number on log scale)')
plt.ylabel('f (Friction Factor)')
plt.show()
