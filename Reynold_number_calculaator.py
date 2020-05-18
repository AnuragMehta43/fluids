#FRICION FACTOR USING NEWTON-RAPHSON AND SECANT METHOD FOR GIVEN REYNOLDS NUMBER
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
            print ("Error! - denominator zero for x1 = ", x1)
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
while True:
        try:
            Re = int(input("Reynolds Number(integer values only):" ))
            break
        except:
            print("That's not a valid option!")
kdv=0.000005
if Re<=0:
  print("Invalid Input")
elif Re<=2300: 
  print("Laminar Flow Friction Factor=",64/Re)
elif Re<=4000 and Re>2300  :
  print("Transition Phase")
elif Re>=4000:
  solution_nr, no_iterations_nr = Newton(f, deriv_f,kdv,Re,64/Re, eps=1.0e-8) 
  print("Turbulent Flow ")
  if no_iterations_nr > 0:    # Solution found
    print("Friction Factor using Newton-Raphson=",solution_nr,"found in",no_iterations_nr,"iterations")
  else:
    print ("Solution not found using Newton-Raphson!")   
  x=64/Re - float(f(kdv,64/Re,Re)/deriv_f(kdv,64/Re,Re)) # second initial value for secant method
  solution_sec, no_iterations_sec = Secant(f,kdv,Re,64/Re,x, eps=1.0e-8) 
  if no_iterations_sec > 0:    # Solution found
    print("Friction Factor using Secant Method=",solution_sec,"found in",no_iterations_sec,"iterations")
  else:
    print ("Solution not found using Secant Method=!")   
else:
  print("Error")