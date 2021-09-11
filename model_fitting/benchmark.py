# ENGSCI263: Lab Exercise 2
# lab2.py

# PURPOSE:
# IMPLEMENT a lumped parameter model and CALIBRATE it to data.

# PREPARATION:
# Review the lumped parameter model notes and use provided data from the kettle experiment.

# SUBMISSION:
# - Show your calibrated LPM to the instructor (in Week 3).

# imports
import numpy as np
from matplotlib import pyplot as plt
import scipy as sp
from scipy.integrate import solve_ivp
from scipy.optimize.nonlin import nonlin_solve
from cm_solver_functions import *


def plot_benchmark():
    ''' Compare analytical and numerical solutions.

        Parameters:
        -----------
        none

        Returns:
        --------
        none

        Notes:
        ------
        This function called within if __name__ == "__main__":

        It should contain commands to obtain analytical and numerical solutions,
        plot these, and either display the plot to the screen or save it to the disk.
        
    '''
    t_analytical = np.zeros(1000)
    pre_analytical = np.zeros(1000)
    tmp_analytical = np.zeros(1000)
    for i in range(0,1000):
        t_analytical[i] = i/100
        pre_analytical[i] = 0.1*np.exp(-i/100)

    for i in range(0, 1000):
        t_analytical[i] = i/100
        tmp_analytical[i]= 1-(0.1105*np.exp(-0.1*np.exp(-i/100)-i/100))

    time, pressure, temp = ode_solve(model, 0, 10 ,0.1 ,0.9 ,[-1,0,1,1,1,1,1,1,1])

    f, ax = plt.subplots(1,1)
    # pressure
    ax.plot(t_analytical, pre_analytical , 'r-', label='Analytical')
    ax.plot(time, pressure, 'bx', label = 'Numerical')

    # temp
    ax.plot(t_analytical, tmp_analytical , 'r-', label='Analytical')
    ax.plot(time, temp, 'bx', label = 'Numerical')

    leg = ax.legend()
    ax.legend(loc = 'upper right', frameon = True)
    plt.show()

    return None

def improved_euler(f, t0, t1, dt, x0, pars, q=None):
    ''' Solve an ODE numerically.

        Parameters:
        -----------
        f : callable
            Function that returns dxdt given variable and parameter inputs.
        t0 : float
            Initial time of solution.
        t1 : float
            Final time of solution.
        dt : float
            Time step length.
        x0 : float
            Initial value of solution.
        pars : array-like
            List of parameters passed to ODE function f.
        q : callable
            function giving the flow at time=t

    '''
    t = np.arange(t0, t1, dt) # create array of times
    x = np.zeros(len(t)) # initialise solution output
    x[0] = x0 # set initial value
    
    if q is None: q = -np.ones(len(t)) # this function gives the value of q at time t

    for i in range(len(t)-1):
        k1 = f(t[i], x[i], q[i], *pars)            # /\
        k2 = f(t[i+1], x[i] + k1*dt, q[i+1], *pars) # equations for the improved Euler method
        x[i+1] = x[i] + dt*(k1+k2)/2                      # \/
    return t,x

    
def plot_kettle_model():
    ''' Plot the kettle LPM over top of the data.

        Parameters:
        -----------
        none

        Returns:
        --------
        none

        Notes:
        ------
        This function called within if __name__ == "__main__":

        It should contain commands to read and plot the experimental data, run and 
        plot the kettle LPM for hard coded parameters, and then either display the 
        plot to the screen or save it to the disk.
    '''
    a = 10.8/(1000*4200*0.005)
    b = 0.0006
    pars = [0,a,b,22]
    time, pressure, temp = ode_solve(ode_model, 0,10,0.1,0,[-1,1,1,1,1,1,1,1,1])

    plt.subplot(1, 1, 1)
    plt.plot(time, temp, label='Model')
    plt.legend()
    plt.xlabel('Time, s')
    plt.ylabel('Temperature, celsius')
    plt.show()

if __name__ == "__main__":
    # plot_kettle_model()
    plot_benchmark()
