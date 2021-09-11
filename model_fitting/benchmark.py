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


def pressure_analytical(t):
    '''
    Analytical pressure solution for parameters [-1,0,1,1,1,1,1,1,1]

    inputs:
    -------
    t : float
        time at which the solution is evaluated

    outputs:
    --------
    P : float
        Pressure at time t
    '''  
    return 0.1*np.exp(-t)

def temperature_analytical(t):
    '''
    Analytical temperature solution for parameters [-1,0,1,1,1,1,1,1,1]

    inputs:
    -------
    t : float
        time at which the solution is evaluated

    outputs:
    --------
    T : float
        Temperature at time t
    '''
    return  1-(0.11051709180756478*np.exp(-0.1*np.exp(-t)-t))


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
    t_analytical = np.linspace(0,10,1000)
    pre_analytical = pressure_analytical(t_analytical)
    tmp_analytical = temperature_analytical(t_analytical)
    
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
