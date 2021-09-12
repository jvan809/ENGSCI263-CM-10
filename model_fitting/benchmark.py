
# imports
import numpy as np
from matplotlib import pyplot as plt
import scipy as sp
from scipy.integrate import solve_ivp
from scipy.optimize.nonlin import nonlin_solve

# Why is this neccessary?
if __name__ == "__main__":
    from cm_solver_functions import *
else:
    from model_fitting.cm_solver_functions import *


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
    ax1 = ax.twinx()
    # pressure
    ax.plot(t_analytical, pre_analytical , 'r-', label='Pressure Analytical sol')
    ax.plot(time, pressure, 'bx', label = 'Pressure Numerical sol')

    # temp
    ax1.plot(t_analytical, tmp_analytical , 'g-', label='Temperature Analytical Sol')
    ax1.plot(time, temp, 'kx', label = 'Temperature Numerical Sol')

    ax.legend()
    ax.legend(bbox_to_anchor=(0.92,1.15), loc = 'upper center', frameon = True)
    ax1.legend(bbox_to_anchor=(0.07,1.15), loc = 'upper center', frameon = True)
    ax.set_xlabel("time")
    ax.set_ylabel("Pressure")
    ax1.set_ylabel("Temperature")
    ax.set_title("Benchmarking")
    plt.show()



if __name__ == "__main__":
    
    plot_benchmark()
