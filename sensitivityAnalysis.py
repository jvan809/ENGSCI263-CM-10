import numpy as np
import matplotlib.pyplot as plt 
from model_fitting.cm_solver_functions import improved_euler
from unit_test import dummyModel
from unit_test import solve_ode_test

def sensitivityAnalysis():
    '''
    Demonstrates the improved_euler() function defined in cm_solver_functions.py
    has instability when the step size is high.
    '''

    # Defining step sizes to calculate and graph
    steps = np.linspace(start = 0.01, stop = 0.5, num = 50)

    # The following loop goes through every step in steps, and runs the 
    # ode_solve function using that step size. The model 
    # the function is solving is borrowed from unit_test.py

    # Creating empty arrays to store the final result of each loop run in
    finalPressureResults = []

    for step in steps:
        # Calling the Improved Euler's function to calculate pressure
        [time, pressure] = improved_euler(dummyModel, t0 = 0, t1 = 3, dt = step, x0 = 50, pars = [1,1,1])
        # Storing the final results
        finalPressureResults = np.append(finalPressureResults, pressure[-1])


    # Plotting
    plt.scatter(steps, finalPressureResults)
    plt.ylabel("Final value")
    plt.xlabel("Step size")
    plt.title("Final value of solution vs step size")
    plt.show()