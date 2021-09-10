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

def ode_model(t, X, q_stm, q_out, Tstm, aP, bP, P0, M0, T0, bT):
    """
    The model giving the derivaitves for the differential equation
    
    inputs:
    -------
    t : float
        time at which derivative is evaluated
    X : numpy array 
        = [P,T]
        P : double
            Pressure of the system (Pa)
        T : double
            Temperature of the system  (deg C)
    q_stm : callable or double
        Flow of steam into the system as a function of time(m^3/s)
    q_out : callable or double
        Flow of water and bitumen out of the system as a function of time(m^3/s)
    Tstm : double
        Injected Steam Temperature (deg C)
    aP, bP : double
        Lumped Parameters for Pressure Equation (units Pa/(m^3) and s^-1 respectively)
    P0 : double
        Ambient pressure of the recharge reservoir  (Pa)
    M0 : double
        Reservoir mass (kg)
    T0 : double
        Ambient Temperature of the recharge reservoir (deg C)
    bT : double
        Lumped Parameter for the Temperature differential equation (s^-1)
    outputs:
    --------
    dPdt : double
        Derivative of Pressure w.r.t time (Pa/s)
    dTdt : double
        Derivative of Temperature w.r.t time (deg C/s)
    """

    # check if q_stm and q_out are functions
    if not callable(q_stm):
        q_stm_function = lambda t: q_stm # if q_stm is a constant make it a function
    else: 
        q_stm_function = q_stm # otherwise leave it as a function

    if not callable(q_out):
        q_out_function = lambda t: q_out # if q_out is a constant make it a function
    else: 
        q_out_function = q_out # otherwise leave it as a function

    P, T = X # unpack array of pressure and temperature

    dPdt = -aP * (q_out_function(t) - q_stm_function(t)) - bP * (P - P0) # pressure differential equation

    Tprime = T if (P>P0) else T0 # changes depending on whether the Pressure is causing fluid to flow in or out of the system

    dTdt = (q_stm_function(t) / M0) * (Tstm - T) - (bP / (aP * M0)) * (P - P0) * (Tprime - T) - bT * (T - T0) # temperature differential equaiton

    dXdt = [dPdt, dTdt]
    return dXdt


def ode_solve(model, t0, t1, Pi, Ti, pars, time_eval = None):
    """
    Finds the solution to the model differential equation
    
    inputs:
    -------
    model : callable
        function giving the derivative at some point
    t0 : float
        initial time
    t1 : float
        final time 
    Pi : float
        Pressure at time = t0
    Ti : float
        Temperature at time = t0
    pars : array
        array containing all the free parameters in the model
    t_eval : array (optional)
        specific times where we want the solution to be found
    outputs:
    --------
    t : array
        times which the solution was evaluated
    P : array
        Solved Pressures at times t
    T : array
        solved temperatures at times t
    """
    Xi = np.array([Pi, Ti]) # put the initial condition into an array

    pars = tuple(pars) # the solve_ivp function needs the pars and times to be in a tuple
    time = tuple([t0,t1])

    # solve the model ode
    # (dense_output = True allows the sol function to be used which lets you evaluate
    # the solution at any time)
    ode = solve_ivp(ode_model, time, Xi, args=pars, dense_output=True, t_eval=None)
    
    t = ode.t # array of times
    X = ode.y # array containing both the temperature and pressure solutions

    P = X[0,:] # pressure solution
    T = X[1,:] # temperature solution

    return t, P, T

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
        pre_analytical[i] = -(1-np.exp(-i/100))

    for i in range(0, 1000):
        t_analytical[i] = i/100
        tmp_analytical[i]= (-1-0.8*(np.exp(-i)))/(-1-0.1*(np.exp(-i)))

    time, pressure, temp = ode_solve(ode_model, 0, 10 ,0.1 ,0 ,[-1,1,1,1,1,1,1,1,1])

    # ode_solve(model, t0, t1, Pi, Ti, pars, time_eval = None):
    # ode_model(t, X, q_stm, q_out, Tstm, aP, bP, P0, M0, T0, bT)

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

def solve_ode(f, t0, t1, dt, x0, pars, q=None):
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
