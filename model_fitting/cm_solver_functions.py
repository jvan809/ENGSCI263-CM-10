import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from scipy.integrate import solve_ivp

def model(t, X, q_stm, q_out, Tstm, aP, bP, P0, M0, T0, bT):
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

    dTdt = (q_stm_function(t)) /(M0) * (Tstm - T) - (bP / (aP * M0)) * (P - P0) * (Tprime - T) - bT * (T - T0) # temperature differential equaiton

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
    ode = solve_ivp(model, time, Xi, args=pars, t_eval=time_eval)
    
    t = ode.t # array of times
    X = ode.y # array containing both the temperature and pressure solutions

    P = X[0,:] # pressure solution
    T = X[1,:] # temperature solution

    

    return t, P, T


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
        k1 = f(t[i], x[i], q[i], *pars)             # /\
        k2 = f(t[i+1], x[i] + k1*dt, q[i+1], *pars) # equations for the improved Euler method
        x[i+1] = x[i] + dt*(k1+k2)/2                # \/
    return t,x


if __name__ == "__main__":
    print(model([1,1],1 ,1,1,1,1,1,1,1,1,1))