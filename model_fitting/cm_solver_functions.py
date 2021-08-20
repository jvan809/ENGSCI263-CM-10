import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from scipy.integrate import solve_ivp

def model(t, X, q_stm, q_out, aP, bP, P0, M0, Tstm, T0, bT):
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
    q_stm : callable
        Flow of steam into the system (m^3/s)
    q_out : callable
        Flow of water and bitumen out of the system (m^3/s)
    aP, bP : double
        Lumped Parameters for Pressure Equation (units Pa/(m^3) and s^-1 respectively)
    P0 : double
        Ambient pressure of the recharge reservoir  (Pa)
    M0 : double
        Reservoir mass (kg)
    Tstm : double
        Injected Steam Temperature (deg C)
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

    P, T = X # unpack array of pressure and temperature

    dPdt = aP*(q_out - q_stm) - bP * (P - P0) # pressure differential equation

    Tprime = T if (P>P0) else T0 # changes depending on whether the    

    dTdt = (q_stm / M0) * (Tstm - T) - (bP / (aP * M0)) * (P - P0) * (Tprime - T) - bT *(T - T0)

    dXdt = [dPdt, dTdt]
    return dXdt


def solve_ode_scipy(model, interval, Pi, Ti, pars, t_eval = None):

    Xi = np.array([Pi, Ti])
    arg = tuple(pars)
    time = tuple(interval)
    ode = solve_ivp(model, time, Xi, dense_output=True)



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


if __name__ == "__main__":
    print(model([1,1],1 ,1,1,1,1,1,1,1,1,1))