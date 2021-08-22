from model_fitting.cm_solver_functions import *
from data.interpolate_data import *
from scipy.optimize import curve_fit


def fit_model(Tstm):
    """
    
    """

    q_oil, q_stm, q_water, t0, X0 = interpolate_data()

    q_out = lambda t: q_oil(t) + q_water(t)

    ti = t0[0]
    Pi = X0[0,0]


    free_model = lambda t, X, aP, bP, P0, M0, T0, bT: solve_and_eval(t, X, q_stm, q_out, Tstm, aP, bP, P0, M0, T0, bT)

    

    p,_ = curve_fit(free_model, )



def solve_and_eval(t0, t, Pi, Ti, q_stm, q_out, Tstm, aP, bP, P0, M0, T0, bT):
    '''
    Solves the differential equation and evalutates it at time = t

    t0 : float
        time at which the initial pressure and temperature is evaluated
    t : float
        time at which the solution is to be solved
    Pi : double
        Initial Pressure of the system (Pa)
    Ti : double
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
    X : array
        array containing the Pressure and temperature evaluated at time t    
    '''
    tt, P, T, sol = ode_solve(model, t0[0], t, Pi, Ti, q_stm, q_out, Tstm, aP, bP, P0, M0, T0, bT)
    
    X = sol(t)
    return X