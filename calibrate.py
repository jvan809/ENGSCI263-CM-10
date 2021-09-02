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
    Ti = X0[1,0]
    t0, index = np.meshgrid(t0, [0,1])
    def free_model(xtuple, aP, bP, P0, M0, T0, bT): 
        (t,index) = xtuple
        return solve_and_eval(ti, t, Pi, Ti, q_stm, q_out, Tstm, aP, bP, P0, M0, T0, bT)

    # choose initial guesses for the parameters
    aP_g = 1
    bP_g = 1
    P0_g = 1
    M0_g = 1
    T0_g = 1
    bT_g = 1
    
    p,_ = curve_fit(free_model, (t0, index), X0.ravel(), p0 = [aP_g, bP_g, P0_g, M0_g, T0_g, bT_g])

    return p

def solve_and_eval(t0, t, Pi, Ti, q_stm, q_out, Tstm, aP, bP, P0, M0, T0, bT):
    '''
    Solves the differential equation and evalutates it at time = t

    t0 : float
        time at which the initial pressure and temperature is evaluated
    t : float/array
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
    index : double

    outputs:
    --------
    X : array
        array containing either the Pressure or the temperature evaluated at time t    
    '''

    #for time in t:
    tt, P, T = ode_solve(model, t0, t[-1,-1], Pi, Ti, [q_stm, q_out, Tstm, aP, bP, P0, M0, T0, bT], time_eval=t[0,:])
    return np.append(P,T)



if __name__ == "__main__":
    p = fit_model(260)
    print(p)