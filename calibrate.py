from model_fitting.cm_solver_functions import *
from data.interpolate_data import *
from scipy.optimize import curve_fit


def fit_model(Tstm):
    """
    Finds an optimal set of parameters 

    Tstm : float
        Temperature of injected steam for the model being calibrated
    

    outputs:
    --------
    p : array
        array with the optimal set of parameters to fit the data
    """

    q_oil, q_stm, q_water, tt0, X0 = interpolate_data() # get the interpolated data

    q_out = lambda t: q_oil(t) + q_water(t) # add the oil and water flow functions

    ti = tt0[0]  #|
    Pi = X0[0,0] #| Initial conditions of equation from data
    Ti = X0[1,0] #|
    
    # Function takes in pressure parameters and evaluates the pressure at time t
    pressure_solution = lambda t, aP, bP, P0: solve_and_eval(ti, t, Pi, Ti, q_stm, q_out, Tstm, aP, bP, P0, 1, 1, 1)[0]

    # initial guess for the parameters (see elsewhere for details on decision)
    aP_g, bP_g, P0_g, M0_g, T0_g, bT_g = (1.21433462e-01, 3.06348597e-02 ,6.63684303e+02 ,5.07597613e+3 ,135, 4.02863923e-2)

    # Find optimal set of parameters
    pressure_pars,pressure_cov = curve_fit(pressure_solution, tt0, X0[0,:], p0 = [aP_g, bP_g, P0_g])

    # function takes in temperature parameters and evaluates the solution at time t, using the optimal set of pressure parameters calculated earlier
    temperature_solution = lambda t, M0, T0, bT: solve_and_eval(ti, t, Pi, Ti, q_stm, q_out, Tstm, *pressure_pars, M0, T0, bT)[1]

    # find the optimal set of temperature parameters
    temp_pars, temp_cov = curve_fit(temperature_solution, tt0, X0[1,:], p0 = [M0_g, T0_g, bT_g])

    # append the parameters together into one list
    p = np.append(pressure_pars,temp_pars)
    
    return p, pressure_cov, temp_cov
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
        array containing Pressures and the temperatures appended and evaluated at the times    
    '''

    #for time in t:
    tt, P, T = ode_solve(model, t0, t[-1], Pi, Ti, [q_stm, q_out, Tstm, aP, bP, P0, M0, T0, bT], time_eval=t)
    X = np.array([P.tolist(), T.tolist()])
    return X



if __name__ == "__main__":
    p, pcov, tcov  = fit_model(260)
    print(p)
    print(pcov)
    print(tcov)


    # p = [1.21156628e-01, 3.06502388e-02, 6.61704220e+02, 5.20398876e+03, 1.44063625e+02, 4.41090440e-02]
    

    q_oil, q_stm, q_water, tt0, X0 = interpolate_data() # get the interpolated data
    q_out = lambda t: q_oil(t) + q_water(t) # add the oil and water flow functions

    model_ = lambda t, X, aP, bP, P0, M0, T0, bT: model(t, X, q_stm, q_out, 260, aP, bP, P0, M0, T0, bT)

    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    
    plt1 = ax1.plot(tt0, X0[0], "k.", label = "Pressure Data")
    plt2 = ax2.plot(tt0, X0[1], "r.", label = "Temperature Data")
    tt, P, T = ode_solve(model_, tt0[0], 400, X0[0,0], X0[1,0], p)
    plt3 = ax1.plot(tt, P, "k", label = "Pressure Numerical Sol")    
    plt4 = ax2.plot(tt, T, "r", label = "Temperature Numerical Sol")
   # ax1.plot(tt, q_out(tt), "g")
    plts = plt1 + plt2 + plt3 + plt4 
    labs = [l.get_label() for l in plts]
    ax1.legend(plts, labs)
    plt.show()
    