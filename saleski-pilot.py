from model_fitting.cm_solver_functions import *

from data.interpolate_data import *

from calibrate import *
def interpolate_saleski_data():
    """
    Interpolates the nz pilot thermal recovery data

    inputs:
    -------
    filename : str
        location of the txt file
    
    outputs:
    -------
    q_oil : callable
        interpolated oil flow 
    q_stm : callable
        interpolated steam flow  
    q_wat : callable
        interpolated water flow         
    tP : array
        times at which the temperature and pressure were evaluated
    Xint : array
        Pressures and temperatures evaluated at the same times
    """
    filenames = ["exp2_oil", "exp2_pressure", "exp2_steam", "exp2_temp", "exp2_water"] # list of files with the needed data
    filenames = ["." + os.sep + "data" + os.sep + "pilot_project_data_2012" + os.sep + f + ".txt" for f in filenames] # parse these into the correct format

    # get the data from each file
    to, qo = get_data(filenames[0])
    tP, P_data = get_data(filenames[1])
    ts, qs = get_data(filenames[2])
    tT, T_data = get_data(filenames[3])
    tw, qw = get_data(filenames[4])


    qs = 1000 * qs / 0.4323625642 # convert from tonnes/day to m^3/day - conversion factor taken from steam table

    # fix units for steam
    #qs = qs*2.4527886845*1000 # from steam calculator for specific volume - tonnes * 2.45 m^3/kg * 1000 kg / tonne
    # assumes constant input pressure (Should this be the same as system pressure??)

    # interpolate all the data
    q_oil = interp1d(to, qo, fill_value = 0, bounds_error=False)
    Pr = interp1d(tP, P_data)
    q_stm = interp1d(ts, qs, fill_value = 0, bounds_error=False)
    Te = interp1d(tT, T_data)
    q_wat = interp1d(tw, qw, fill_value = 0, bounds_error=False)
    
    T_inter = Te(tP) # evaluate the temperature at the time when the pressure was evaluated
    Xint = np.array([P_data, T_inter]) # combine the pressure and temperature into one array

    

    return q_oil, q_stm, q_wat, tP, Xint

def fit_model_saleski():
    """
    Finds an optimal set of parameters using the initial condition as variable parameters as well

    Tstm : float
        Temperature of injected steam for the model being calibrated
    

    outputs:
    --------
    p : array
        array with the optimal set of parameters to fit the data
    """

    q_oil, q_stm, q_water, tt0, X0 = interpolate_saleski_data() # get the interpolated data

    q_out = lambda t: q_oil(t) + q_water(t) # add the oil and water flow functions

    ti = tt0[0] 

    Tstm = 230


    
    # Function takes in pressure parameters as well as the pressure initial condition and evaluates the pressure at time t
    pressure_solution = lambda t, Pi, aP, bP, P0: solve_and_eval(ti, t, Pi, X0[1,0], q_stm, q_out, Tstm, aP, bP, P0, 1e11, 1, 0)[0]

    # initial guess for the parameters (see elsewhere for details on decision)
    Pi_g, Ti_g, aP_g, bP_g, P0_g, M0_g, T0_g, bT_g = (X0[0,0], X0[1,0], 7.20508397e-05, 3.09508220e-02, 120, 7.08928148e+06, 130, 1.77635973e-2)

    # Find optimal set of parameters
    pressure_pars,pressure_cov = curve_fit(pressure_solution, tt0[:45], X0[0,:45], p0 = [Pi_g, aP_g, bP_g, P0_g])

    p_pars = pressure_pars[1:]

    # function takes in temperature parameters and evaluates the solution at time t, using the optimal set of pressure parameters calculated earlier
    temperature_solution = lambda t, Ti, M0, T0, bT: solve_and_eval(ti, t, pressure_pars[0], Ti, q_stm, q_out, Tstm, *p_pars, M0, T0, bT)[1]

    # find the optimal set of temperature parameters
    temp_pars, temp_cov = curve_fit(temperature_solution, tt0[:45], X0[1,:45], p0 = [Ti_g, M0_g, T0_g, bT_g])

    # append the parameters together into one list
    p = np.append(pressure_pars,temp_pars)
    
    return p




if __name__ == "__main__":
    p = fit_model_saleski()
    print(p)

    Pi = p[0]
    Ti = p[4]
    p = np.delete(p, [0,4])
    

    q_oil, q_stm, q_water, tt0, X0 = interpolate_saleski_data() # get the interpolated data
    q_out = lambda t: q_oil(t) + q_water(t) # add the oil and water flow functions

   
    # p= [7.20508397e-05, 3.09508220e-02, 120, 7.08928148e+06, 130, 1.77635973e-2] 
    # Pi = X0[0,0]
    # Ti = X0[1,0]

    model_ = lambda t, X, aP, bP, P0, M0, T0, bT: model(t, X, q_stm, q_out, 230, aP, bP, P0, M0, T0, bT)

    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()

    tt, P, T = ode_solve(model_, tt0[0], tt0[-1]+10, Pi, Ti, p,time_eval = tt0)

    dP = (P - X0[0,:])/(1+X0[0,:])
    dT = (T - X0[1,:])/(1+X0[1,:])

    plt1 = ax1.plot(tt0, X0[0], "k.", label = "Pressure Data")
    plt2 = ax2.plot(tt0, X0[1], "r.", label = "Temperature Data")
    plt3 = ax1.plot(tt, P, "k", label = "Pressure Numerical Sol")    
    plt4 = ax2.plot(tt, T, "r", label = "Temperature Numerical Sol")
    plt5 = ax1.vlines(tt0[45],0,2700,"b")
    plts = plt1 + plt2 + plt3 + plt4 
    labs = [l.get_label() for l in plts]
    ax1.legend(plts, labs)
    plt.show()

    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()   

    plt1 = ax1.plot(tt0, dP, "kx", label = "Pressure Data")
    plt2 = ax2.plot(tt0, dT, "rx", label = "Temperature Data")
    plt.show()