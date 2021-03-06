from sensitivityAnalysis import sensitivityAnalysis
from prediction.uncertanity import *
from calibrate import *
from prediction.cycle import *
from model_fitting.benchmark import *
from model_fitting.cm_solver_functions import *
from data.interpolate_data import *
from matplotlib.pyplot import *
from plot_data import *
from sensitivityAnalysis import *
from saleski_pilot import *



if __name__ == "__main__":



    # plot the pilot data given to us
    plot_given_TP_data() # temperature and pressure data
    plot_given_q_data() # inflow and outflow data
    
    # Calling convergence analysis function (this function plots the convegence analysis)
    sensitivityAnalysis()

    # plot the benchmarking
    plot_benchmark()

    # calibrate to the saleski model and check the error
    saleski_model()
    
    
    # optimal set of parameters without varying Pi and Ti
    p = [1.21156628e-01, 3.06502388e-02, 6.61704220e+02, 5.20398876e+03, 1.44063625e+02, 4.41090440e-02]
    
    # with Pi and Ti free
    # parameter values found from calibrate.py - neglected to be ran here due to time to run
    pp = [1.20508397e-01, 3.09508220e-02, 6.56020856e+02, 5.08928148e+03, 1.45644569e+02, 4.77635973e-02]
    Pi = 1.44320757e+03
    Ti = 1.92550478e+02
    pars = [1.44320757e+03, 1.20508397e-01, 3.09508220e-02, 6.56020856e+02, 1.92550478e+02, 5.08928148e+03, 1.45644569e+02, 4.77635973e-02]
    
    q_oil, q_stm, q_water, tt0, X0 = interpolate_data() # get the interpolated data

    endTime = tt0[-1] # used for starting extrapolation

    q_out = lambda t: q_oil(t) + q_water(t) # add the oil and water flow functions

    # interpolated model
    model_ = lambda t, X, aP, bP, P0, M0, T0, bT: model(t, X, q_stm, q_out, 260, aP, bP, P0, M0, T0, bT)


    # initial pass numerical solution
    fig, ax1 = plt.subplots() # create subplots
    #fig, (ax1,ax1_) = plt.subplots(2) # create subplots
    ax2 = ax1.twinx()

    plt1 = ax1.plot(tt0, X0[0], "k.", label = "Pressure Data") # pressure data
    plt2 = ax2.plot(tt0, X0[1], "r.", label = "Temperature Data") # temperature data
    tt, P, T = ode_solve(model_, tt0[0], 400, Pi, Ti, pp, time_eval=np.linspace(tt0[0],400, 1000))
    plt3 = ax1.plot(tt, P, "k", label = "Pressure Numerical Sol") # 
    plt4 = ax2.plot(tt, T, "r", label = "Temperature Numerical Sol")
    plts = plt1 + plt2 + plt3 + plt4 
    ax1.set_xlabel("time (days)")
    ax1.set_ylabel("Pressure (kPa)")
    ax2.set_ylabel("Temperature ($^\circ$C)")
    labs = [l.get_label() for l in plts]
    ax1.legend(plts, labs)
    plt.title("model fitted against the data")
    plt.show()



    
    # constant maximum in/out flow
    stm1 = lambda t: const_flow(t, 60, 150, endTime, 1000, 0)
    out1 = lambda t: const_flow(t, 60, 150, endTime, 250, 1)
    models = [lambda t, X, aP, bP, P0, M0, T0, bT: model(t, X, stm1, out1, 260, aP, bP, P0, M0, T0, bT)]
    names = ["1000 tonnes / day constant"]

    # less conservative linear model
    stm2 = lambda t: interp_flow(t, 150, vols = [1000,850,0,0], times = [0,60,60.1,150], offset=endTime) # first cycle * 1000/460 
    out2 = lambda t: interp_flow(t, 150, vols = [0,0,250,200], times = [0,60,60.1,150], offset=endTime) # guess
    models.append( lambda t, X, aP, bP, P0, M0, T0, bT: model(t, X, stm2, out2, 260, aP, bP, P0, M0, T0, bT) )
    names.append("Extension of pilot injection")

    # constant 900 tonnes / day
    stm4 = lambda t: const_flow(t, 60, 150, endTime, 900, 0) # 2200000
    out4 = lambda t: const_flow(t, 60, 150, endTime, 225, 1)
    models.append( lambda t, X, aP, bP, P0, M0, T0, bT: model(t, X, stm4, out4, 260, aP, bP, P0, M0, T0, bT) )
    names.append("900 tonnes/day constant")


    # # misfit of first model
    # fig = figure()
    # fig.subplots_adjust()
    # ax1 = fig.add_subplot()
    # ax2 = ax1.twinx()

    # ax2.set_ylabel("Temperature (Deg C)")
    # ax1.set_ylabel("Pressure (Pa)")
    # ax1.set_xlabel("Time (Days)")

    # tt, P, T = ode_solve(model_, tt0[0], 400, Pi, Ti, pp, time_eval=tt0)

    # misP = P - X0[0]
    # misT = T - X0[1]

    # plts = ax1.plot(tt0, misP, "kx", label = "Pressure Misfit")
    # plts += ax2.plot(tt0, misT, "rx", label = "Temperature Misfit")
   
    # labs = [l.get_label() for l in plts]
    # ax1.legend(plts, labs)
    # plt.title("Misfit of Temp and Pressure, First model")
    # plt.show()


    # uniform error without the initial condition as a free parm
    uniform_error(p, title = "Uniform error of Temp and Pressure, First model")

    # plots the 2nd pass model - wrapped up in function to make for less repetition
    def plotBestFit():
        fig = figure()
        fig.subplots_adjust()
        ax1 = fig.add_subplot()
        ax2 = ax1.twinx()

        ax2.set_ylabel("Temperature (Deg C)")
        ax1.set_ylabel("Pressure (Pa)")
        ax1.set_xlabel("Time (Days)")

        # initial values based on params
        tt, P, T = ode_solve(model_, tt0[0], tt0[-1], Pi, Ti, pp, tt0)

        # interpolation
        plts = ax1.plot(tt0, X0[0], "k.", label = "Pressure Data")
        plts += ax2.plot(tt0, X0[1], "r.", label = "Temperature Data")
        plts += ax1.plot(tt, P, "k", label = "Pressure Numerical Sol")    
        plts += ax2.plot(tt, T, "r", label = "Temperature Numerical Sol")

        
        return P, T, plts, ax2


    P, T, plts, ax2 =  plotBestFit()
    labs = [l.get_label() for l in plts]
    ax2.legend(plts, labs)
    plt.title("With Initial Value Parameters")
    plt.show()

    # # new misfit for 2nd pass
    # fig = figure()
    # fig.subplots_adjust()
    # ax1 = fig.add_subplot()
    # ax2 = ax1.twinx()

    # ax2.set_ylabel("Temperature (Deg C)")
    # ax1.set_ylabel("Pressure (Pa)")
    # ax1.set_xlabel("Time (Days)")

    # misP = P - X0[0]
    # misT = T - X0[1]

    # plts = ax1.plot(tt0, misP, "kx", label = "Pressure Misfit")
    # plts += ax2.plot(tt0, misT, "rx", label = "Temperature Misfit")
   
    # labs = [l.get_label() for l in plts]
    # ax1.legend(plts, labs)
    # plt.title("Misfit of Temp and Pressure, with IV params")
    # plt.show()
    
    
    # find the uniform error between the model and the given data
    uniform_error(pars, init_is_pars=True, title = "Uniform error of Temp and Pressure, with IV params")

    # add predictions for scenarios
    def plotPred(P, T, plts, ax2):
        # predictions start from where interpolation left off
        Pie = P[-1]
        Tie = T[-1]

        # predictions to plot
        for i in range(3):
            m = models[i]
            col = ['c','g','y'][i]
            tt, P, T = ode_solve(m, tt0[-1], finalTime, Pie, Tie, pp)
            plts += ax2.plot(tt, T, col, label = "Temperature Prediction " + names[i])


    
    finalTime = 650

    P, T, plts, ax2 =plotBestFit()
    plotPred(P, T, plts, ax2)
    plt.title("Prediction")
    labs = [l.get_label() for l in plts]
    ax2.legend(plts, labs)
    plt.show()

    # adding ensemble of models
    P, T, plts, ax2 =plotBestFit()
    plotPred(P, T, plts, ax2)

    labs = [l.get_label() for l in plts]
    ax2.legend(plts, labs)
    plt.title("Prediction with Ensemble of Models")

    # covariance matricies also generated from calibrate.py - also neglected here due to time to run
    pcov = [[ 7.85744566e-05, -1.08497587e-08, -4.06830883e-01],
             [-1.08497587e-08,  1.19555937e-09, -5.37910126e-04],
             [-4.06830883e-01, -5.37910126e-04 , 2.40294395e+03]]
   
   
    tcov = [[ 4.02138041e+03, -1.75419958e+01,  1.08953251e-01],
             [-1.75419958e+01,  5.30715709e+00,  1.05492717e-03],
             [ 1.08953251e-01,  1.05492717e-03 , 3.40814200e-06]]

    # note: fixed seed in uncertanity.py
    normModels = normSample(pp, pcov, tcov, 100)

    maxTemps = [[],[],[]]

    for params in normModels:
        _, p, t = ode_solve(model_, tt0[0], tt0[-1], Pi, Ti, params, tt0)
        ax1.plot(tt0,p,'k-', lw=0.25,alpha=0.2)
        ax2.plot(tt0,t,'r-', lw=0.25,alpha=0.2)
        Pie = p[-1]
        Tie = t[-1]
        for i in range(3):
            col = ['c','g','y'][i]
            tt, p, t = ode_solve(models[i], tt0[-1], finalTime, Pie, Tie, params)
            ax2.plot(tt,t,col, lw=0.25,alpha=0.2)
            maxTemps[i].append(max(t))



    labs = [l.get_label() for l in plts]
    ax1.legend(plts, labs)


    plt.show()

    # plots the pdf histogram

    fig, (ax0,ax1,ax2) = plt.subplots(3,sharex=True)
    colours = ["b", "tab:orange", "g"]
    
    plot0 = ax0.hist(maxTemps[0], density = True, label = names[0], color = colours[0])
    plot1 = ax1.hist(maxTemps[1], density = True, label = names[1], color = colours[1])
    plot2 = ax2.hist(maxTemps[2], density = True, label = names[2], color = colours[2])

    ax0.vlines(np.percentile(maxTemps[0], 5),0,0.5, color = "r", linestyles = "-.")
    ax0.vlines(np.percentile(maxTemps[0], 95),0,0.5, color = "r", linestyles = "-.")
    ax1.vlines(np.percentile(maxTemps[1], 5),0,0.5, color = "r", linestyles = "-.")
    ax1.vlines(np.percentile(maxTemps[1], 95),0,0.5, color = "r", linestyles = "-.")
    ax2.vlines(np.percentile(maxTemps[2], 5),0,0.5, color = "r", linestyles = "-.")
    ax2.vlines(np.percentile(maxTemps[2], 95),0,0.5, color = "r", linestyles = "-.")

    
    for i in range(3):
        print(names[i] + " 95% CI maximum temperature: ", round(np.percentile(maxTemps[i], 95), 1))
        
    plots = plot0 + plot1 + plot2
    ax0.set_title(names[0])
    ax1.set_title(names[1])
    ax2.set_title(names[2])
    ax2.set_xlabel("max temperature $^\circ$C")
    plt.show()




