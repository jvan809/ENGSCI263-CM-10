from prediction.uncertanity import *
from calibrate import *
from prediction.cycle import *
from model_fitting.cm_solver_functions import *
from data.interpolate_data import *
from matplotlib.pyplot import *
from plot_data import *
from sensitivityAnalysis import *
# this file should be able to produce all figures included in reports **without modification**
# Graphs Used in template report:
#   production & pressure - could be production / injection + temp?
#   best fit of first version of model + misfit over time
#   as before with improved model (+ inital allowed to vary)
#   prediction with varied scenarios (unchanged, constant max, similar to before scaled up??)
#   as before with ensemble modelling
#   posterior dist of one parameter (possibly not needed for this application)
#   - I'd like to add dist of max temp for max all the time scenario





if __name__ == "__main__":
    # plot the pilot data given to us
    plot_given_TP_data() # temperature and pressure data
    plot_given_q_data() # inflow and outflow data
    
    
    pp = [1.20508397e-01, 3.09508220e-02, 6.56020856e+02, 5.08928148e+03, 1.45644569e+02, 4.77635973e-02]
    Pi = 1.44320757e+03
    Ti = 1.92550478e+02
    pars = [1.44320757e+03, 1.20508397e-01, 3.09508220e-02, 6.56020856e+02, 1.92550478e+02, 5.08928148e+03, 1.45644569e+02, 4.77635973e-02]
    
    q_oil, q_stm, q_water, tt0, X0 = interpolate_data() # get the interpolated data

    endTime = tt0[-1]

    q_out = lambda t: q_oil(t) + q_water(t) # add the oil and water flow functions

    # interpolated model
    model_ = lambda t, X, aP, bP, P0, M0, T0, bT: model(t, X, q_stm, q_out, 260, aP, bP, P0, M0, T0, bT)


    # constant maximum in/out flow
    stm1 = lambda t: const_flow(t, 60, 150, endTime, 1000, 0) # 2450000
    out1 = lambda t: const_flow(t, 60, 150, endTime, 250, 1)
    models = [lambda t, X, aP, bP, P0, M0, T0, bT: model(t, X, stm1, out1, 260, aP, bP, P0, M0, T0, bT)]
    names = ["1000 tonnes / day constant"]

    # less conservative linear model
    stm2 = lambda t: interp_flow(t, 150, vols = [1000,800,0,0], times = [0,60,60.1,150], offset=endTime) # first cycle * 1000/460 approx - 2500000,2100000
    out2 = lambda t: interp_flow(t, 150, vols = [0,0,250,200], times = [0,60,60.1,150], offset=endTime) # guess
    models.append( lambda t, X, aP, bP, P0, M0, T0, bT: model(t, X, stm2, out2, 260, aP, bP, P0, M0, T0, bT) )
    names.append("Extension of pilot injection")

    # constant 900 tonnes / day
    stm4 = lambda t: const_flow(t, 60, 150, endTime, 900, 0) # 2200000
    out4 = lambda t: const_flow(t, 60, 150, endTime, 225, 1)
    models.append( lambda t, X, aP, bP, P0, M0, T0, bT: model(t, X, stm4, out4, 260, aP, bP, P0, M0, T0, bT) )
    names.append("900 tonnes/day constant")

    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()

    tt, P, T = ode_solve(model_, tt0[0], tt0[-1], Pi, Ti, pp)

    # interpolation
    plts = ax1.plot(tt0, X0[0], "k.", label = "Pressure Data")
    plts += ax2.plot(tt0, X0[1], "r.", label = "Temperature Data")
    plts += ax1.plot(tt, P, "k", label = "Pressure Numerical Sol")    
    plts += ax2.plot(tt, T, "r", label = "Temperature Numerical Sol")

    # predictions start from where interpolation left off
    Pie = P[-1]
    Tie = T[-1]

    finalTime = 650
    # predictions to plot
    for i in range(3):
        m = models[i]
        col = ['c','g','y'][i]
        tt, P, T = ode_solve(m, tt0[-1], finalTime, Pie, Tie, pp)
        plts += ax2.plot(tt, T, col, label = "Temperature Prediction " + names[i])


    #plt.show()




    pcov = [[ 7.85744566e-05, -1.08497587e-08, -4.06830883e-01],
             [-1.08497587e-08,  1.19555937e-09, -5.37910126e-04],
             [-4.06830883e-01, -5.37910126e-04 , 2.40294395e+03]]
   
   
    tcov = [[ 4.02138041e+03, -1.75419958e+01,  1.08953251e-01],
             [-1.75419958e+01,  5.30715709e+00,  1.05492717e-03],
             [ 1.08953251e-01,  1.05492717e-03 , 3.40814200e-06]]

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
    fig, ax1 = plt.subplots()
    for i in range(3):
        plt.hist(maxTemps[i], density = True, label = names[i])
        
    plt.legend(names)
    plt.title("Max Temp reached")
    plt.show()


    plots = [0,0,0]
    fig, (ax0,ax1,ax2) = plt.subplots(3)
    colours = ["b", "tab:orange", "g"]
    
    plot0 = ax0.hist(maxTemps[0], density = True, label = names[0], color = colours[0])[2]
    plot1 = ax1.hist(maxTemps[1], density = True, label = names[1], color = colours[1])[2]
    plot2 = ax2.hist(maxTemps[2], density = True, label = names[2], color = colours[2])[2]

        
    plots = plot0 + plot1 + plot2
    labs = [l.get_label() for l in plots]
    ax0.legend(names)
    ax0.legend(plots, labs)
    plt.show()

    # Calling convergence analysis function (this function plots the convegence analysis)
    sensitivityAnalysis()


