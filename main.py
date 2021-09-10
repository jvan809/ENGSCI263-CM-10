from prediction.uncertanity import *
from calibrate import *
from prediction.cycle import *
from model_fitting.cm_solver_functions import *
from data.interpolate_data import *
from scipy.optimize import curve_fit









if __name__ == "__main__":
    pp = [1.20508397e-01, 3.09508220e-02, 6.56020856e+02, 5.08928148e+03, 1.45644569e+02, 4.77635973e-02]
    Pi = 1.44320757e+03
    Ti = 1.92550478e+02
    pars = [1.44320757e+03, 1.20508397e-01, 3.09508220e-02, 6.56020856e+02, 1.92550478e+02, 5.08928148e+03, 1.45644569e+02, 4.77635973e-02]
    
    q_oil, q_stmo, q_water, tt0, X0 = interpolate_data() # get the interpolated data

    q_outo = lambda t: q_oil(t) + q_water(t) # add the oil and water flow functions
    q_stm = lambda t: q_stmo(t) if (t < tt0[-1]) else const_flow(t, offset = tt0[-1])
    q_out = lambda t: q_outo(t) if (t < tt0[-1]) else const_flow(t,vol = 250 ,isProduction = 1, offset = tt0[-1]) # note: total guess

    model_ = lambda t, X, aP, bP, P0, M0, T0, bT: model(t, X, q_stm, q_out, 260, aP, bP, P0, M0, T0, bT)

    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()

    tt, P, T = ode_solve(model_, tt0[0], 400, Pi, Ti, pp)


    plt1 = ax1.plot(tt0, X0[0], "k.", label = "Pressure Data")
    plt2 = ax2.plot(tt0, X0[1], "r.", label = "Temperature Data")
    plt3 = ax1.plot(tt, P, "k", label = "Pressure Numerical Sol")    
    plt4 = ax2.plot(tt, T, "r", label = "Temperature Numerical Sol")
    plts = plt1 + plt2 + plt3 + plt4 
    labs = [l.get_label() for l in plts]
    ax1.legend(plts, labs)


    pcov = [[ 7.85744566e-05, -1.08497587e-08, -4.06830883e-01],
             [-1.08497587e-08,  1.19555937e-09, -5.37910126e-04],
             [-4.06830883e-01, -5.37910126e-04 , 2.40294395e+03]]
   
   
    tcov = [[ 4.02138041e+03, -1.75419958e+01,  1.08953251e-01],
             [-1.75419958e+01,  5.30715709e+00,  1.05492717e-03],
             [ 1.08953251e-01,  1.05492717e-03 , 3.40814200e-06]]

    normModels = normSample(pp, pcov, tcov, 100)


    for params in normModels:
        _, p, t = ode_solve(model_, tt0[0], 400, Pi, Ti, params, tt)
        ax1.plot(tt,p,'k-', lw=0.25,alpha=0.2)
        ax2.plot(tt,t,'r-', lw=0.25,alpha=0.2)
    
    ax1.plot([],[],'k-', lw=0.5,alpha=0.4, label='model ensemble (P)')
    ax2.plot([],[],'r-', lw=0.5,alpha=0.4, label='model ensemble (T)')



    plt.show()
