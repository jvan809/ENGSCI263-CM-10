# ENGSCI263-CM-10

#calibrate.py

  Countains the functions:
1.    		misfit
	       		Calculates the misfit function of a given model
1.    		fit_model
          		Finds an optimal set of parameters
1.    		fit_model_initialcondit_as_pars
          		Finds an optimal set of parameters using the initial condition as variable parameters as well
1.    		solve_and_eval
          		Solves the differential equation and evalutates it at time = t

#main.py
1.	Countians a main function

#plot_data.py
	Contians the functions:
1.		plot_given_TP_data
			creates a plot for the given data in temp and pressure
1.		plot_given_q_data
			creates a plot for the given flow data
#sensitivityAnalysis.py
	Contains the function:
1.		sensitivityAnalysis
			Demonstrates the improved_euler() function defined in cm_solver_functions.py has instability when the step size is high.


#unit_test.py
	Contains the functions:
1.		dummyModel
			 A dummy model simplified enough to obtain an analytic solution for testing purposes. This model is based on the ODE for pressure.
1.		solve_ode_test
			Tests: solve_ode(model, t0, t1, Pi, Ti, pars, time_eval = None) From : cm_solver_functions.py
1.		model_test
			Tests: model(t, X, q_stm, q_out, Tstm, aP, bP, P0, M0, T0, bT) From : cm_solver_functions.py
