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
			
#cycle.py
	Contains the functions:
1.		const_flow
			gives injection / output volume at specified time 
1.		interp_flow
			
#uncertanity.py 
	Contains the functions:
1.		findGridBorders

1.		gridSample

1.		normSample

#benchmark.py
	Contains the functions:
1.		plot_benchmark
			Analytical pressure solution for parameters [-1,0,1,1,1,1,1,1,1]
1.		temperature_analytical
			Analytical temperature solution for parameters [-1,0,1,1,1,1,1,1,1]
1.		pressure_analytical
			Compare analytical and numerical solutions.

#cm_solver_functions.py 
	Contains the functions:
1.		model
			The model giving the derivaitves for the differential equation
1.		ode_solve
			Finds the solution to the model differential equation
1.		improved_euler
			Solve an ODE numerically.
			
#interpolate_data.py
	Contains the functions:
1.		get_data
			gets data from a txt file like the ones given for the nz pilot data
1.		interpolate_data
			Interpolates the nz pilot thermal recovery data
