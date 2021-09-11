
from model_fitting.cm_solver_functions import *

def dummyModel(t, x, q, a, b, x0):
    """
    A dummy model simplified enough to obtain an analytic solution
    for testing purposes. This model is based on the ODE for pressure.
    """
    
    # pressure differential equation:
    dPdt = a*q - b*(x-x0)

    return dPdt

def solve_ode_test():
    """
    Tests: solve_ode(model, t0, t1, Pi, Ti, pars, time_eval = None)
    From : cm_solver_functions.py
    """
    # Note: using the improved euler's method, not the Scipy function.
    [t, x] = improved_euler(dummyModel, t0 = 0, t1 = 3, dt = 1, x0 = 50, pars = [1,1,1])
    assert (x == [50.0, 25.0, 12.50]).all(), "function solve_ode producting incorrect output (improved Euler's method)"

def model_test():
    """
    Tests: model(t, X, q_stm, q_out, Tstm, aP, bP, P0, M0, T0, bT)
    From : cm_solver_functions.py
    """
    [dPdt, dTdt] = model(t=1, X=[1,1], q_stm=1, q_out=1, Tstm=1, aP=1, bP=1, P0=2, M0=1, T0=3, bT=1)
    assert [dPdt, dTdt] == [1, 4], "model function producing incorrect ouput"


if __name__ == "__main__":
    solve_ode_test()
    model_test()