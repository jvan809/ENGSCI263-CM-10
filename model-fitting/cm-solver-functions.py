import numpy as np
import matplotlib.pyplot as plt
import scipy as sp


def model(P, T, q_stm, q_out, aP, bP, P0, M0, Tstm, T0, bT):
    """
    The model giving the derivaitves for the differential equation
    
    inputs:
    -------
    P : double
        Pressure of the system (Pa)
    T : double
        Temperature of the system  (deg C)
    q_stm : double
        Flow of steam into the system (m^3/s)
    q_out : double
        Flow of water and bitumen out of the system (m^3/s)
    aP, bP : double
        Lumped Parameters for Pressure Equation (units Pa/(m^3) and s^-1 respectively)
    P0 : double
        Ambient pressure of the recharge reservoir  (Pa)
    M0 : double
        Reservoir mass (kg)
    Tstm : double
        Injected Steam Temperature (deg C)
    T0 : double
        Ambient Temperature of the recharge reservoir (deg C)
    bT : double
        Lumped Parameter for the Temperature differential equation (s^-1)

    outputs:
    --------
    dPdt : double
        Derivative of Pressure w.r.t time (Pa/s)
    dTdt : double
        Derivative of Temperature w.r.t time (deg C/s)
    """

    dPdt = -aP * (q_out - q_stm) - bP * (P - P0)

    Tprime = T if (P>P0) else T0

    dTdt = (q_stm / M0) * (Tstm - T) - (bP / (aP * M0)) * (P - P0) * (Tprime - T) - bT *(T - T0)

    return dPdt, dTdt