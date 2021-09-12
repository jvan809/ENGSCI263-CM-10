import numpy as np
from operator import xor
from scipy.interpolate import interp1d


def const_flow(time, injtime = 60, tottime = 150, offset = 0, vol = 1000, isProduction = 0):
    # gives injection / output volume at specified time 
    #
    # Inputs:
    #   time: float - time at which flow is evaluated (must be consistent with inj and total time)
    #   injtime: float - amount of time (same units as time) that steam is injected
    #   tottime: float - length (in same units) of full injection + production cycle (i.e. amount of time between starting of one injection and the next)
    #   vol: float - constant value to be applied when in appropirate phase
    #   isProduction: bool - 0 = first phase has non-zero, 1 = second phase has non-zero
    #
    # Outputs:
    #   flow: amount of steam injected / output at time

    # find where each time is in the full cycle
    x = (time - offset) % tottime

    flow = vol*xor(x <= injtime, isProduction)

    return flow

def interp_flow(time, tottime = 150, offset = 0, vols = [1000,1000,0,0], times = [0,60,60.1,150]):
    # 
    # Inputs:
    #   time: time to evaluate injection volume
    #   tottime: float - length of the full cycle, after which behaviour repeats (same units as time, otherwise not important)
    #   vols: list - y values for points to interpolate
    #   times: list same length as vols - x values for same
    #
    # Outputs:
    #   flow - amount of flow at the given time 
    
    x = (time - offset) % tottime

    f = interp1d(times, vols, kind = 'linear')

    return f(x)


if __name__ == "__main__":
    # unit tests
    assert const_flow(30) == 1000
    assert const_flow(90) == 0
    assert const_flow(90, isProduction=1) == 1000
    assert const_flow(30, isProduction=1) == 0

    assert interp_flow(0) == 1000
    assert interp_flow(60) == 1000
    assert interp_flow(90) == 0