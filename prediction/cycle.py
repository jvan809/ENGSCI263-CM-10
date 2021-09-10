import numpy as np
from operator import xor

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
    #   flow: amount of steam injected / outpu t at time

    # find where each time is in the full cycle
    x = (time - offset) % tottime

    flow = vol*xor(x <= injtime, isProduction)

    return flow

def linear_vol(time, injtime = 60, tottime = 150, injstart = 1000, injend = 800):
    # 
    # Inputs:
    #   time: time to evaluate injection volume
    # 
    #   injstart: Injection volume at the start of the injection phase
    #   injend: Injection volume at end of injection phase (just before production) 
    
    x = time % tottime

    if x > injtime:
        return 0
    
    injection = (1 - x / injtime)*(injstart - injend) + injend

    return injection


if __name__ == "__main__":
    pass