from operator import delitem
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import glob as gl
import os 
from scipy.interpolate import interp1d


def get_data(filename):
    """
    gets data from a txt file like the ones given for the nz pilot data

    inputs:
    -------
    filename : str
        location of the txt file
    
    outputs:
    -------
    t : array
        array of times the data is collected
    data : array
        values taken from the txt file
    """

    A = np.genfromtxt(filename, skip_header = 1, delimiter=",")
    t, data = A.T
    return t, data


def interpolate_data():
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
    filenames = ["tr_oil", "tr_p", "tr_steam", "tr_T", "tr_water"] # list of files with the needed data
    filenames = ["." + os.sep + "data" + os.sep + f + ".txt" for f in filenames] # parse these into the correct format

    # get the data from each file
    to, qo = get_data(filenames[0])
    tP, P_data = get_data(filenames[1])
    ts, qs = get_data(filenames[2])
    tT, T_data = get_data(filenames[3])
    tw, qw = get_data(filenames[4])

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



if __name__ == "__main__":
    interpolate_data()