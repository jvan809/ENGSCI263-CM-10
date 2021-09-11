from re import L
import matplotlib.pyplot as plt
from data.interpolate_data import *


def plot_given_TP_data():
    #       creates a plot for the given data in temp and pressure
    # Inputs: None
    # Outputs:
    #   shows a plot to the screen with the given temp and pressure data

    # q_oil, q_stm, q_water, tt0, X0 = interpolate_data() # get the interpolated data
    filenames = ["tr_p", "tr_T"] # list of files with the needed data
    filenames = ["." + os.sep + "data" + os.sep + f + ".txt" for f in filenames] # parse these into the correct format

    # get the data from each file
    tP, P_data = get_data(filenames[0])
    tT, T_data = get_data(filenames[1])
    
    # Create structure of plots
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx() # give both plots the same x scale

    # plot data
    plt1 = ax1.plot(tP, P_data, "k.", label = "Pressure Data")
    plt2 = ax2.plot(tT, T_data, "r.", label = "Temperature Data")

    # add labels / legend
    plts = plt1 + plt2
    labs = [l.get_label() for l in plts]
    ax1.legend(plts, labs)
    ax1.set_ylabel("Pressure (kPa)")
    ax2.set_ylabel("Temperature ($^\circ$C)")
    ax1.set_xlabel("time (days)")

    plt.show()


def plot_given_q_data():
    #   creates a plot for the given flow data
    # Inputs: none
    # Outputs:
    #   shows a plot to the screen
    
     # q_oil, q_stm, q_water, tt0, X0 = interpolate_data() # get the interpolated data
    filenames = ["tr_oil", "tr_steam",  "tr_water"] # list of files with the needed data
    filenames = ["." + os.sep + "data" + os.sep + f + ".txt" for f in filenames] # parse these into the correct format

    # get the data from each file
    to, qo = get_data(filenames[0])
    
    ts, qs = get_data(filenames[1])
    
    tw, qw = get_data(filenames[2])

    # create structure of plots
    fig, ax1 = plt.subplots(3)

    qs = qs*2.4527886845*1000 # from steam calculator for specific volume - tonnes * 2.45 m^3/kg * 1000 kg / tonne
    # assumes constant input pressure (Should this be the same as system pressure??)

    # add data
    plt1 = ax1[0].plot(to, qo, "go", label = "Oil Outflow rate")
    plt2 = ax1[1].plot(tw, qw, "co", label = "Water Outflow rate")
    plt3 = ax1[2].plot(ts, qs, "mo", label = "Steam Injection rate")

    # add labels / cosmetic things
    plts = plt1 + plt2 + plt3
    labs = [l.get_label() for l in plts]
    ax1[0].legend(plts, labs)
    [ax1[i].set_ylabel("Flow rate ($m^3$/day)") for i in range(3)]
    [ax1[i].set_xlabel("time (days)") for i in range(3)]
    [ax1[i].set_xlim([0,225]) for i in range(3)] # set x range to be the same for all 3 plots so that they line up

    plt.show()

if __name__ == "__main__":
    plot_given_TP_data()
    plot_given_q_data() 