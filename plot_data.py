from re import L
import matplotlib.pyplot as plt
from data.interpolate_data import *


def plot_given_data():
    q_oil, q_stm, q_water, tt0, X0 = interpolate_data() # get the interpolated data

    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()

    plt1 = ax1.plot(tt0, X0[0], "k.", label = "Pressure Data")
    plt2 = ax2.plot(tt0, X0[1], "r.", label = "Temperature Data")

    plts = plt1 + plt2
    labs = [l.get_label() for l in plts]
    ax1.legend(plts, labs)
    ax1.set_ylabel("Pressure (kPa)")
    ax2.set_ylabel("Temperature ($^\circ$C)")
    ax1.set_xlabel("time (days)")
    plt.show()


if __name__ == "__main__":
    plot_given_data()