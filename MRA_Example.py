#!/usr/bin/env python

"""
This module contains the discrete wavelet based multiresolution analysis (MRA),
used to determine circadian rhythmicities
"""

__author__ = "Christoph Schmal"
__license__ = "GPL"
__maintainer__ = "Christoph Schmal"
__email__ = "christoph.schmal@hu-berlin.de"
__repository__ = "https://github.com/cschmal/TNBC_ToD"

# plotting libraries
from pylab import*

# general libraries
import pandas as pd

#
import pyboat
from pyboat import sinc_smooth
# wavelet libraries
import pywt
# we make use of the modwt.py library from https://github.com/pistonly/modwtpy
from modwt import*

# Function that plots results of the MRA analysis
def PlotMRA(t, signal, MRA, dt=30./60, NumLevels=7, t_edge=12.):
    f, axarr = plt.subplots(NumLevels+2, sharex=True, figsize=(6/1.7, 12/1.7))
    c_max = ceil(max(signal.flatten()))

    i_edge = round(t_edge/dt)

    OverallEnergy = sum([sum(n[i_edge:-i_edge]**2) for n in MRA[:-1]])/len(t[i_edge:-i_edge])

    for i, j in enumerate(MRA):     # 7 elements 6
        axarr[i].plot(t, j, color="gray", linestyle="--")

        # calculate energy in detail j
        DetailEnergy = sum(j[i_edge:-i_edge]**2) / len(t[i_edge:-i_edge]) / OverallEnergy

        axarr[i].plot(t[i_edge:-i_edge], j[i_edge:-i_edge], label = str(round(DetailEnergy*100,2)) + r"% in " + r"["+str(round((2.**(i+1)*dt)))+"h, "+ str(round((2.**(i+2)*dt))) +"h[")
        axarr[i].set_xticks(arange(0, t[-1]+24., 24.))
        axarr[i].set_xlim(0, t[-1])
        axarr[i].set_ylabel("D$_{"+str(i+1) +"}$")
        leg2 = axarr[i].legend(loc="upper right", prop={'size':10}, fancybox=True, handlelength=0, handletextpad=0)
        leg2.get_frame().set_alpha(0.5)
    axarr[i+1].plot(t, signal, "k-", linewidth=2, label="Signal")
    axarr[i+1].plot(t, sum(MRA, axis=0), linestyle="--", label="$ D_{\sum} + S$")
    axarr[i+1].set_ylim(-c_max, c_max)
    axarr[i+1].set_yticks([-c_max, 0, c_max])
    leg1 = axarr[i+1].legend(loc="upper right", prop={'size':8}, fancybox=True, ncol=2)
    leg1.get_frame().set_alpha(0.5)

    axarr[i+1].set_xlabel("Time in h")
    f.tight_layout()

# Function that plots results of the MRA analysis with details merged into
# noise, infradian, circadian and ultradian regimes
def PlotMRA_CoarseGrained(t, signal, MRA, dt=30./60, NumLevels=7, t_edge=12.):
    f, axarr = plt.subplots(5, sharex=True, figsize=(6/1.7, 12/2.5))
    c_max = ceil(max(signal.flatten()))

    i_edge = round(t_edge/dt)

    OverallEnergy = sum([sum(n[i_edge:-i_edge]**2) for n in MRA[:-1]])/len(t[i_edge:-i_edge])
    D1Energy = sum(MRA[0][i_edge:-i_edge]**2) / len(t[i_edge:-i_edge]) / OverallEnergy
    D2Energy = sum(MRA[1][i_edge:-i_edge]**2) / len(t[i_edge:-i_edge]) / OverallEnergy
    D3Energy = sum(MRA[2][i_edge:-i_edge]**2) / len(t[i_edge:-i_edge]) / OverallEnergy
    D4Energy = sum(MRA[3][i_edge:-i_edge]**2) / len(t[i_edge:-i_edge]) / OverallEnergy
    D5Energy = sum(MRA[4][i_edge:-i_edge]**2) / len(t[i_edge:-i_edge]) / OverallEnergy
    D6Energy = sum(MRA[5][i_edge:-i_edge]**2) / len(t[i_edge:-i_edge]) / OverallEnergy
    D7Energy = sum(MRA[6][i_edge:-i_edge]**2) / len(t[i_edge:-i_edge]) / OverallEnergy
    D8Energy = sum(MRA[7][i_edge:-i_edge]**2) / len(t[i_edge:-i_edge]) / OverallEnergy


    D12  = MRA[0] + MRA[1]
    D34  = MRA[2] + MRA[3]
    D5   = MRA[4]
    D678 = MRA[5] + MRA[6] + MRA[7]

    c_label  = ["Noise", "Ultradian", "Circadian", "Infradian"]
    c_energy = [D1Energy+D2Energy, D3Energy+D4Energy, D5Energy, D6Energy+D7Energy+D8Energy]
    for i, j in enumerate([D12, D34, D5, D678]):
        axarr[i].plot(t, j, color="gray", linestyle="--")

        c_maxtmp = max(abs(j[i_edge:-i_edge]))
        axarr[i].plot(t[i_edge:-i_edge], j[i_edge:-i_edge], label = str(round(c_energy[i]*100,2)))
        axarr[i].set_xticks(arange(0, t[-1]+24., 24.))
        axarr[i].set_xlim(0, t[-1])
        axarr[i].set_ylim(-c_maxtmp, c_maxtmp)
        axarr[i].set_ylabel(c_label[i])
        leg2 = axarr[i].legend(loc="upper right", prop={'size':10}, fancybox=True, handlelength=0, handletextpad=0)
        leg2.get_frame().set_alpha(0.5)


    axarr[4].plot(t, signal, "k-", linewidth=2, label="Signal")
    axarr[4].plot(t, sum(MRA, axis=0), linestyle="--", label="$ D_{\sum} + S$")
    axarr[4].set_ylim(-c_max, c_max)
    axarr[4].set_yticks([-c_max, 0, c_max])
    leg1 = axarr[4].legend(loc="upper right", prop={'size':8}, fancybox=True, ncol=2)
    leg1.get_frame().set_alpha(0.5)

    axarr[4].set_xlabel("Time in h")
    f.tight_layout()


# Load Data
Data = pd.read_csv("./MRA_ExampleData.csv", sep=";", skipinitialspace=True)

# We analyze the first 6.5 days of the data recording
Data = Data.loc[Data["Time"] < 156]

# Load data into numpy arrays
t           = Data["Time"].values
signal      = Data.iloc[:,1].values

#
# Detrending of data
#

# sampling interval
dt = 10./60     # [h]
T_cut_off = 36  # cut-off period for sinc-filter


trend = sinc_smooth(signal, T_cut_off, dt)
detr_signal = signal - trend

#
# Start MRA
#

# Downsample data to obtain a proper circadian frequency band
t      = t[::3]
signal = detr_signal[::3]

# DWT parameters
c_wavelet = "db20"      # wavelet
NumLevels = 7           # number of details considered

# Apply MRA
LastOne = len(t)
t = array(t).copy()[:LastOne]
level = pywt.swt_max_level(LastOne)
c_modwt = modwt(signal, c_wavelet, NumLevels)
c_wtmra = modwtmra(c_modwt, c_wavelet)

# Plot MRA results
PlotMRA(t, signal, c_wtmra, dt=30./60, NumLevels=7, t_edge=12.)

# Plot MRA results with details merged into noise, infradian, circadian
# and ultradian regimes
PlotMRA_CoarseGrained(t, signal, c_wtmra, dt=30./60, NumLevels=7, t_edge=12.)

show()
