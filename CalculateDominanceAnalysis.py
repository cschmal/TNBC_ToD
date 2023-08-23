#!/usr/bin/env python

"""
Apply dominance analysis to obtain data
underlying Fig. 5i-j of the manuscript
"""

__author__ = "Christoph Schmal"
__license__ = "GPL"
__maintainer__ = "Christoph Schmal"
__email__ = "christoph.schmal@hu-berlin.de"
__repository__ = "https://github.com/cschmal/TNBC_ToD"


# import plotting libraries
from pylab import*
import seaborn as sns

# general libraries
import pickle

# we make use of the dominance analysis package found at
# https://dominance-analysis.github.io/dominance-analysis/
from dominance import*

# define plotting paramerers
c_annot = True
c_annot_fmt = ".2"
c_annot_size = 7

# Load pickled data
file2load = 'ToD_Data.xlsx'

with open('PickledData'+file2load.split('.')[0]+'.pickle', 'rb') as f:
    data = pickle.load(f)

Circadian_Clock, Growth_Properties, Drug_Sensitivity, Time_of_Day_Sensitivity = data

# choose which parameters to load into data frame
df = pd.DataFrame()
df1 = pd.DataFrame()
df2 = pd.DataFrame()
df3 = pd.DataFrame()
df4 = pd.DataFrame()

clock_reporter = "Bmal1"
#clock_reporter = "Per2"
#clock_reporter = "Mixed"

# load clock parameters
count1 = 0
for i in Circadian_Clock.keys():
    if i == clock_reporter:
        print(Circadian_Clock[i].keys()[3:])
        for j in Circadian_Clock[i].keys()[3:]:
            print("_".join( (i,j.replace(" ", "_")) ))
            df["_".join( (i,j.replace(" ", "_")))] = Circadian_Clock[i][j]
            df1["_".join( (i,j.replace(" ", "_")))] = Circadian_Clock[i][j]
            count1 += 1

c = "Cell Number"
count2 = 0


Parameters = list(Drug_Sensitivity.keys())
DrugsSensitivity = list(Drug_Sensitivity[Parameters[0]][c].keys())[3:]
for j in Parameters:
    for i in DrugsSensitivity:
        df[i+"_"+j] = Drug_Sensitivity[j][c][i]
        df2[i+"_"+j] = Drug_Sensitivity[j][c][i]


#alternative sorting
Drugs = Time_of_Day_Sensitivity.keys()
Parameters = Time_of_Day_Sensitivity["5FU"][c].keys()[3:]
for j in Parameters:
    for i in Drugs:
        if i != None:
            df[i+"_"+j] = Time_of_Day_Sensitivity[i][c][j]
            df3[i+"_"+j] = Time_of_Day_Sensitivity[i][c][j]
            count2 += 1

# load-in growth properties
df["DoublingTime"] = Growth_Properties["Doubling Time [h]"][c]
df4["DoublingTime"] = Growth_Properties["Doubling Time [h]"][c]

df["GrowthRateExpFit"] = Growth_Properties["Growth Rate ExpFit"][c]
df4["GrowthRateExpFit"] = Growth_Properties["Growth Rate ExpFit"][c]

df.rename(index=Circadian_Clock["Mixed"]["Cellline"], inplace=True)
df1.rename(index=Circadian_Clock["Mixed"]["Cellline"], inplace=True)
df2.rename(index=Circadian_Clock["Mixed"]["Cellline"], inplace=True)
df3.rename(index=Circadian_Clock["Mixed"]["Cellline"], inplace=True)
df4.rename(index=Circadian_Clock["Mixed"]["Cellline"], inplace=True)

# drop second (unreliable) HCC1937_2 cell line data
df = df.drop(index="HCC1937_2")
df1 = df1.drop(index="HCC1937_2")
df2 = df2.drop(index="HCC1937_2")
df3 = df3.drop(index="HCC1937_2")
df4 = df4.drop(index="HCC1937_2")

circadian_keys = df1.keys()


''' Drug_Sensitivity Regressions '''

drug_keys = ["Paclitaxel_GRinf_ExpFit", "5FU_GRinf_ExpFit", "Alpelisib_GRinf_ExpFit", "Adavosertib_GRinf_ExpFit", "Doxorubicin_GRinf_ExpFit", "Alisertib_GRinf_ExpFit", "Torin2_GRinf_ExpFit", "Cisplatin_GRinf_ExpFit", "Olaparib_GRinf_ExpFit"]

# define list of drugs to be tested
Drug_List = ["Paclitaxel", "5FU", "Alpelisib", "Adavosertib", "Doxorubicin", "Alisertib", "Torin2", "Cisplatin"]

# generate
#y_parameter, y_key_label, c_title = "maximum range", "MaxRange", "maximum range"
#y_parameter, y_key_label, c_title = "cos_amplitude", "CosAmp", "cosinus amplitude"
y_parameter, y_key_label, c_title = "maxrange_spline", "MRSpline", "maxrange spline"

y_keys = [i + "_" + y_parameter for i in Drug_List]

DrugParameters = ["GEC50_ExpFit", "GRinf_ExpFit", "Hill Coefficient_ExpFit"]



circadian_keys = ["circadian_mean", "period_median_2d", "amplitude_median_5d"]
circadian_keys = [clock_reporter+"_"+i for i in circadian_keys ]

# Calculate importance
y_container = zeros(len(y_keys)*(len(DrugParameters) + 1 + len(circadian_keys))).reshape(len(y_keys),len(DrugParameters) + 1 + len(circadian_keys))
for n, y_key in enumerate(y_keys):
    # read in drug to analyze
    c_drug = y_key.split("_")[0]
    drug_keys = [c_drug + "_" + k for k in DrugParameters]
    combined_keys = list(circadian_keys) + drug_keys + ["GrowthRateExpFit"]

    # select subset based on keys
    X = df[combined_keys]
    y = df[y_key].values

    # filter data based on nans in X
    X_no_nan = X[~X.isnull().any(axis=1)]
    y_no_nan = y[~X.isnull().any(axis=1)]

    # filter data based on nans in y
    X_no_nan = X_no_nan[invert(isnan(y_no_nan))]
    y_no_nan = y_no_nan[invert(isnan(y_no_nan))]

    X_final = X_no_nan[:]
    X_final[y_key] = y_no_nan

    dominance_regression=Dominance(data=X_final,target=y_key,objective=1)
    incr_variable_rsquare=dominance_regression.incremental_rsquare()
    y_container[n] = list(incr_variable_rsquare.values())


x_keys = list(circadian_keys) + DrugParameters +  ["GrowthRateExpFit"]
data = x_keys, y_keys, y_container

x_keys, y_keys, y_container = data
y_keys = [i.split("_")[0] + " ($R^2=$" + str(round(j, 2)) + ")" for i, j in zip(y_keys, sum(y_container, axis=1))]

# plot total dominance
df = pd.DataFrame(y_container, columns=x_keys, index=y_keys)
g = sns.clustermap(df, figsize=(6.4, 4.8), cmap="viridis", standard_scale=None, annot=c_annot, annot_kws={"size": c_annot_size}, fmt=c_annot_fmt)
x0, _y0, _w, _h = g.cbar_pos
g.ax_cbar.set_position([0.75, 0.05, g.ax_row_dendrogram.get_position().width, 0.2])
g.ax_cbar.set_title('total dominance')
g.ax_cbar.tick_params(axis='x', length=10)
g.fig.suptitle(c_title, x=0.98, horizontalalignment="right")


#normalize by diving max value per row to get relative dominance
y_container = (y_container.T / sum(y_container, axis=1)).T

df = pd.DataFrame(y_container, columns=x_keys, index=y_keys)
g = sns.clustermap(df, figsize=(6.4, 4.8), cmap="viridis", standard_scale=None, annot=c_annot, annot_kws={"size": c_annot_size}, fmt=c_annot_fmt)
x0, _y0, _w, _h = g.cbar_pos
g.ax_cbar.set_position([0.75, 0.05, g.ax_row_dendrogram.get_position().width, 0.2])
g.ax_cbar.set_title('relative dominance')
g.ax_cbar.tick_params(axis='x', length=10)
g.fig.suptitle(c_title, x=0.98, horizontalalignment="right")


show()
