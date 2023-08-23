#!/usr/bin/env python

"""
Calculate correlations between circadian parameters
"""

__author__ = "Christoph Schmal"
__license__ = "GPL"
__maintainer__ = "Christoph Schmal"
__email__ = "christoph.schmal@hu-berlin.de"


# plotting libraries
from pylab import*
import seaborn as sns

# general libraries
import pandas as pd
import pickle
from scipy.stats import pearsonr

# plotting parameters
c_annot = True
c_annot_fmt = ".2"
c_annot_size = 7

# defines a function to calculate correlations and corresponding p-values
# of a data frame
def r_pvalues(df):
    cols = pd.DataFrame(columns=df.columns)
    p = cols.transpose().join(cols, how='outer')
    correlation = p.copy()
    #print(p)
    for r in df.columns:
        for c in df.columns:
            tmp = df[df[r].notnull() & df[c].notnull()]
            p[r][c] = round(pearsonr(tmp[r], tmp[c])[1], 4)
            correlation[r][c] = round(pearsonr(tmp[r], tmp[c])[0], 4)
    return correlation, p

file2load = 'ToD_Data.xlsx'

with open('PickledData'+file2load.split('.')[0]+'.pickle', 'rb') as f:
    # The protocol version used is detected automatically, so we do not
    # have to specify it.
    data = pickle.load(f)

Circadian_Clock, Growth_Properties, Drug_Sensitivity, Time_of_Day_Sensitivity = data

# intialize data frames
df = pd.DataFrame()
df1 = pd.DataFrame()
df2 = pd.DataFrame()
df3 = pd.DataFrame()

clock_reporter = "Bmal1"
clock_reporter = "Per2"

# load clock parameters
count1 = 0
for i in Circadian_Clock.keys():
    if i == clock_reporter:
        for j in Circadian_Clock[i].keys()[3:]:
            df["_".join( (i,j.replace(" ", "_")))] = Circadian_Clock[i][j]
            df1["_".join( (i,j.replace(" ", "_")))] = Circadian_Clock[i][j]
            count1 += 1

c = "Cell Number"
count2 = 0

#alternative sorting
#print(Drug_Sensitivity.keys())
Parameters = list(Drug_Sensitivity.keys())
#print(list(Drug_Sensitivity[Parameters[0]][c].keys())[3:])
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


''' Drug_Sensitivity Regressions '''

#drug_keys = df2.keys()
drug_keys = ["Paclitaxel_GRinf_ExpFit", "5FU_GRinf_ExpFit", "Alpelisib_GRinf_ExpFit", "Adavosertib_GRinf_ExpFit", "Doxorubicin_GRinf_ExpFit", "Alisertib_GRinf_ExpFit", "Torin2_GRinf_ExpFit", "Cisplatin_GRinf_ExpFit", "Olaparib_GRinf_ExpFit"]

# define list of drugs to be tested
Drug_List = ["Paclitaxel", "5FU", "Alpelisib", "Adavosertib", "Doxorubicin", "Alisertib", "Torin2", "Cisplatin"]

# generate
#y_parameter, y_key_label = "maximum range", "MaxRange"
y_parameter, y_key_label = "cos_amplitude", "CosAmp"
#y_parameter, y_key_label = "cos_Rsq", "CosRsq"

y_keys = [i + "_" + y_parameter for i in Drug_List]


DrugParameters = ["GRinf_ExpFit", "GR50_ExpFit", "GEC50_ExpFit", "Hill Coefficient_ExpFit", "EC50", "AOC_ExpFit"]

circadian_keys = ["circadian_mean", "noise_mean", "period_median_2d", "period_coeffvar_2d", "period_median_5d", "period_stdev_5d", "amplitude_median_5d", "amplitude_coeffvar_5d", "autocorr_period_median", "autocorr_peak_median", "ridgelength_median"]

if clock_reporter == "Per2":
    circadian_keys = ["circadian_mean", "noise_mean", "period_median_2d", "period_coeffvar_2d", "period_median_5d", "period_stdev_5d", "amplitude_median_5d", "amplitude_coeffvar_5d", "autocorr_period", "autocorr_peak", "ridgelength_median"]
elif clock_reporter == "Mixed":
    circadian_keys = ["circadian_mean", "noise_mean", "period_median_2d", "period_coeffvar_2d", "amplitude_median_5d", "amplitude_coeffvar_5d", "autocorr_period", "autocorr_peak", "ridgelength_median"]

circadian_keys = [clock_reporter+"_"+i for i in circadian_keys ]

X = df[circadian_keys]

print(X.corr())

g = sns.clustermap(X.corr(), figsize=(6.4, 6.4), vmin=-1, vmax=1, cmap="vlag", standard_scale=None, annot=c_annot, annot_kws={"size": c_annot_size}, fmt=c_annot_fmt)
x0, _y0, _w, _h = g.cbar_pos
g.ax_cbar.set_position([0.75, 0.05, g.ax_row_dendrogram.get_position().width, 0.2])
g.ax_cbar.set_title('correlation')
g.ax_cbar.tick_params(axis='x', length=10)

#savefig("./ClockParameterCorrelations/"+clock_reporter+".png" )
#savefig("./ClockParameterCorrelations/"+clock_reporter+".svg" )


# Write correlation values and p-values to a csv
X_alt, p_alt = r_pvalues(X)
reordered_row_indices = g.dendrogram_row.reordered_ind
reordered_col_indices = g.dendrogram_col.reordered_ind
p_alt = p_alt.loc[array(p_alt.index.to_list())[reordered_row_indices]]
p_alt = p_alt[array(p_alt.columns.to_list())[reordered_col_indices]]
X_alt = X_alt.loc[array(X_alt.index.to_list())[reordered_row_indices]]
X_alt = X_alt[array(X_alt.columns.to_list())[reordered_col_indices]]
#X_alt.to_csv("./ClockParameterCorrelations/Data_Correlation_"+clock_reporter+".csv", header=True, index=True)
#p_alt.to_csv("."./ClockParameterCorrelations/Data_PValue_"+clock_reporter+".csv", header=True, index=True)



show()
