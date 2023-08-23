#!/usr/bin/env python

"""
Calculate correlations between drug sensitivity parameters as depicted in
Fig. 3o of the manuscript
"""

__author__ = "Christoph Schmal"
__license__ = "GPL"
__maintainer__ = "Christoph Schmal"
__email__ = "christoph.schmal@hu-berlin.de"
__repository__ = "https://github.com/cschmal/TNBC_ToD"


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

# Load data
file2load = 'ToD_Data.xlsx'

with open('PickledData'+file2load.split('.')[0]+'.pickle', 'rb') as f:
    data = pickle.load(f)
Circadian_Clock, Growth_Properties, Drug_Sensitivity, Time_of_Day_Sensitivity = data



# choose which parameters to load into data frame
df = pd.DataFrame()
df1 = pd.DataFrame()
df2 = pd.DataFrame()
df3 = pd.DataFrame()

# choose clock reporter data to be analyzed
clock_reporter = "Bmal1"
#clock_reporter = "Per2"
#clock_reporter = "Mixed"

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
Parameters = list(Drug_Sensitivity.keys())
DrugsSensitivity = list(Drug_Sensitivity[Parameters[0]][c].keys())[3:]
for j in Parameters:
    for i in DrugsSensitivity:
        df[i+"_"+j] = Drug_Sensitivity[j][c][i]
        df2[i+"_"+j] = Drug_Sensitivity[j][c][i]


#alternative sorting
Drugs = Time_of_Day_Sensitivity.keys()
#print(Drugs)
Parameters = Time_of_Day_Sensitivity["5FU"][c].keys()[3:]
for j in Parameters:
    for i in Drugs:
        if i != None:
            df[i+"_"+j] = Time_of_Day_Sensitivity[i][c][j]
            df3[i+"_"+j] = Time_of_Day_Sensitivity[i][c][j]
            count2 += 1


# select drug sensitivity parameters to be further analyzed
drug_keys = ["Paclitaxel_GRinf_ExpFit", "5FU_GRinf_ExpFit", "Alpelisib_GRinf_ExpFit", "Adavosertib_GRinf_ExpFit", "Doxorubicin_GRinf_ExpFit", "Alisertib_GRinf_ExpFit", "Torin2_GRinf_ExpFit", "Cisplatin_GRinf_ExpFit", "Olaparib_GRinf_ExpFit"]

# define list of drugs to be tested
Drug_List = ["Paclitaxel", "5FU", "Alpelisib", "Adavosertib", "Doxorubicin", "Alisertib", "Torin2", "Cisplatin"]

# select ToD observable to be considered
y_parameter, y_key_label = "cos_amplitude", "CosAmp"

y_keys = [i + "_" + y_parameter for i in Drug_List]

# choose drug parameters considered for analysis
DrugParameters = ["GRinf_ExpFit", "GR50_ExpFit", "GEC50_ExpFit", "Hill Coefficient_ExpFit", "EC50", "AOC_ExpFit"]

X_container = zeros(len(DrugParameters)*len(DrugParameters)).reshape(len(DrugParameters), len(DrugParameters))
for i, Drug in enumerate(Drug_List):
    drug_keys = [Drug + "_" + i for i in DrugParameters]

    X = df[drug_keys]

    X_alt, p_alt = r_pvalues(X)

    g = sns.clustermap(X.corr(), figsize=(6.4, 6.4), vmin=-1, vmax=1, cmap="vlag", standard_scale=None, annot=c_annot, annot_kws={"size": c_annot_size}, fmt=c_annot_fmt)

    X_container += X.corr().values
    x0, _y0, _w, _h = g.cbar_pos
    g.ax_cbar.set_position([0.75, 0.05, g.ax_row_dendrogram.get_position().width, 0.2])
    g.ax_cbar.set_title('correlation')
    g.ax_cbar.tick_params(axis='x', length=10)

    #savefig("./DrugCorrelations_Fig3o/"+Drug+".svg" )

    # re-order correlation and p-value matrix by dendongram ordering
    reordered_row_indices = g.dendrogram_row.reordered_ind
    reordered_col_indices = g.dendrogram_col.reordered_ind

    p_alt = p_alt.loc[array(p_alt.index.to_list())[reordered_row_indices]]
    p_alt = p_alt[array(p_alt.columns.to_list())[reordered_col_indices]]
    X_alt = X_alt.loc[array(X_alt.index.to_list())[reordered_row_indices]]
    X_alt = X_alt[array(X_alt.columns.to_list())[reordered_col_indices]]

# plot average correlation
df_X_container = pd.DataFrame(X_container/len(Drug_List), index=DrugParameters, columns=DrugParameters)
g = sns.clustermap(df_X_container, figsize=(6.4, 6.4), vmin=-1, vmax=1, cmap="vlag", standard_scale=None, annot=c_annot, annot_kws={"size": c_annot_size}, fmt=c_annot_fmt)
x0, _y0, _w, _h = g.cbar_pos
g.ax_cbar.set_position([0.75, 0.05, g.ax_row_dendrogram.get_position().width, 0.2])
g.ax_cbar.set_title('correlation')
g.ax_cbar.tick_params(axis='x', length=10)
#savefig("./DrugCorrelations_Fig3o/Dendogram_AverageDrugCorrelations.svg" )

show()
