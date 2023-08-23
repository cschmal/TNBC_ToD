#!/usr/bin/env python

"""
Calculate correlations between circadian or growth or drug sensitivity parameters and
drug dependent time-of-day responses together with p-values
as depicted in Fig. 5d-f of the manuscript
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
df4 = pd.DataFrame()

# choose which parameters to load into data frame
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
Parameters = list(Drug_Sensitivity.keys())
DrugsSensitivity = list(Drug_Sensitivity[Parameters[0]][c].keys())[3:]
for j in Parameters:
    for i in DrugsSensitivity:
        df[i+"_"+j] = Drug_Sensitivity[j][c][i]
        df2[i+"_"+j] = Drug_Sensitivity[j][c][i]
        count2 += 1


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

color_index = {k:i for i, k in enumerate(df.index.values.tolist())}

circadian_keys = df1.keys()

df = df.drop(index="HCC1937_2")
df1 = df1.drop(index="HCC1937_2")
df2 = df2.drop(index="HCC1937_2")
df3 = df3.drop(index="HCC1937_2")
df4 = df4.drop(index="HCC1937_2")

Drug_List = ["Paclitaxel", "5FU", "Alpelisib", "Adavosertib", "Doxorubicin", "Alisertib", "Torin2", "Cisplatin"]

y_dependentvariable = "maxrange_spline"

y_list = [ y_drug + "_" + y_dependentvariable for y_drug in Drug_List]

# select reduced drug list
df3_reduced = df3[df3.columns.intersection(y_list)]

if clock_reporter != "Mixed":
    df1 = df1.drop(columns= clock_reporter + "_" + "period_median_5d" )
    df1 = df1.drop(columns= clock_reporter + "_" + "period_stdev_5d" )
df1 = df1.drop(columns= clock_reporter + "_" + "period_coeffvar_2d" )
df1 = df1.drop(columns= clock_reporter + "_" + "amplitude_coeffvar_5d" )

count1 = len(df1.keys())

djoin = pd.concat([df1, df3_reduced], axis=1, join='inner')

corr = djoin.corr()
corr = corr.iloc[0:count1,count1:]

#corr.to_csv("./Correlation_Circadian_"+clock_reporter+"_"+y_dependentvariable+".csv", header=True, index=True)

g = sns.clustermap(corr.T, figsize=(6.4, 6.4), vmin=-1, vmax=1, cmap="vlag", standard_scale=None, annot=c_annot, annot_kws={"size": c_annot_size}, fmt=c_annot_fmt)
x0, _y0, _w, _h = g.cbar_pos
g.ax_cbar.set_position([0.75, 0.05, g.ax_row_dendrogram.get_position().width, 0.2])
g.ax_cbar.set_title('correlation')
g.ax_cbar.tick_params(axis='x', length=10)

tight_layout()

#savefig("./Circadian_correlation_"+clock_reporter+"_"+y_dependentvariable+".png" )
#savefig("./Circadian_correlation_"+clock_reporter+"_"+y_dependentvariable+".svg" )

g = sns.clustermap(corr, figsize=(6.4, 6.4), vmin=-1, vmax=1, cmap="vlag", standard_scale=None, annot=c_annot, annot_kws={"size": c_annot_size}, fmt=c_annot_fmt)


def r_pvalues(df):
    cols = pd.DataFrame(columns=df.columns)
    p = cols.transpose().join(cols, how='outer')
    correlation = p.copy()
    for r in df.columns:
        for c in df.columns:
            tmp = df[df[r].notnull() & df[c].notnull()]
            p[r][c] = round(pearsonr(tmp[r], tmp[c])[1], 4)
            correlation[r][c] = round(pearsonr(tmp[r], tmp[c])[0], 4)
    return correlation, p

c_matrix, p_matrix = r_pvalues(djoin)
p_matrix = p_matrix.iloc[0:count1,count1:]
p_matrix = p_matrix.apply(pd.to_numeric)

c_matrix = c_matrix.iloc[0:count1,count1:]
c_matrix = c_matrix.apply(pd.to_numeric)


reordered_row_indices = g.dendrogram_row.reordered_ind
reordered_col_indices = g.dendrogram_col.reordered_ind

p_matrix = p_matrix.loc[array(corr.index.to_list())[reordered_row_indices]]
p_matrix = p_matrix[array(corr.columns.to_list())[reordered_col_indices]]

#p_matrix.to_csv("./Correlation_Circadian_"+clock_reporter+"_"+y_dependentvariable+"_pvalues.csv", header=True, index=True)

f, ax = plt.subplots(figsize=(6.4, 6.4))
g2 = sns.heatmap(p_matrix.T, cmap="spring", center=0,
            square=True, linewidths=.5, cbar_kws={"shrink": .5}, annot=c_annot, annot_kws={"size": c_annot_size}, fmt=c_annot_fmt, vmin=0, vmax=1)

tight_layout()

#savefig("./Circadian_pvalue_"+clock_reporter+"_"+y_dependentvariable+".png" )
#savefig("./Circadian_pvalue_"+clock_reporter+"_"+y_dependentvariable+".svg" )

DrugParameters = ["GRinf_ExpFit", "GR50_ExpFit", "GEC50_ExpFit", "Hill Coefficient_ExpFit", "AOC_ExpFit"]

corr_matrix = zeros(len(Drug_List)*len(DrugParameters)).reshape(len(Drug_List), len(DrugParameters))
p_matrix = zeros(len(Drug_List)*len(DrugParameters)).reshape(len(Drug_List), len(DrugParameters))

for n, c_drug in enumerate(Drug_List):
    c_index = [c_drug + "_" + i for i in DrugParameters]
    df2_tmp = df2[c_index]
    y = df3_reduced[c_drug + "_" + y_dependentvariable]
    c_corr_tmp = zeros( len(df2_tmp.keys()) )
    c_p_tmp = zeros( len(df2_tmp.keys()) )
    for m, c_param in enumerate(df2_tmp.keys()):
        x = df2_tmp[c_param]

        x_no_nan = x[~x.isnull()]
        y_no_nan = y[~x.isnull()]

        # filter data based on nans in y
        x_no_nan = x_no_nan[invert(isnan(y_no_nan))]
        y_no_nan = y_no_nan[invert(isnan(y_no_nan))]

        c_corr = round(pearsonr(x_no_nan, y_no_nan)[0], 6)
        c_p = round(pearsonr(x_no_nan, y_no_nan)[1], 6)

        c_corr_tmp[m] = copy(c_corr)
        c_p_tmp[m]    = copy(c_p)
    corr_matrix[n] = c_corr_tmp
    p_matrix[n] = c_p_tmp

df_corr = pd.DataFrame(corr_matrix, index=Drug_List, columns=DrugParameters)
df_p = pd.DataFrame(p_matrix, index=Drug_List, columns=DrugParameters)

#df_corr.to_csv("./Correlation_DrugResponse_"+y_dependentvariable+".csv", header=True, index=True)
#df_p.to_csv("./Correlation_DrugResponse_"+y_dependentvariable+"_pvalue.csv", header=True, index=True)



g = sns.clustermap(df_corr.T, figsize=(6.4, 4.4), vmin=-1, vmax=1, cmap="vlag", standard_scale=None, annot=c_annot, annot_kws={"size": c_annot_size}, fmt=c_annot_fmt)
x0, _y0, _w, _h = g.cbar_pos
g.ax_cbar.set_position([0.75, 0.05, g.ax_row_dendrogram.get_position().width, 0.1])
g.ax_cbar.set_title('correlation')
g.ax_cbar.tick_params(axis='x', length=10)

tight_layout()

#savefig("./DrugResponse_correlation"+"_"+y_dependentvariable+".png" )
#savefig("./DrugResponse_correlation"+"_"+y_dependentvariable+".svg" )

g = sns.clustermap(df_corr, figsize=(6.4, 4.4), vmin=-1, vmax=1, cmap="vlag", standard_scale=None, annot=c_annot, annot_kws={"size": c_annot_size}, fmt=c_annot_fmt)


reordered_row_indices = g.dendrogram_row.reordered_ind
reordered_col_indices = g.dendrogram_col.reordered_ind

df_p = df_p.loc[array(df_corr.index.to_list())[reordered_row_indices]]
df_p = df_p[array(df_corr.columns.to_list())[reordered_col_indices]]


f, ax = plt.subplots(figsize=(6.4, 4.4))
g2 = sns.heatmap(df_p.T, cmap="spring", center=0,
            square=True, linewidths=.5, cbar_kws={"shrink": .5}, annot=c_annot, annot_kws={"size": c_annot_size}, fmt=c_annot_fmt, vmin=0, vmax=1)

tight_layout()

#savefig("./DrugResponse_pvalue"+"_"+y_dependentvariable+".png" )
#savefig("./DrugResponse_pvalue"+"_"+y_dependentvariable+".svg" )


djoin = pd.concat([df4, df3_reduced], axis=1, join='inner')

count4 = len(df4.keys())

corr = djoin.corr()
corr = corr.iloc[0:count4,count4:]

#corr.to_csv("./Correlation_GrowthProperties_"+y_dependentvariable+".csv", header=True, index=True)

g = sns.clustermap(corr.T, figsize=(6.4, 6.4), vmin=-1, vmax=1, cmap="vlag", standard_scale=None, annot=c_annot, annot_kws={"size": c_annot_size}, fmt=c_annot_fmt)
x0, _y0, _w, _h = g.cbar_pos
g.ax_cbar.set_position([0.75, 0.05, g.ax_row_dendrogram.get_position().width, 0.2])
g.ax_cbar.set_title('correlation')
g.ax_cbar.tick_params(axis='x', length=10)

tight_layout()

#savefig("./GrowthProperties_correlation"+"_"+y_dependentvariable+".png" )
#savefig("./GrowthProperties_correlation"+"_"+y_dependentvariable+".svg" )

c_matrix, p_matrix = r_pvalues(djoin)
p_matrix = p_matrix.iloc[0:count4,count4:]
p_matrix = p_matrix.apply(pd.to_numeric)
c_matrix = c_matrix.iloc[0:count4,count4:]
c_matrix = c_matrix.apply(pd.to_numeric)

reordered_row_indices = g.dendrogram_row.reordered_ind
reordered_col_indices = g.dendrogram_col.reordered_ind

p_matrix = p_matrix.loc[array(corr.index.to_list())[reordered_col_indices]]
p_matrix = p_matrix[array(corr.columns.to_list())[reordered_row_indices]]

#corr.to_csv("./Correlation_GrowthProperties_"+y_dependentvariable+"_pvalues.csv", header=True, index=True)

f, ax = plt.subplots(figsize=(6.4, 6.4))
g2 = sns.heatmap(p_matrix.T, cmap="spring", center=0,
            square=True, linewidths=.5, cbar_kws={"shrink": .5}, annot=c_annot, annot_kws={"size": c_annot_size}, fmt=c_annot_fmt, vmin=0, vmax=1)


tight_layout()

#savefig("./GrowthProperties_pvalue"+"_"+y_dependentvariable+".png" )
#savefig("./GrowthProperties_pvalue"+"_"+y_dependentvariable+".svg" )



show()
