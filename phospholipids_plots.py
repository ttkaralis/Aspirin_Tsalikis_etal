import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

data = pd.read_csv('Negative_Mode_TICnormnalised_Lipidomics_THanasis_300625.txt', sep = '\t', index_col = 0)
data

data.columns

data = data[['QC negative 1: Log2', 'QC negative 2: Log2', 'QC negative 3: Log2', 'QC negative 4: Log2',

             'WT 1to1 negative 1: Log2','WT 1to1 negative 2: Log2', 'WT 1to1 negative 3: Log2', 'WT 1to1 negative 4: Log2', 'WT 1to1 negative 5: Log2',
             'WT 1to2 negative 1: Log2', 'WT 1to2 negative 2: Log2', 'WT 1to2 negative 3: Log2', 'WT 1to2 negative 4: Log2', 'WT 1to2 negative 5: Log2',
             'WT 10to1 negative 1: Log2', 'WT 10to1 negative 2: Log2', 'WT 10to1 negative 3: Log2', 'WT 10to1 negative 4: Log2', 'WT 10to1 negative 5: Log2',

             'H10 1to1 negative 1: Log2', 'H10 1to1 negative 2: Log2', 'H10 1to1 negative 3: Log2', 'H10 1to1 negative 4: Log2', 'H10 1to1 negative 5: Log2',
             'H10 1to2 negative 1: Log2', 'H10 1to2 negative 2: Log2', 'H10 1to2 negative 3: Log2', 'H10 1to2 negative 4: Log2', 'H10 1to2 negative 5: Log2',
             'H10 10to1 negative 1: Log2', 'H10 10to1 negative 2: Log2', 'H10 10to1 negative 3: Log2', 'H10 10to1 negative 4: Log2', 'H10 10to1 negative 5: Log2']]

data = data.drop(columns = ['QC negative 1: Log2',
       'QC negative 2: Log2', 'QC negative 3: Log2', 'QC negative 4: Log2',])
data

data.columns = data.columns.str.replace(' negative', '').str.replace(' ', ('_')).str.replace(':_Log2', '')
data

mutation = [
            'WT', 'WT', 'WT', 'WT', 'WT', 'WT', 'WT', 'WT', 'WT', 'WT', 'WT', 'WT', 'WT', 'WT', 'WT',
            'H1047R', 'H1047R', 'H1047R', 'H1047R', 'H1047R', 'H1047R', 'H1047R', 'H1047R', 'H1047R', 'H1047R', 'H1047R', 'H1047R', 'H1047R', 'H1047R', 'H1047R',]
diet = [
        '1_1', '1_1', '1_1', '1_1', '1_1',
        '1_2', '1_2', '1_2', '1_2', '1_2',
        '10_1', '10_1', '10_1', '10_1', '10_1',

        '1_1', '1_1', '1_1', '1_1', '1_1',
        '1_2', '1_2', '1_2', '1_2', '1_2',
        '10_1', '10_1', '10_1', '10_1', '10_1']

data = data.T
data

data['PIK3CA_mutation'] = mutation
data['Diet'] = diet
data

import statsmodels.api as sm
from statsmodels.formula.api import ols

anova_results = []

for compound in data.columns[:-2]:
    compound_df = data[['PIK3CA_mutation', 'Diet', compound]].copy()

    cleaned_compound = compound.replace(' ', '_').replace(':', '_').replace('/', '_').replace('(', '_').replace(')', '_')
    formula = f'{cleaned_compound} ~ C(PIK3CA_mutation) + C(Diet) + C(PIK3CA_mutation):C(Diet)'

    compound_df = compound_df.rename(columns={compound: cleaned_compound})

    model = ols(formula, data=compound_df).fit()

    anova_table = sm.stats.anova_lm(model, typ=2)

    p_mutation = anova_table['PR(>F)']['C(PIK3CA_mutation)']
    p_diet = anova_table['PR(>F)']['C(Diet)']
    p_interaction = anova_table['PR(>F)']['C(PIK3CA_mutation):C(Diet)']

    anova_results.append({'Compound': compound,
                          'p_value_mutation': p_mutation,
                          'p_value_diet': p_diet,
                          'p_value_interaction': p_interaction})

anova_results_df = pd.DataFrame(anova_results)
display(anova_results_df)

from statsmodels.sandbox.stats.multicomp import multipletests

reject_mutation, pvals_corrected_mutation, _, _ = multipletests(anova_results_df['p_value_mutation'], method='fdr_i')
anova_results_df['p_value_mutation_fdr'] = pvals_corrected_mutation

reject_diet, pvals_corrected_diet, _, _ = multipletests(anova_results_df['p_value_diet'], method='fdr_i')
anova_results_df['p_value_diet_fdr'] = pvals_corrected_diet


reject_interaction, pvals_corrected_interaction, _, _ = multipletests(anova_results_df['p_value_interaction'], method='fdr_i')
anova_results_df['p_value_interaction_fdr'] = pvals_corrected_interaction


display(anova_results_df.head())

anova_results_df_filtered = anova_results_df[anova_results_df['p_value_diet_fdr'] < 0.05]
anova_results_df_filtered

data_sig = data[anova_results_df_filtered.Compound]
data_sig = pd.concat([data_sig, data[['PIK3CA_mutation', 'Diet']]], axis = 1)
data_sig

lipids = ['20:4', '20:5', '22:6']
select_cols = []
for i in data_sig.columns[:-2]:
  for j in lipids:
    if j in i:
      if i not in select_cols:
        select_cols.append(i)
len(select_cols)

select_cols

data_sig = data_sig[select_cols]
data_sig = pd.concat([data_sig, data[['PIK3CA_mutation', 'Diet']]], axis = 1)
data_sig

data_sig.columns

from sklearn.preprocessing import StandardScaler

data_scaled = StandardScaler().fit_transform(data_sig.drop(columns = ['PIK3CA_mutation', 'Diet']))
data_scaled = pd.DataFrame(data_scaled, columns = data_sig.drop(columns = ['PIK3CA_mutation', 'Diet']).columns, index=data_sig.index)
data_scaled

from matplotlib.colors import LinearSegmentedColormap

df_plot = data_scaled.transpose()
colors = ["#00008B", "white", "#8B0000"]  # blue → white → red
cmap = LinearSegmentedColormap.from_list("my_cmap", colors)

plt.figure(figsize = (10, 10))
sns.clustermap(df_plot, cmap = cmap,col_cluster=False)

metadata = data_sig[['Diet', 'PIK3CA_mutation']]
metadata = metadata.replace({"1_1": "1:1",
    "1_2": "1:2",
    "10_1": "10:1"})
metadata

metadata = data_sig[['Diet', 'PIK3CA_mutation']]
metadata = metadata.replace({"1_1": "1:1",
    "1_2": "1:2",
    "10_1": "10:1"})
metadata

metabolite_names = data_scaled.columns.to_list()
metabolite_metadata = {}
for i in metabolite_names:
  if '20:4' in i:
    metabolite_metadata[i] = 'AA'
  else:
    metabolite_metadata[i] = 'EPA/DHA'
metabolite_metadata = pd.DataFrame.from_dict(metabolite_metadata, orient = 'index', columns = ['Lipid sidechains'])
metabolite_metadata

lut = {
    "1:1": "#4B6CB7",    # medium blue
    "1:2": "#8B0000",    # dark red
    "10:1": "#B23A48",   # warm muted red
    "WT": "#B0B0B0",     # neutral grey
    "H1047R": "#1E3A8A"  # deep navy blue
}

lut_1 = {
    "AA": "#2F4F4F",     # dark slate gray (neutral, cool)
    "EPA/DHA": "#C0A080" # warm beige-tan (slightly warm to echo the red side)
}

# Create a DataFrame for col_colors with both 'Diet' and 'PIK3CA_mutation'
col_colors = metadata.apply(lambda x: x.map(lut))
row_colors = metabolite_metadata.apply(lambda x: x.map(lut_1))

g = sns.clustermap(df_plot, cmap = cmap,col_cluster=False, col_colors = col_colors, row_colors = row_colors, linewidths=0)
plt.ylabel('Metabolite levels')
g.ax_heatmap.set_xlabel('')
g.ax_heatmap.set_ylabel('')
g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)

# Create separate legends for Diet and PIK3CA_mutation
import matplotlib.patches as mpatches

diet_patches = [mpatches.Patch(color=lut[label], label=label) for label in metadata['Diet'].unique()]
mutation_patches = [mpatches.Patch(color=lut[label], label=label) for label in metadata['PIK3CA_mutation'].unique()]
lipid_patches = [mpatches.Patch(color=lut_1[label], label=label) for label in metabolite_metadata['Lipid sidechains'].unique()]


# Add the legends to the figure
g.fig.legend(handles=diet_patches, title="Diet Ratios", loc="upper left", bbox_to_anchor=(1, 1), ncol=1)
g.fig.legend(handles=mutation_patches, title="PIK3CA Mutation", loc="upper right", bbox_to_anchor=(1, 1), ncol=1)
g.fig.legend(handles=lipid_patches, title="Lipid Sidechains", loc="upper left", bbox_to_anchor=(1, 0.8), ncol=1)

plt.savefig('phospholipids_heatmap.svg', bbox_inches = 'tight')

plt.show()

lut = {
    "1:1": "#457B9D",    # muted steel blue (cool)
    "1:2": "#E63946",    # rich coral red (warm)
    "10:1": "#B5838D",   # dusty rose (warm but soft)
    "WT": "#A8A8A8",     # neutral mid-grey
    "H1047R": "#264653"  # deep teal-blue (cool but distinct)
}


lut_1 = {
    "AA": "#4A4E69",     # muted indigo-grey (neutral cool)
    "EPA/DHA": "#F4A261" # warm desaturated orange-beige
}

# Create a DataFrame for col_colors with both 'Diet' and 'PIK3CA_mutation'
col_colors = metadata.apply(lambda x: x.map(lut))
row_colors = metabolite_metadata.apply(lambda x: x.map(lut_1))

g = sns.clustermap(df_plot, cmap = cmap,col_cluster=False, col_colors = col_colors, row_colors = row_colors, linewidths=0)
plt.ylabel('Metabolite levels')
g.ax_heatmap.set_xlabel('')
g.ax_heatmap.set_ylabel('')
g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)

# Create separate legends for Diet and PIK3CA_mutation
import matplotlib.patches as mpatches

diet_patches = [mpatches.Patch(color=lut[label], label=label) for label in metadata['Diet'].unique()]
mutation_patches = [mpatches.Patch(color=lut[label], label=label) for label in metadata['PIK3CA_mutation'].unique()]
lipid_patches = [mpatches.Patch(color=lut_1[label], label=label) for label in metabolite_metadata['Lipid sidechains'].unique()]


# Add the legends to the figure
g.fig.legend(handles=diet_patches, title="Diet Ratios", loc="upper left", bbox_to_anchor=(1, 1), ncol=1)
g.fig.legend(handles=mutation_patches, title="PIK3CA Mutation", loc="upper right", bbox_to_anchor=(1, 1), ncol=1)
g.fig.legend(handles=lipid_patches, title="Lipid Sidechains", loc="upper left", bbox_to_anchor=(1, 0.8), ncol=1)

plt.savefig('phospholipids_heatmap_1.svg', bbox_inches = 'tight')

plt.show()

min_val = np.min(df_plot)
max_val = np.max(df_plot)
data_scaled = 2 * ((df_plot - min_val) / (max_val - min_val)) - 1
print(np.max(data_scaled))
print(np.min(data_scaled))
data_scaled.head(2)

colors = ["#00008B", "white", "#8B0000"]  # blue → white → red
cmap = LinearSegmentedColormap.from_list("my_cmap", colors)

lut = {
    "1:1": "#323232",    # muted steel blue (cool)
    "1:2": "#0072A1",    # rich coral red (warm)
    "10:1": "#D3A111",   # dusty rose (warm but soft)
    "WT": "#921A1D",     # neutral mid-grey
    "H1047R": "#2B3283"  # deep teal-blue (cool but distinct)
}


lut_1 = {
    "AA": "#7B6423",     # muted indigo-grey (neutral cool)
    "EPA/DHA": "#17556F" # warm desaturated orange-beige
}

# Create a DataFrame for col_colors with both 'Diet' and 'PIK3CA_mutation'
col_colors = metadata.apply(lambda x: x.map(lut))
row_colors = metabolite_metadata.apply(lambda x: x.map(lut_1))

g = sns.clustermap(data_scaled, cmap = cmap,col_cluster=False, col_colors = col_colors, row_colors = row_colors, linewidths=0)
plt.ylabel('Metabolite levels')
g.ax_heatmap.set_xlabel('')
g.ax_heatmap.set_ylabel('')
g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)

# Create separate legends for Diet and PIK3CA_mutation
import matplotlib.patches as mpatches

diet_patches = [mpatches.Patch(color=lut[label], label=label) for label in metadata['Diet'].unique()]
mutation_patches = [mpatches.Patch(color=lut[label], label=label) for label in metadata['PIK3CA_mutation'].unique()]
lipid_patches = [mpatches.Patch(color=lut_1[label], label=label) for label in metabolite_metadata['Lipid sidechains'].unique()]


# Add the legends to the figure
g.fig.legend(handles=diet_patches, title="Diet Ratios", loc="upper left", bbox_to_anchor=(1, 1), ncol=1)
g.fig.legend(handles=mutation_patches, title="PIK3CA Mutation", loc="upper right", bbox_to_anchor=(1, 1), ncol=1)
g.fig.legend(handles=lipid_patches, title="Lipid Sidechains", loc="upper left", bbox_to_anchor=(1, 0.8), ncol=1)

plt.savefig('phospholipids_heatmap_2.svg', bbox_inches = 'tight')

plt.show()

