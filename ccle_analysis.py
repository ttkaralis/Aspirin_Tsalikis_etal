#Import modules

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

#Data

rna = pd.read_csv("CCLE_24Q4/Batch_corrected_Expression_Public_24Q4_subsetted.csv", index_col = 0)
rna

crispr = pd.read_csv("CCLE_24Q4/CRISPR_(DepMap_Public_24Q4+Score,_Chronos)_subsetted.csv", index_col = 0)
crispr

crispr = crispr.rename(columns = {'BRAF' : 'BRAF_crispr',
                                  'KRAS' : 'KRAS_crispr',
                                  'MYC' : 'MYC_crispr',
                                  'PIK3CA' : 'PIK3CA_crispr',
                                  'PTEN' : 'PTEN_crispr',
                                  'TP53' : 'TP53_crispr',
                                  'PLA2G4A' : 'PLA2G4A_crispr'})
crispr.head(3)

mutations = pd.read_csv("CCLE_24Q4/Hotspot_Mutations_subsetted.csv", index_col = 0)
mutations

def get_mutation_status(value):
  if value == 0:
    return "WT"
  else:
    return "MUT"

mutations_genes = mutations.iloc[:, 6:]
mutations_genes = mutations_genes.map(get_mutation_status)
mutations = pd.concat([mutations.iloc[:, 0:6], mutations_genes], axis=1)
mutations.head(3)

mutations = mutations.rename(columns = {'BRAF' : 'BRAF_mutation',
                                        'KRAS' : 'KRAS_mutation',
                                        'MYC' : 'MYC_mutation',
                                        'PIK3CA' : 'PIK3CA_mutation',
                                        'PTEN' : 'PTEN_mutation',
                                        'TP53' : 'TP53_mutation'})
mutations.head(3)

cnv = pd.read_csv("CCLE_24Q4/Omics_Absolute_CN_Gene_Public_24Q4_subsetted.csv", index_col = 0)
cnv

def get_cnv_status(value):
  if value == 2.0:
    return "diploid"
  elif value < 2.0:
    return "loss"
  else:
    return "amplification"

cnv_genes = cnv.iloc[:, 6:]
cnv_genes = cnv_genes.map(get_cnv_status)
cnv = pd.concat([cnv.iloc[:, 0:6], cnv_genes], axis=1)
cnv.head(3)

cnv = cnv.rename(columns = {'BRAF' : 'BRAF_cnv',
                            'KRAS' : 'KRAS_cnv',
                            'MYC' : 'MYC_cvn',
                            'PIK3CA' : 'PIK3CA_cnv',
                            'PTEN' : 'PTEN_cnv',
                            'TP53' : 'TP53_cnv',
                            'PLA2G4A' : 'PLA2G4A_cnv'})
cnv.head(3)

drugs = pd.read_csv('CCLE_24Q4/Drug_sensitivity_AUC_(PRISM_Repurposing_Secondary_Screen)_subsetted.csv', index_col = 0)
drugs

for i in drugs.columns:
  if i.upper().startswith('ASPIRIN'):
    print(i)

aspirin = drugs[['ASPIRIN (BRD:BRD-K11433652-001-17-0)']]
aspirin = aspirin.rename(columns = {'ASPIRIN (BRD:BRD-K11433652-001-17-0)' : 'Aspirin'})
aspirin

"""#Merging data"""

crispr['lineage_1'].value_counts()

crispr

common_index = set(crispr.index).intersection(set(mutations.index))
common_index = common_index.intersection(set(cnv.index))
common_index = list(common_index)
len(common_index)

crispr = crispr.loc[common_index]
cnv = cnv.loc[common_index]
mutations = mutations.loc[common_index]

all_data = pd.concat([crispr[['lineage_1', 'PLA2G4A_crispr']],
                      mutations[['BRAF_mutation', 'KRAS_mutation', 'MYC_mutation', 'PIK3CA_mutation', 'PTEN_mutation', 'TP53_mutation']],
                      cnv[['BRAF_cnv', 'PTEN_cnv', 'TP53_cnv', 'PIK3CA_cnv', 'MYC_cvn', 'KRAS_cnv', 'PLA2G4A_cnv']]],
                     axis = "columns")
all_data

all_data

all_data_altered_non_altered = all_data.copy()

all_data_altered_non_altered['BRAF_status'] = all_data_altered_non_altered['BRAF_mutation'] + str('_') +  all_data_altered_non_altered['BRAF_cnv']
braf = []
for i in all_data_altered_non_altered.BRAF_status:
  if 'amplification' in i or 'MUT' in i:
    braf.append('altered')
  else:
    braf.append('non_altered')
all_data_altered_non_altered['BRAF_status'] = braf


all_data_altered_non_altered['KRAS_status'] = all_data_altered_non_altered['KRAS_mutation'] + str('_') +  all_data_altered_non_altered['KRAS_cnv']
kras = []
for i in all_data_altered_non_altered.KRAS_status:
  if 'amplification' in i or 'MUT' in i:
    kras.append('altered')
  else:
    kras.append('non_altered')
all_data_altered_non_altered['KRAS_status'] = kras


all_data_altered_non_altered['MYC_status'] = all_data_altered_non_altered['MYC_mutation'] + str('_') +  all_data_altered_non_altered['MYC_cvn']
myc = []
for i in all_data_altered_non_altered.MYC_status:
  if 'amplification' in i or 'MUT' in i:
    myc.append('altered')
  else:
    myc.append('non_altered')
all_data_altered_non_altered['MYC_status'] = myc


all_data_altered_non_altered['PIK3CA_status'] = all_data_altered_non_altered['PIK3CA_mutation'] + str('_') +  all_data_altered_non_altered['PIK3CA_cnv']
pik3ca = []
for i in all_data_altered_non_altered.PIK3CA_status:
  if 'amplification' in i or 'MUT' in i:
    pik3ca.append('altered')
  else:
    pik3ca.append('non_altered')
all_data_altered_non_altered['PIK3CA_status'] = pik3ca

all_data_altered_non_altered['PTEN_status'] = all_data_altered_non_altered['PTEN_mutation'] + str('_') +  all_data_altered_non_altered['PTEN_cnv']
pten = []
for i in all_data_altered_non_altered.PTEN_status:
  if 'loss' in i or 'MUT' in i:
    pten.append('altered')
  else:
    pten.append('non_altered')
all_data_altered_non_altered['PTEN_status'] = pten


all_data_altered_non_altered['TP53_status'] = all_data_altered_non_altered['TP53_mutation'] + str('_') +  all_data_altered_non_altered['TP53_cnv']
tp53 = []
for i in all_data_altered_non_altered.TP53_status:
  if 'loss' in i or 'MUT' in i:
    tp53.append('altered')
  else:
    tp53.append('non_altered')
all_data_altered_non_altered['TP53_status'] = tp53


all_data_altered_non_altered

all_data_altered_non_altered.to_csv('all_data_CCLE_cPLA2_KO.csv')

def get_altered_plots(data, gene_name):
  sns.set_context("talk", font_scale=1.3)
  sns.boxplot(data = data, x = f'{gene_name}', y = 'PLA2G4A_crispr', color="grey", order = ['non_altered', 'altered'])
  ax = sns.swarmplot(data = data, x = f'{gene_name}', y = 'PLA2G4A_crispr', dodge=True, color="black", order = ['non_altered', 'altered'], size = 2)
  plt.savefig(f'drive/My Drive/Colab Notebooks/GP_lab_data/Nature_cancer_revisions/results_altered_non_altered/{gene_name}_plot.svg', bbox_inches = 'tight')
  plt.show()

status_genes = [x for x in all_data_altered_non_altered.columns if 'status' in x]
for i in status_genes:
  get_altered_plots(data = all_data_altered_non_altered, gene_name = i)

aspirin_common = list(set(common_index).intersection(set(aspirin.index)))
len(aspirin_common)
aspirin = aspirin.loc[aspirin_common]
aspirin

aspirin = pd.concat([aspirin,
                     mutations[['lineage_1', 'BRAF_mutation', 'KRAS_mutation', 'MYC_mutation', 'PIK3CA_mutation', 'PTEN_mutation', 'TP53_mutation']].loc[aspirin_common],
                     cnv[['BRAF_cnv', 'PTEN_cnv', 'TP53_cnv', 'PIK3CA_cnv', 'MYC_cvn', 'KRAS_cnv', 'PLA2G4A_cnv']].loc[aspirin_common]],
                     axis = "columns")
aspirin

"""#Get data for each tissue (Breast/bowel)"""

breast = all_data[all_data['lineage_1'] == 'Breast']
breast

bowel = all_data[all_data['lineage_1'] == 'Bowel']
bowel

aspirin_breast = aspirin[aspirin['lineage_1'] == 'Breast']
aspirin_breast

aspirin_bowel = aspirin[aspirin['lineage_1'] == 'Bowel']
aspirin_bowel

"""#Plots"""

import scipy.stats as stats

def get_cnv_plots(df_tissue, tissue_name, gene_cnv):
  #Get plot
  sns.set_context("talk", font_scale=1.3)
  sns.boxplot(data = df_tissue, x = f'{gene_cnv}', y = 'PLA2G4A_crispr', color="grey", order = ['diploid', 'amplification', 'loss'])
  ax = sns.swarmplot(data = df_tissue, x = f'{gene_cnv}', y = 'PLA2G4A_crispr', dodge=True, color="black", order = ['diploid', 'amplification', 'loss'])
  plt.title(tissue_name)
  plt.savefig(f'results/{tissue_name}_{gene_cnv}_plot.svg', bbox_inches = 'tight')
  plt.show()


  #Get t-test statistics
  group_1 = df_tissue[df_tissue[f'{gene_cnv}'] == 'diploid']['PLA2G4A_crispr']
  group_2 = df_tissue[df_tissue[f'{gene_cnv}'] == 'amplification']['PLA2G4A_crispr']
  group_3 = df_tissue[df_tissue[f'{gene_cnv}'] == 'loss']['PLA2G4A_crispr']

  group_names = ['diploid', 'amplification', 'loss']

  results = []

  for i in range(len(group_names)):
      for j in range(i + 1, len(group_names)):
          group1_name = group_names[i]
          group2_name = group_names[j]

          group1_data = df_tissue[df_tissue[f'{gene_cnv}'] == group1_name]['PLA2G4A_crispr']
          group2_data = df_tissue[df_tissue[f'{gene_cnv}'] == group2_name]['PLA2G4A_crispr']

          t_statistic, p_value = stats.ttest_ind(group1_data, group2_data)

          results.append([group1_name, group2_name, t_statistic, p_value])


  results_df = pd.DataFrame(results, columns=['Group 1', 'Group 2', 'T-statistic', 'P-value'])
  results_df.to_csv(f'results/{tissue_name}_{gene_cnv}_t_test.csv')

bowel.columns

cnv_list = ['BRAF_cnv', 'PTEN_cnv', 'TP53_cnv', 'PIK3CA_cnv', 'MYC_cvn', 'KRAS_cnv']
for gene in cnv_list:
 get_cnv_plots(df_tissue = bowel, tissue_name = 'Bowel', gene_cnv = gene)

for gene in cnv_list:
 get_cnv_plots(df_tissue = breast, tissue_name = 'Breast', gene_cnv = gene)

def get_cnv_plots_grid(df_tissue, tissue_name, cnv_list):
    sns.set_context("talk", font_scale=1.1)

    # Set up the grid
    n_genes = len(cnv_list)
    n_cols = 3
    n_rows = (n_genes + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(18, 6 * n_rows))
    axes = axes.flatten()

    results = []

    for i, gene_cnv in enumerate(cnv_list):
        ax = axes[i]

        # Plot CNV box + swarm
        sns.boxplot(
            data=df_tissue,
            x=gene_cnv,
            y='PLA2G4A_crispr',
            color="lightgrey",
            order=['diploid', 'amplification', 'loss'],
            ax=ax
        )
        sns.swarmplot(
            data=df_tissue,
            x=gene_cnv,
            y='PLA2G4A_crispr',
            dodge=True,
            color="black",
            order=['diploid', 'amplification', 'loss'],
            ax=ax
        )

        ax.set_title(gene_cnv.replace('_cnv', ''))
        ax.set_xlabel('')
        ax.set_ylabel('PLA2G4A_crispr')

        # Compute t-tests
        group_names = ['diploid', 'amplification', 'loss']
        for g1 in range(len(group_names)):
            for g2 in range(g1 + 1, len(group_names)):
                group1_name = group_names[g1]
                group2_name = group_names[g2]
                g1_data = df_tissue[df_tissue[gene_cnv] == group1_name]['PLA2G4A_crispr']
                g2_data = df_tissue[df_tissue[gene_cnv] == group2_name]['PLA2G4A_crispr']
                if len(g1_data) > 1 and len(g2_data) > 1:
                    t_stat, p_val = stats.ttest_ind(g1_data, g2_data, nan_policy='omit')
                    results.append([gene_cnv, group1_name, group2_name, t_stat, p_val])

    # Remove empty subplots
    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])

    plt.suptitle(f"{tissue_name} CNV vs PLA2G4A_crispr", fontsize=20)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(
        f"drive/My Drive/Colab Notebooks/GP_lab_data/Nature_cancer_revisions/results/{tissue_name}_CNV_grid.svg",
        bbox_inches="tight"
    )
    plt.show()

    # Save t-test results
    results_df = pd.DataFrame(results, columns=['Gene', 'Group 1', 'Group 2', 'T-statistic', 'P-value'])
    results_df.to_csv(
        f"results/{tissue_name}_CNV_t_tests.csv",
        index=False
    )

# Example usage
cnv_list = ['BRAF_cnv', 'PTEN_cnv', 'TP53_cnv', 'PIK3CA_cnv', 'MYC_cvn', 'KRAS_cnv']
get_cnv_plots_grid(df_tissue=bowel, tissue_name='Bowel', cnv_list=cnv_list)
get_cnv_plots_grid(df_tissue = breast, tissue_name = 'Breast', cnv_list = cnv_list)

def get_mutation_plots(df_tissue, tissue_name, gene_mutation):
  #Get plot
  sns.set_context("talk", font_scale=1.3)
  sns.boxplot(data = df_tissue, x = f'{gene_mutation}', y = 'PLA2G4A_crispr', color="grey", order = ['WT', 'MUT'])
  ax = sns.swarmplot(data = df_tissue, x = f'{gene_mutation}', y = 'PLA2G4A_crispr', order = ['WT', 'MUT'], dodge=True, color="black")
  plt.title(tissue_name)
  plt.savefig(f'results/{tissue_name}_{gene_mutation}_plot.svg', bbox_inches = 'tight')
  plt.show()


  #Get t-test statistics
  group_1 = df_tissue[df_tissue[f'{gene_mutation}'] == 'WT']['PLA2G4A_crispr']
  group_2 = df_tissue[df_tissue[f'{gene_mutation}'] == 'MUT']['PLA2G4A_crispr']

  group_names = ['WT', 'MUT', 'loss']

  results = []

  for i in range(len(group_names)):
      for j in range(i + 1, len(group_names)):
          group1_name = group_names[i]
          group2_name = group_names[j]

          group1_data = df_tissue[df_tissue[f'{gene_mutation}'] == group1_name]['PLA2G4A_crispr']
          group2_data = df_tissue[df_tissue[f'{gene_mutation}'] == group2_name]['PLA2G4A_crispr']

          t_statistic, p_value = stats.ttest_ind(group1_data, group2_data)

          results.append([group1_name, group2_name, t_statistic, p_value])


  results_df = pd.DataFrame(results, columns=['Group 1', 'Group 2', 'T-statistic', 'P-value'])
  results_df.to_csv(f'results/{tissue_name}_{gene_mutation}_t_test.csv')

mutations.columns

mutations_list = [ 'BRAF_mutation', 'KRAS_mutation', 'MYC_mutation', 'PIK3CA_mutation', 'PTEN_mutation', 'TP53_mutation']
for gene in mutations_list:
 get_mutation_plots(df_tissue = bowel, tissue_name = 'Bowel', gene_mutation = gene)

for gene in mutations_list:
 get_mutation_plots(df_tissue = breast, tissue_name = 'Breast', gene_mutation = gene)

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.stats as stats

def get_mutation_plots_grid(df_tissue, tissue_name, mutation_list):
    sns.set_context("talk", font_scale=1.1)

    # Set up grid (2x3, 3x3, etc. depending on how many genes)
    n_genes = len(mutation_list)
    n_cols = 3
    n_rows = (n_genes + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(18, 6 * n_rows))
    axes = axes.flatten()

    results = []

    for i, gene_mut in enumerate(mutation_list):
        ax = axes[i]

        # Box + swarm plot
        sns.boxplot(
            data=df_tissue,
            x=gene_mut,
            y='PLA2G4A_crispr',
            color="lightgrey",
            order=['WT', 'MUT'],
            ax=ax
        )
        sns.swarmplot(
            data=df_tissue,
            x=gene_mut,
            y='PLA2G4A_crispr',
            color="black",
            order=['WT', 'MUT'],
            dodge=True,
            ax=ax
        )

        ax.set_title(gene_mut.replace('_mutation', ''))
        ax.set_xlabel('')
        ax.set_ylabel('PLA2G4A_crispr')

        # Compute t-test (WT vs MUT only)
        wt = df_tissue[df_tissue[gene_mut] == 'WT']['PLA2G4A_crispr']
        mut = df_tissue[df_tissue[gene_mut] == 'MUT']['PLA2G4A_crispr']

        if len(wt) > 1 and len(mut) > 1:
            t_stat, p_val = stats.ttest_ind(wt, mut, nan_policy='omit')
            results.append([gene_mut, 'WT', 'MUT', t_stat, p_val])

    # Remove unused subplots
    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])

    plt.suptitle(f"{tissue_name} Mutations vs PLA2G4A_crispr", fontsize=20)
    plt.tight_layout(rect=[0, 0, 1, 0.96])

    # Save grid figure
    plt.savefig(
        f"results/{tissue_name}_MUT_grid.svg",
        bbox_inches="tight"
    )
    plt.show()

    # Save t-test results
    results_df = pd.DataFrame(results, columns=['Gene', 'Group 1', 'Group 2', 'T-statistic', 'P-value'])
    results_df.to_csv(
        f"results/{tissue_name}_MUT_t_tests.csv",
        index=False
    )
mutation_list = [ 'BRAF_mutation', 'KRAS_mutation', 'MYC_mutation', 'PIK3CA_mutation', 'PTEN_mutation', 'TP53_mutation']
get_mutation_plots_grid(df_tissue=bowel, tissue_name='Bowel', mutation_list=mutation_list)
get_mutation_plots_grid(df_tissue=bowel, tissue_name='Breast', mutation_list=mutation_list)

def get_cnv_plots_1(df_tissue, tissue_name, gene_cnv):
  #Get plot
  sns.set_context("talk", font_scale=1.3)
  sns.boxplot(data = df_tissue, x = f'{gene_cnv}', y = 'PLA2G4A_crispr', color="grey", order = ['diploid', 'amplification', 'loss'])
  ax = sns.swarmplot(data = df_tissue, x = f'{gene_cnv}', y = 'PLA2G4A_crispr', dodge=True, color="black", order = ['diploid', 'amplification', 'loss'], size = 2)
  plt.title(tissue_name)
  plt.savefig(f'results/{tissue_name}_{gene_cnv}_plot.svg', bbox_inches = 'tight')
  plt.show()


  #Get t-test statistics
  group_1 = df_tissue[df_tissue[f'{gene_cnv}'] == 'diploid']['PLA2G4A_crispr']
  group_2 = df_tissue[df_tissue[f'{gene_cnv}'] == 'amplification']['PLA2G4A_crispr']
  group_3 = df_tissue[df_tissue[f'{gene_cnv}'] == 'loss']['PLA2G4A_crispr']

  group_names = ['diploid', 'amplification', 'loss']

  results = []

  for i in range(len(group_names)):
      for j in range(i + 1, len(group_names)):
          group1_name = group_names[i]
          group2_name = group_names[j]

          group1_data = df_tissue[df_tissue[f'{gene_cnv}'] == group1_name]['PLA2G4A_crispr']
          group2_data = df_tissue[df_tissue[f'{gene_cnv}'] == group2_name]['PLA2G4A_crispr']

          t_statistic, p_value = stats.ttest_ind(group1_data, group2_data)

          results.append([group1_name, group2_name, t_statistic, p_value])


  results_df = pd.DataFrame(results, columns=['Group 1', 'Group 2', 'T-statistic', 'P-value'])
  results_df.to_csv(f'results/{tissue_name}_{gene_cnv}_t_test.csv')

def get_mutation_plots_1(df_tissue, tissue_name, gene_mutation):
  #Get plot
  sns.set_context("talk", font_scale=1.3)
  sns.boxplot(data = df_tissue, x = f'{gene_mutation}', y = 'PLA2G4A_crispr', color="grey", order = ['WT', 'MUT'])
  ax = sns.swarmplot(data = df_tissue, x = f'{gene_mutation}', y = 'PLA2G4A_crispr', order = ['WT', 'MUT'], dodge=True, color="black", size = 2)
  plt.title(tissue_name)
  plt.savefig(f'results/{tissue_name}_{gene_mutation}_plot.svg', bbox_inches = 'tight')
  plt.show()


  #Get t-test statistics
  group_1 = df_tissue[df_tissue[f'{gene_mutation}'] == 'WT']['PLA2G4A_crispr']
  group_2 = df_tissue[df_tissue[f'{gene_mutation}'] == 'MUT']['PLA2G4A_crispr']

  group_names = ['WT', 'MUT', 'loss']

  results = []

  for i in range(len(group_names)):
      for j in range(i + 1, len(group_names)):
          group1_name = group_names[i]
          group2_name = group_names[j]

          group1_data = df_tissue[df_tissue[f'{gene_mutation}'] == group1_name]['PLA2G4A_crispr']
          group2_data = df_tissue[df_tissue[f'{gene_mutation}'] == group2_name]['PLA2G4A_crispr']

          t_statistic, p_value = stats.ttest_ind(group1_data, group2_data)

          results.append([group1_name, group2_name, t_statistic, p_value])


  results_df = pd.DataFrame(results, columns=['Group 1', 'Group 2', 'T-statistic', 'P-value'])
  results_df.to_csv(f'results/{tissue_name}_{gene_mutation}_t_test.csv')

for gene in cnv_list:
 get_cnv_plots_1(df_tissue = all_data, tissue_name = 'All tumours', gene_cnv = gene)

for gene in mutations_list:
 get_mutation_plots_1(df_tissue = all_data, tissue_name = 'All tumours', gene_mutation = gene)

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.stats as stats

def get_cnv_plots_grid_1(df_tissue, tissue_name, cnv_list):
    sns.set_context("talk", font_scale=1.1)

    n_genes = len(cnv_list)
    n_cols = 3
    n_rows = (n_genes + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(18, 6 * n_rows))
    axes = axes.flatten()

    results = []

    for i, gene_cnv in enumerate(cnv_list):
        ax = axes[i]

        # Box + swarm
        sns.boxplot(
            data=df_tissue,
            x=gene_cnv,
            y='PLA2G4A_crispr',
            color="lightgrey",
            order=['diploid', 'amplification', 'loss'],
            ax=ax
        )
        sns.swarmplot(
            data=df_tissue,
            x=gene_cnv,
            y='PLA2G4A_crispr',
            color="black",
            order=['diploid', 'amplification', 'loss'],
            dodge=True,
            size=2,
            ax=ax
        )

        ax.set_title(gene_cnv.replace('_cnv', ''))
        ax.set_xlabel('')
        ax.set_ylabel('PLA2G4A_crispr')

        # Compute t-tests
        group_names = ['diploid', 'amplification', 'loss']
        for g1 in range(len(group_names)):
            for g2 in range(g1 + 1, len(group_names)):
                g1_name = group_names[g1]
                g2_name = group_names[g2]
                g1_data = df_tissue[df_tissue[gene_cnv] == g1_name]['PLA2G4A_crispr']
                g2_data = df_tissue[df_tissue[gene_cnv] == g2_name]['PLA2G4A_crispr']
                if len(g1_data) > 1 and len(g2_data) > 1:
                    t_stat, p_val = stats.ttest_ind(g1_data, g2_data, nan_policy='omit')
                    results.append([gene_cnv, g1_name, g2_name, t_stat, p_val])

    # Remove empty plots
    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])

    plt.suptitle(f"{tissue_name} CNV vs PLA2G4A_crispr", fontsize=20)
    plt.tight_layout(rect=[0, 0, 1, 0.96])

    plt.savefig(
        f"results/{tissue_name}_CNV_grid_1.svg",
        bbox_inches="tight"
    )
    plt.show()

    results_df = pd.DataFrame(results, columns=['Gene', 'Group 1', 'Group 2', 'T-statistic', 'P-value'])
    results_df.to_csv(
        f"results/{tissue_name}_CNV_t_tests_1.csv",
        index=False
    )

def get_mutation_plots_grid_1(df_tissue, tissue_name, mutation_list):
    sns.set_context("talk", font_scale=1.1)

    n_genes = len(mutation_list)
    n_cols = 3
    n_rows = (n_genes + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(18, 6 * n_rows))
    axes = axes.flatten()

    results = []

    for i, gene_mut in enumerate(mutation_list):
        ax = axes[i]

        sns.boxplot(
            data=df_tissue,
            x=gene_mut,
            y='PLA2G4A_crispr',
            color="lightgrey",
            order=['WT', 'MUT'],
            ax=ax
        )
        sns.swarmplot(
            data=df_tissue,
            x=gene_mut,
            y='PLA2G4A_crispr',
            color="black",
            order=['WT', 'MUT'],
            dodge=True,
            size=2,
            ax=ax
        )

        ax.set_title(gene_mut.replace('_mutation', ''))
        ax.set_xlabel('')
        ax.set_ylabel('PLA2G4A_crispr')

        # Compute t-test (WT vs MUT)
        wt = df_tissue[df_tissue[gene_mut] == 'WT']['PLA2G4A_crispr']
        mut = df_tissue[df_tissue[gene_mut] == 'MUT']['PLA2G4A_crispr']

        if len(wt) > 1 and len(mut) > 1:
            t_stat, p_val = stats.ttest_ind(wt, mut, nan_policy='omit')
            results.append([gene_mut, 'WT', 'MUT', t_stat, p_val])

    # Remove empty plots
    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])

    plt.suptitle(f"{tissue_name} Mutations vs PLA2G4A_crispr", fontsize=20)
    plt.tight_layout(rect=[0, 0, 1, 0.96])

    plt.savefig(
        f"results/{tissue_name}_MUT_grid_1.svg",
        bbox_inches="tight"
    )
    plt.show()

    results_df = pd.DataFrame(results, columns=['Gene', 'Group 1', 'Group 2', 'T-statistic', 'P-value'])
    results_df.to_csv(
        f"results/{tissue_name}_MUT_t_tests_1.csv",
        index=False
    )

get_cnv_plots_grid_1(df_tissue=bowel, tissue_name='All tumours', cnv_list=cnv_list)
get_mutation_plots_grid_1(df_tissue=bowel, tissue_name='All tumours', mutation_list=mutation_list)

plt.title('All tumours')
sns.histplot(all_data.PLA2G4A_crispr)

def get_cnv_plots_aspirin(df_tissue, tissue_name, gene_cnv):
  #Get plot
  sns.set_context("talk", font_scale=1.3)
  sns.boxplot(data = df_tissue, x = f'{gene_cnv}', y = 'Aspirin', color="grey", order = ['diploid', 'amplification', 'loss'])
  ax = sns.swarmplot(data = df_tissue, x = f'{gene_cnv}', y = 'Aspirin', dodge=True, color="black", order = ['diploid', 'amplification', 'loss'])
  plt.title(tissue_name)
  plt.savefig(f'results/{tissue_name}_{gene_cnv}_plot.svg', bbox_inches = 'tight')
  plt.show()


  #Get t-test statistics
  group_1 = df_tissue[df_tissue[f'{gene_cnv}'] == 'diploid']['Aspirin']
  group_2 = df_tissue[df_tissue[f'{gene_cnv}'] == 'amplification']['Aspirin']
  group_3 = df_tissue[df_tissue[f'{gene_cnv}'] == 'loss']['Aspirin']

  group_names = ['diploid', 'amplification', 'loss']

  results = []

  for i in range(len(group_names)):
      for j in range(i + 1, len(group_names)):
          group1_name = group_names[i]
          group2_name = group_names[j]

          group1_data = df_tissue[df_tissue[f'{gene_cnv}'] == group1_name]['Aspirin']
          group2_data = df_tissue[df_tissue[f'{gene_cnv}'] == group2_name]['Aspirin']

          t_statistic, p_value = stats.ttest_ind(group1_data, group2_data)

          results.append([group1_name, group2_name, t_statistic, p_value])


  results_df = pd.DataFrame(results, columns=['Group 1', 'Group 2', 'T-statistic', 'P-value'])
  results_df.to_csv(f'results/{tissue_name}_{gene_cnv}_t_test.csv')

for gene in cnv_list:
 get_cnv_plots_aspirin(df_tissue = aspirin_breast, tissue_name = 'Aspirin_Breast', gene_cnv = gene)

for gene in cnv_list:
 get_cnv_plots_aspirin(df_tissue = aspirin_bowel, tissue_name = 'Aspirin_Bowel', gene_cnv = gene)

def get_cnv_plots_aspirin_grid(df_tissue, tissue_name, cnv_list):
    sns.set_context("talk", font_scale=1.1)

    n_genes = len(cnv_list)
    n_cols = 3
    n_rows = (n_genes + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(18, 6 * n_rows))
    axes = axes.flatten()

    results = []

    for i, gene_cnv in enumerate(cnv_list):
        ax = axes[i]

        # Box + swarm plot for Aspirin
        sns.boxplot(
            data=df_tissue,
            x=gene_cnv,
            y='Aspirin',
            color="lightgrey",
            order=['diploid', 'amplification', 'loss'],
            ax=ax
        )
        sns.swarmplot(
            data=df_tissue,
            x=gene_cnv,
            y='Aspirin',
            color="black",
            order=['diploid', 'amplification', 'loss'],
            dodge=True,
            ax=ax
        )

        ax.set_title(gene_cnv.replace('_cnv', ''))
        ax.set_xlabel('')
        ax.set_ylabel('Aspirin')

        # Compute t-tests
        group_names = ['diploid', 'amplification', 'loss']
        for g1 in range(len(group_names)):
            for g2 in range(g1 + 1, len(group_names)):
                g1_name = group_names[g1]
                g2_name = group_names[g2]
                g1_data = df_tissue[df_tissue[gene_cnv] == g1_name]['Aspirin']
                g2_data = df_tissue[df_tissue[gene_cnv] == g2_name]['Aspirin']
                if len(g1_data) > 1 and len(g2_data) > 1:
                    t_stat, p_val = stats.ttest_ind(g1_data, g2_data, nan_policy='omit')
                    results.append([gene_cnv, g1_name, g2_name, t_stat, p_val])

    # Remove unused subplots
    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])

    plt.suptitle(f"{tissue_name} CNV vs Aspirin response", fontsize=20)
    plt.tight_layout(rect=[0, 0, 1, 0.96])

    # Save figure
    plt.savefig(
        f"results/{tissue_name}_CNV_Aspirin_grid.svg",
        bbox_inches="tight"
    )
    plt.show()

    # Save t-test results
    results_df = pd.DataFrame(results, columns=['Gene', 'Group 1', 'Group 2', 'T-statistic', 'P-value'])
    results_df.to_csv(
        f"results/{tissue_name}_CNV_Aspirin_t_tests.csv",
        index=False
    )

get_cnv_plots_aspirin_grid(df_tissue = aspirin_breast, tissue_name = 'Aspirin_Breast', cnv_list=cnv_list)
get_cnv_plots_aspirin_grid(df_tissue = aspirin_bowel, tissue_name = 'Aspirin_Bowel', cnv_list=cnv_list)

def get_mutation_plots_aspirin(df_tissue, tissue_name, gene_mutation):
  #Get plot
  sns.set_context("talk", font_scale=1.3)
  sns.boxplot(data = df_tissue, x = f'{gene_mutation}', y = 'Aspirin', color="grey", order = ['WT', 'MUT'])
  ax = sns.swarmplot(data = df_tissue, x = f'{gene_mutation}', y = 'Aspirin', order = ['WT', 'MUT'], dodge=True, color="black")
  plt.title(tissue_name)
  plt.savefig(f'results/{tissue_name}_{gene_mutation}_plot.svg', bbox_inches = 'tight')
  plt.show()


  #Get t-test statistics
  group_1 = df_tissue[df_tissue[f'{gene_mutation}'] == 'WT']['Aspirin']
  group_2 = df_tissue[df_tissue[f'{gene_mutation}'] == 'MUT']['Aspirin']

  group_names = ['WT', 'MUT', 'loss']

  results = []

  for i in range(len(group_names)):
      for j in range(i + 1, len(group_names)):
          group1_name = group_names[i]
          group2_name = group_names[j]

          group1_data = df_tissue[df_tissue[f'{gene_mutation}'] == group1_name]['Aspirin']
          group2_data = df_tissue[df_tissue[f'{gene_mutation}'] == group2_name]['Aspirin']

          t_statistic, p_value = stats.ttest_ind(group1_data, group2_data)

          results.append([group1_name, group2_name, t_statistic, p_value])


  results_df = pd.DataFrame(results, columns=['Group 1', 'Group 2', 'T-statistic', 'P-value'])
  results_df.to_csv(f'results/{tissue_name}_{gene_mutation}_t_test.csv')

for gene in mutations_list:
 get_mutation_plots_aspirin(df_tissue = aspirin_breast, tissue_name = 'Aspirin_Breast', gene_mutation = gene)

for gene in mutations_list:
 get_mutation_plots_aspirin(df_tissue = aspirin_bowel, tissue_name = 'Aspirin_Bowel', gene_mutation = gene)

def get_mutation_plots_aspirin_grid(df_tissue, tissue_name, mutation_list):
    sns.set_context("talk", font_scale=1.1)

    n_genes = len(mutation_list)
    n_cols = 3
    n_rows = (n_genes + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(18, 6 * n_rows))
    axes = axes.flatten()

    results = []

    for i, gene_mut in enumerate(mutation_list):
        ax = axes[i]

        # Box + swarm plot
        sns.boxplot(
            data=df_tissue,
            x=gene_mut,
            y='Aspirin',
            color="lightgrey",
            order=['WT', 'MUT'],
            ax=ax
        )
        sns.swarmplot(
            data=df_tissue,
            x=gene_mut,
            y='Aspirin',
            color="black",
            order=['WT', 'MUT'],
            dodge=True,
            ax=ax
        )

        ax.set_title(gene_mut.replace('_mutation', ''))
        ax.set_xlabel('')
        ax.set_ylabel('Aspirin')

        # t-test WT vs MUT
        wt = df_tissue[df_tissue[gene_mut] == 'WT']['Aspirin']
        mut = df_tissue[df_tissue[gene_mut] == 'MUT']['Aspirin']

        if len(wt) > 1 and len(mut) > 1:
            t_stat, p_val = stats.ttest_ind(wt, mut, nan_policy='omit')
            results.append([gene_mut, 'WT', 'MUT', t_stat, p_val])

    # Remove unused subplots
    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])

    plt.suptitle(f"{tissue_name} Mutations vs Aspirin response", fontsize=20)
    plt.tight_layout(rect=[0, 0, 1, 0.96])

    # Save figure
    plt.savefig(
        f"results/{tissue_name}_MUT_Aspirin_grid.svg",
        bbox_inches="tight"
    )
    plt.show()

    # Save t-test results
    results_df = pd.DataFrame(results, columns=['Gene', 'Group 1', 'Group 2', 'T-statistic', 'P-value'])
    results_df.to_csv(
        f"results/{tissue_name}_MUT_Aspirin_t_tests.csv",
        index=False
    )


get_mutation_plots_aspirin_grid(df_tissue = aspirin_bowel, tissue_name = 'Aspirin_Bowel', mutation_list=mutation_list)
get_mutation_plots_aspirin_grid(df_tissue = aspirin_breast, tissue_name = 'Aspirin_Breast', mutation_list=mutation_list)

def get_cnv_plots_aspirin_1(df_tissue, tissue_name, gene_cnv):
  #Get plot
  sns.set_context("talk", font_scale=1.3)
  sns.boxplot(data = df_tissue, x = f'{gene_cnv}', y = 'Aspirin', color="grey", order = ['diploid', 'amplification', 'loss'])
  ax = sns.swarmplot(data = df_tissue, x = f'{gene_cnv}', y = 'Aspirin', dodge=True, color="black", order = ['diploid', 'amplification', 'loss'], size = 2)
  plt.title(tissue_name)
  plt.savefig(f'results/{tissue_name}_{gene_cnv}_plot.svg', bbox_inches = 'tight')
  plt.show()


  #Get t-test statistics
  group_1 = df_tissue[df_tissue[f'{gene_cnv}'] == 'diploid']['Aspirin']
  group_2 = df_tissue[df_tissue[f'{gene_cnv}'] == 'amplification']['Aspirin']
  group_3 = df_tissue[df_tissue[f'{gene_cnv}'] == 'loss']['Aspirin']

  group_names = ['diploid', 'amplification', 'loss']

  results = []

  for i in range(len(group_names)):
      for j in range(i + 1, len(group_names)):
          group1_name = group_names[i]
          group2_name = group_names[j]

          group1_data = df_tissue[df_tissue[f'{gene_cnv}'] == group1_name]['Aspirin']
          group2_data = df_tissue[df_tissue[f'{gene_cnv}'] == group2_name]['Aspirin']

          t_statistic, p_value = stats.ttest_ind(group1_data, group2_data)

          results.append([group1_name, group2_name, t_statistic, p_value])


  results_df = pd.DataFrame(results, columns=['Group 1', 'Group 2', 'T-statistic', 'P-value'])
  results_df.to_csv(f'drive/My Drive/Colab Notebooks/GP_lab_data/Nature_cancer_revisions/results/{tissue_name}_{gene_cnv}_t_test.csv')

def get_mutation_plots_aspirin_1(df_tissue, tissue_name, gene_mutation):
  #Get plot
  sns.set_context("talk", font_scale=1.3)
  sns.boxplot(data = df_tissue, x = f'{gene_mutation}', y = 'Aspirin', color="grey", order = ['WT', 'MUT'])
  ax = sns.swarmplot(data = df_tissue, x = f'{gene_mutation}', y = 'Aspirin', order = ['WT', 'MUT'], dodge=True, color="black", size = 2)
  plt.title(tissue_name)
  plt.savefig(f'results/{tissue_name}_{gene_mutation}_plot.svg', bbox_inches = 'tight')
  plt.show()


  #Get t-test statistics
  group_1 = df_tissue[df_tissue[f'{gene_mutation}'] == 'WT']['Aspirin']
  group_2 = df_tissue[df_tissue[f'{gene_mutation}'] == 'MUT']['Aspirin']

  group_names = ['WT', 'MUT', 'loss']

  results = []

  for i in range(len(group_names)):
      for j in range(i + 1, len(group_names)):
          group1_name = group_names[i]
          group2_name = group_names[j]

          group1_data = df_tissue[df_tissue[f'{gene_mutation}'] == group1_name]['Aspirin']
          group2_data = df_tissue[df_tissue[f'{gene_mutation}'] == group2_name]['Aspirin']

          t_statistic, p_value = stats.ttest_ind(group1_data, group2_data)

          results.append([group1_name, group2_name, t_statistic, p_value])


  results_df = pd.DataFrame(results, columns=['Group 1', 'Group 2', 'T-statistic', 'P-value'])
  results_df.to_csv(f'results/{tissue_name}_{gene_mutation}_t_test.csv')

for gene in cnv_list:
 get_cnv_plots_aspirin_1(df_tissue = aspirin, tissue_name = 'Aspirin All tumours', gene_cnv = gene)

for gene in mutations_list:
 get_mutation_plots_aspirin_1(df_tissue = aspirin, tissue_name = 'Aspirin All tumours', gene_mutation = gene)

def get_cnv_plots_aspirin_grid_1(df_tissue, tissue_name, cnv_list):
    sns.set_context("talk", font_scale=1.1)

    n_genes = len(cnv_list)
    n_cols = 3
    n_rows = (n_genes + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(18, 6 * n_rows))
    axes = axes.flatten()
    results = []

    for i, gene_cnv in enumerate(cnv_list):
        ax = axes[i]

        sns.boxplot(
            data=df_tissue,
            x=gene_cnv,
            y='Aspirin',
            color="lightgrey",
            order=['diploid', 'amplification', 'loss'],
            ax=ax
        )
        sns.swarmplot(
            data=df_tissue,
            x=gene_cnv,
            y='Aspirin',
            color="black",
            order=['diploid', 'amplification', 'loss'],
            dodge=True,
            size=2,
            ax=ax
        )

        ax.set_title(gene_cnv.replace('_cnv', ''))
        ax.set_xlabel('')
        ax.set_ylabel('Aspirin')

        # t-tests
        group_names = ['diploid', 'amplification', 'loss']
        for g1 in range(len(group_names)):
            for g2 in range(g1 + 1, len(group_names)):
                g1_name = group_names[g1]
                g2_name = group_names[g2]
                g1_data = df_tissue[df_tissue[gene_cnv] == g1_name]['Aspirin']
                g2_data = df_tissue[df_tissue[gene_cnv] == g2_name]['Aspirin']
                if len(g1_data) > 1 and len(g2_data) > 1:
                    t_stat, p_val = stats.ttest_ind(g1_data, g2_data, nan_policy='omit')
                    results.append([gene_cnv, g1_name, g2_name, t_stat, p_val])

    # Remove empty subplots
    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])

    plt.suptitle(f"{tissue_name} CNV vs Aspirin response", fontsize=20)
    plt.tight_layout(rect=[0, 0, 1, 0.96])

    plt.savefig(
        f"results/{tissue_name}_CNV_Aspirin_grid_1.svg",
        bbox_inches="tight"
    )
    plt.show()

    results_df = pd.DataFrame(results, columns=['Gene', 'Group 1', 'Group 2', 'T-statistic', 'P-value'])
    results_df.to_csv(
        f"results/{tissue_name}_CNV_Aspirin_t_tests_1.csv",
        index=False
    )

def get_mutation_plots_aspirin_grid_1(df_tissue, tissue_name, mutation_list):
    sns.set_context("talk", font_scale=1.1)

    n_genes = len(mutation_list)
    n_cols = 3
    n_rows = (n_genes + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(18, 6 * n_rows))
    axes = axes.flatten()
    results = []

    for i, gene_mut in enumerate(mutation_list):
        ax = axes[i]

        sns.boxplot(
            data=df_tissue,
            x=gene_mut,
            y='Aspirin',
            color="lightgrey",
            order=['WT', 'MUT'],
            ax=ax
        )
        sns.swarmplot(
            data=df_tissue,
            x=gene_mut,
            y='Aspirin',
            color="black",
            order=['WT', 'MUT'],
            dodge=True,
            size=2,
            ax=ax
        )

        ax.set_title(gene_mut.replace('_mutation', ''))
        ax.set_xlabel('')
        ax.set_ylabel('Aspirin')

        # t-test WT vs MUT
        wt = df_tissue[df_tissue[gene_mut] == 'WT']['Aspirin']
        mut = df_tissue[df_tissue[gene_mut] == 'MUT']['Aspirin']
        if len(wt) > 1 and len(mut) > 1:
            t_stat, p_val = stats.ttest_ind(wt, mut, nan_policy='omit')
            results.append([gene_mut, 'WT', 'MUT', t_stat, p_val])

    # Remove empty subplots
    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])

    plt.suptitle(f"{tissue_name} Mutations vs Aspirin response", fontsize=20)
    plt.tight_layout(rect=[0, 0, 1, 0.96])

    plt.savefig(
        f"results/{tissue_name}_MUT_Aspirin_grid_1.svg",
        bbox_inches="tight"
    )
    plt.show()

    results_df = pd.DataFrame(results, columns=['Gene', 'Group 1', 'Group 2', 'T-statistic', 'P-value'])
    results_df.to_csv(
        f"results/{tissue_name}_MUT_Aspirin_t_tests_1.csv",
        index=False
    )


get_cnv_plots_aspirin_grid_1(df_tissue = aspirin, tissue_name = 'Aspirin All tumours', cnv_list = cnv_list)
get_mutation_plots_aspirin_grid_1(df_tissue = aspirin, tissue_name = 'Aspirin All tumours', mutation_list = mutation_list)