import pandas as pd
import scipy.stats as stats
from more_itertools import locate
import numpy as np
import cobra
from cobra import Reaction
from cobra.flux_analysis import flux_variability_analysis

def fdr(p_vals):

    p_vals = np.array(p_vals,dtype = 'float64')
    ranked_p_values = stats.rankdata(p_vals)
    fdr = p_vals * len(p_vals) / ranked_p_values
    fdr[fdr > 1] = 1

    return fdr

metabolome_biolog_map = pd.read_csv('Calbicans_LRhamnosus_EpithelialCells_Interaction/data/metabolome_biolog_map.csv')

pm1_1 = pd.read_csv('Calbicans_LRhamnosus_EpithelialCells_Interaction/data/calb_biolog/rep1/PM1/PM1_SC5314.csv')
pm2_1 = pd.read_csv('Calbicans_LRhamnosus_EpithelialCells_Interaction/data/calb_biolog/rep1/PM2/PM2_SC5314.csv')
pm3_1 = pd.read_csv('Calbicans_LRhamnosus_EpithelialCells_Interaction/data/calb_biolog/rep1/PM3/PM3_SC5314.csv')
pm4_1 = pd.read_csv('Calbicans_LRhamnosus_EpithelialCells_Interaction/data/calb_biolog/rep1/PM4/PM4_SC5314.csv')

pm1_2 = pd.read_csv('Calbicans_LRhamnosus_EpithelialCells_Interaction/data/calb_biolog/rep2/PM1/PM1_SC5314.csv')
pm2_2 = pd.read_csv('Calbicans_LRhamnosus_EpithelialCells_Interaction/data/calb_biolog/rep2/PM2/PM2_SC5314.csv')
pm3_2 = pd.read_csv('Calbicans_LRhamnosus_EpithelialCells_Interaction/data/calb_biolog/rep2/PM3/PM3_SC5314.csv')
pm4_2 = pd.read_csv('Calbicans_LRhamnosus_EpithelialCells_Interaction/data/calb_biolog/rep2/PM4/PM4_SC5314.csv')

#Table of final growth signals for ANOVA
#pm1
pm1_data = pd.DataFrame(columns=['nc','substrate'],index=['rep1','rep2'])
pm1_data['nc'][0] = pm1_1['  A01'][pm1_1.shape[0]-1]
pm1_data['nc'][1] = pm1_2['  A01'][pm1_2.shape[0]-1]

#pm2
pm2_data = pd.DataFrame(columns=['nc','substrate'],index=['rep1','rep2'])
pm2_data['nc'][0] = pm2_1['  A01'][pm2_1.shape[0]-1]
pm2_data['nc'][1] = pm2_2['  A01'][pm2_2.shape[0]-1]

#pm3
pm3_data = pd.DataFrame(columns=['nc','substrate'],index=['rep1','rep2'])
pm3_data['nc'][0] = pm3_1['  A01'][pm3_1.shape[0]-1]
pm3_data['nc'][1] = pm3_2['  A01'][pm3_2.shape[0]-1]

#pm4_p
pm4_p_data = pd.DataFrame(columns=['nc','substrate'],index=['rep1','rep2'])
pm4_p_data['nc'][0] = pm4_1['  A01'][pm4_1.shape[0]-1]
pm4_p_data['nc'][1] = pm4_2['  A01'][pm4_2.shape[0]-1]

#pm4_s
pm4_s_data = pd.DataFrame(columns=['nc','substrate'],index=['rep1','rep2'])
pm4_s_data['nc'][0] = pm4_1['  F01'][pm4_1.shape[0]-1]
pm4_s_data['nc'][1] = pm4_2['  F01'][pm4_2.shape[0]-1]

source = list()
metabolite = list()
rel_growths = list()
pvalues = list()#anova
wilcoxes = list()
ttests = list()
for i in range(metabolome_biolog_map.shape[0]):

    if not pd.isna(metabolome_biolog_map['biolog_plate'][i]):

        if metabolome_biolog_map['biolog_plate'][i]=='PM01':

            source.append('carbon')
            metabolite.append(metabolome_biolog_map['metabolome_met'][i])

            pm1_data['substrate'][0] = pm1_1['  ' + metabolome_biolog_map['biolog_well_id'][i]][pm1_1.shape[0]-1]
            pm1_data['substrate'][1] = pm1_2['  ' + metabolome_biolog_map['biolog_well_id'][i]][pm1_2.shape[0]-1]

            rel_growth = np.log2(pm1_data.sum()[1]/pm1_data.sum()[0])
            rel_growths.append(rel_growth)

            #if rel_growth < 1:

            fvalue, pvalue = stats.f_oneway(pm1_data['nc'], pm1_data['substrate'])
            wilcox = stats.wilcoxon(pm1_data['nc'], pm1_data['substrate'])
            ttest = stats.ttest_ind(pm1_data['nc'], pm1_data['substrate'])
            pvalues.append(pvalue)
            wilcoxes.append(wilcox[1])
            ttests.append(ttest[1])
            #else:
            #    pvalues.append('')

        if metabolome_biolog_map['biolog_plate'][i]=='PM02':

            source.append('carbon')
            metabolite.append(metabolome_biolog_map['metabolome_met'][i])

            pm2_data['substrate'][0] = pm2_1['  ' + metabolome_biolog_map['biolog_well_id'][i]][pm2_1.shape[0]-1]
            pm2_data['substrate'][1] = pm2_2['  ' + metabolome_biolog_map['biolog_well_id'][i]][pm2_2.shape[0]-1]

            rel_growth = np.log2(pm2_data.sum()[1]/pm2_data.sum()[0])
            rel_growths.append(rel_growth)

            #if rel_growth < 1:

            fvalue, pvalue = stats.f_oneway(pm2_data['nc'], pm2_data['substrate'])
            wilcox = stats.wilcoxon(pm2_data['nc'], pm2_data['substrate'])
            ttest = stats.ttest_ind(pm2_data['nc'], pm2_data['substrate'])
            pvalues.append(pvalue)
            wilcoxes.append(wilcox[1])
            ttests.append(ttest[1])
            #else:
            #    pvalues.append('')

        if metabolome_biolog_map['biolog_plate'][i]=='PM03':

            source.append('nitrogen')
            metabolite.append(metabolome_biolog_map['metabolome_met'][i])

            pm3_data['substrate'][0] = pm3_1['  ' + metabolome_biolog_map['biolog_well_id'][i]][pm3_1.shape[0]-1]
            pm3_data['substrate'][1] = pm3_2['  ' + metabolome_biolog_map['biolog_well_id'][i]][pm3_2.shape[0]-1]

            rel_growth = np.log2(pm3_data.sum()[1]/pm3_data.sum()[0])
            rel_growths.append(rel_growth)

            #if rel_growth < 1:

            fvalue, pvalue = stats.f_oneway(pm3_data['nc'], pm3_data['substrate'])
            wilcoxon = stats.wilcoxon(pm3_data['nc'], pm3_data['substrate'])
            ttest = stats.ttest_ind(pm3_data['nc'], pm3_data['substrate'])
            pvalues.append(pvalue)
            wilcoxes.append(wilcox[1])
            ttests.append(ttest[1])
            #else:
            #    pvalues.append('')

        if metabolome_biolog_map['biolog_plate'][i]=='PM04' and metabolome_biolog_map['biolog_well_id'][i][0] in 'ABCDE':

            source.append('phosphor')
            metabolite.append(metabolome_biolog_map['metabolome_met'][i])

            pm4_p_data['substrate'][0] = pm4_1['  ' + metabolome_biolog_map['biolog_well_id'][i]][pm4_1.shape[0]-1]
            pm4_p_data['substrate'][1] = pm4_2['  ' + metabolome_biolog_map['biolog_well_id'][i]][pm4_2.shape[0]-1]

            rel_growth = np.log2(pm4_p_data.sum()[1]/pm4_p_data.sum()[0])
            rel_growths.append(rel_growth)

            #if rel_growth < 1:

            fvalue, pvalue = stats.f_oneway(pm4_p_data['nc'], pm4_p_data['substrate'])
            wilcox = stats.wilcoxon(pm4_p_data['nc'], pm4_p_data['substrate'])
            ttest = stats.ttest_ind(pm4_p_data['nc'], pm4_p_data['substrate'])
            pvalues.append(pvalue)
            wilcoxes.append(wilcox[1])
            ttests.append(ttest[1])
            #else:
            #    pvalues.append('')

        if metabolome_biolog_map['biolog_plate'][i]=='PM04' and metabolome_biolog_map['biolog_well_id'][i][0] in 'FGH':

            source.append('sulfur')
            metabolite.append(metabolome_biolog_map['metabolome_met'][i])

            pm4_s_data['substrate'][0] = pm4_1['  ' + metabolome_biolog_map['biolog_well_id'][i]][pm4_1.shape[0]-1]
            pm4_s_data['substrate'][1] = pm4_2['  ' + metabolome_biolog_map['biolog_well_id'][i]][pm4_2.shape[0]-1]

            rel_growth = np.log2(pm4_s_data.sum()[1]/pm4_s_data.sum()[0])
            rel_growths.append(rel_growth)

            #if rel_growth < 1:

            fvalue, pvalue = stats.f_oneway(pm4_s_data['nc'], pm4_s_data['substrate'])
            wilcox = stats.wilcoxon(pm4_s_data['nc'], pm4_s_data['substrate'])
            ttest = stats.ttest_ind(pm4_s_data['nc'], pm4_s_data['substrate'])
            pvalues.append(pvalue)
            wilcoxes.append(wilcox[1])
            ttests.append(ttest[1])
            #else:
            #    pvalues.append('')


unq_metabolites = list(set(metabolite))
anova_c = ['']*len(unq_metabolites)
anova_n = ['']*len(unq_metabolites)
anova_p = ['']*len(unq_metabolites)
anova_s = ['']*len(unq_metabolites)
ttest_c = ['']*len(unq_metabolites)
ttest_n = ['']*len(unq_metabolites)
ttest_p = ['']*len(unq_metabolites)
ttest_s = ['']*len(unq_metabolites)
wilcox_c = ['']*len(unq_metabolites)
wilcox_n = ['']*len(unq_metabolites)
wilcox_p = ['']*len(unq_metabolites)
wilcox_s = ['']*len(unq_metabolites)
rel_growths_c = ['']*len(unq_metabolites)
rel_growths_n = ['']*len(unq_metabolites)
rel_growths_p = ['']*len(unq_metabolites)
rel_growths_s = ['']*len(unq_metabolites)

for i in range(len(metabolite)):
    loc = unq_metabolites.index(metabolite[i])
    if source[i]=='carbon':
        anova_c[loc] = pvalues[i]
        wilcox_c[loc] = wilcoxes[i]
        ttest_c[loc] = ttests[i]
        rel_growths_c[loc] = rel_growths[i]
    if source[i]=='nitrogen':
        anova_n[loc] = pvalues[i]
        wilcox_n[loc] = wilcoxes[i]
        ttest_n[loc] = ttests[i]
        rel_growths_n[loc] = rel_growths[i]
    if source[i]=='phosphor':
        anova_p[loc] = pvalues[i]
        wilcox_p[loc] = wilcoxes[i]
        ttest_p[loc] = ttests[i]
        rel_growths_p[loc] = rel_growths[i]
    if source[i]=='sulfur':
        anova_s[loc] = pvalues[i]
        wilcox_s[loc] = wilcoxes[i]
        ttest_s[loc] = ttests[i]
        rel_growths_s[loc] = rel_growths[i]
        
#pvalue adjustment
all_anova = np.array(anova_c+anova_n+anova_p+anova_s)
locs = np.array(list(locate(all_anova, lambda a: a != '')))
pvals = all_anova[locs]
pvals = fdr(pvals)
all_anova[locs] = pvals
anova_adj_c = all_anova[:len(unq_metabolites)]
anova_adj_n = all_anova[len(unq_metabolites):2*len(unq_metabolites)]
anova_adj_p = all_anova[2*len(unq_metabolites):3*len(unq_metabolites)]
anova_adj_s = all_anova[3*len(unq_metabolites):]

all_wilcox = np.array(wilcox_c+wilcox_n+wilcox_p+wilcox_s)
locs = np.array(list(locate(all_wilcox, lambda a: a != '')))
pvals = all_wilcox[locs]
pvals = fdr(pvals)
all_wilcox[locs] = pvals
wilcox_adj_c = all_wilcox[:len(unq_metabolites)]
wilcox_adj_n = all_wilcox[len(unq_metabolites):2*len(unq_metabolites)]
wilcox_adj_p = all_wilcox[2*len(unq_metabolites):3*len(unq_metabolites)]
wilcox_adj_s = all_wilcox[3*len(unq_metabolites):]

all_ttest = np.array(ttest_c+ttest_n+ttest_p+ttest_s)
locs = np.array(list(locate(all_ttest, lambda a: a != '')))
pvals = all_ttest[locs]
pvals = fdr(pvals)
all_ttest[locs] = pvals
ttest_adj_c = all_ttest[:len(unq_metabolites)]
ttest_adj_n = all_ttest[len(unq_metabolites):2*len(unq_metabolites)]
ttest_adj_p = all_ttest[2*len(unq_metabolites):3*len(unq_metabolites)]
ttest_adj_s = all_ttest[3*len(unq_metabolites):]


########################
##### Metabolomics analysis
#improt data
metab_data = pd.read_csv('Calbicans_LRhamnosus_EpithelialCells_Interaction/data/metabolomics_data.csv',index_col=0)

#select the metabolites matched with biolog data
metab_data = metab_data[metab_data.index.isin(unq_metabolites)]

#take related data
lr_ec_6 = ['Blank_6 h','Lr_EC_6 h','Lr_EC_6 h.1','Lr_EC_6 h.2','Lr_EC_6 h.3','Lr_EC_6 h.4']
lr_ec_12 = ['Blank_12 h','Lr_EC_12 h','Lr_EC_12 h.1','Lr_EC_12 h.2','Lr_EC_12 h.3','Lr_EC_12 h.4']
ca_ec_6 = ['Blank_6 h','Ca_EC_6 h','Ca_EC_6 h.1','Ca_EC_6 h.2','Ca_EC_6 h.3','Ca_EC_6 h.4']
ca_ec_12 = ['Blank_12 h','Ca_EC_12 h','Ca_EC_12 h.1','Ca_EC_12 h.2','Ca_EC_12 h.3','Ca_EC_12 h.4']
lr_ec_ca_6 = ['Blank_6 h','Ca_LR_EC_6 h','Ca_LR_EC_6 h.1','Ca_LR_EC_6 h.2','Ca_LR_EC_6 h.3','Ca_LR_EC_6 h.4']
lr_ec_ca_12 = ['Blank_12 h','Ca_LR_EC_12 h','Ca_LR_EC_12 h.1','Ca_LR_EC_12 h.2','Ca_LR_EC_12 h.3','Ca_LR_EC_12 h.4']
lr_6 = ['Blank_6 h','Lr_6 h','Lr_6 h.1','Lr_6 h.2','Lr_6 h.3','Lr_6 h.4']
lr_12 = ['Blank_12 h','Lr_12 h','Lr_12 h.1','Lr_12 h.2','Lr_12 h.3','Lr_12 h.4']
ca_6 = ['Blank_6 h','Ca_6 h','Ca_6 h.1','Ca_6 h.2','Ca_6 h.3','Ca_6 h.4']
ca_12 = ['Blank_12 h','Ca_12 h','Ca_12 h.1','Ca_12 h.2','Ca_12 h.3','Ca_12 h.4']
ec_6 = ['Blank_6 h','EC_6 h','EC_6 h.1','EC_6 h.2','EC_6 h.3','EC_6 h.4']
ec_12 = ['Blank_12 h','EC_12 h','EC_12 h.1','EC_12 h.2','EC_12 h.3','EC_12 h.4']



#comparison to blank
pvalues_all = list()
wilcoxon_all = list()
mean_values_all = list()
for c in [lr_6, lr_12, ca_6, ca_12, ec_6, ec_12, lr_ec_6, lr_ec_12, ca_ec_6, ca_ec_12, lr_ec_ca_6, lr_ec_ca_12]:
    sub_df = metab_data.iloc[:,metab_data.columns.isin(c)]
    sub_df = sub_df.fillna(0)#replace NA values with 0
    sub_df = sub_df + 1 #for logarithmic transformation in case the average is zero
    sub_df_log = pd.DataFrame(columns = c[1:],index = sub_df.index)
    #stats.ttest_ind()
    pvalues = list()
    wilcoxons = list()
    for i in range(sub_df.shape[0]):
        for j in range(1,sub_df.shape[1]):
            sub_df_log.iloc[i,j-1] = np.log2(sub_df.iloc[i,j]/sub_df.iloc[i,0])
                
        pval = stats.ttest_1samp(sub_df.iloc[i,:],popmean=0)[1]
        wilcox = stats.wilcoxon(sub_df.iloc[i,:])[1]
        pvalues.append(pval)
        wilcoxons.append(wilcox)

    mean_values_all.append(list(sub_df_log.mean(axis=1)))
        
    pvalues_all.append(pvalues)
    wilcoxon_all.append(wilcoxons)

del lr_ec_6[0]
del lr_ec_12[0]
del lr_ec_ca_6[0]
del lr_ec_ca_12[0]
del lr_6[0]
del lr_12[0]
del ca_6[0]
del ca_12[0]
del ec_6[0]
del ec_12[0]
del ca_ec_6[0]
del ca_ec_12[0]

#ec_lr vs ec_lr_ca 6 and 12 h
for c in [lr_ec_6+lr_ec_ca_6, lr_ec_12+lr_ec_ca_12, ec_6+lr_ec_6, ec_12+lr_ec_12,
          lr_6+lr_ec_6, lr_12+lr_ec_12, ca_6+ca_ec_6, ca_12+ca_ec_12,ec_6+lr_ec_ca_6,ec_12+lr_ec_ca_12,ca_ec_6+lr_ec_ca_6,ca_ec_12+lr_ec_ca_12]:
    sub_df = metab_data.iloc[:,metab_data.columns.isin(c)]
    sub_df = sub_df.fillna(0)#replace NA values with 0
    sub_df = sub_df + 1 #for logarithmic transformation in case the average is zero
    #sub_df_log = pd.DataFrame(columns = c[1:],index = sub_df.index)
    pvalues = list()
    wilcoxons = list()
    mean_values = list()
    for i in range(sub_df.shape[0]):
        pvalues.append(stats.ttest_ind(sub_df.iloc[i,:5],sub_df.iloc[i,5:])[1])
        try:
            wilcoxons.append(stats.wilcoxon(sub_df.iloc[i,:5],sub_df.iloc[i,5:])[1])
        except ValueError:
            wilcoxons.append(1)
        mean_values.append(np.log2(sub_df.iloc[i,5:].mean()/sub_df.iloc[i,:5].mean()))
        
    
    wilcoxon_all.append(wilcoxons)
    pvalues_all.append(pvalues)
    mean_values_all.append(mean_values)

pvalues_all_fdr = list()
for i in pvalues_all:
    pvalues_all_fdr = pvalues_all_fdr+i

wilcoxon_all_fdr = list()
for i in wilcoxon_all:
    wilcoxon_all_fdr = wilcoxon_all_fdr+i
    
pvalues_all_corrected = fdr(pvalues_all_fdr)

wilcoxon_all_corrected = fdr(wilcoxon_all_fdr)





#results.to_csv('Calbicans_LRhamnosus_EpithelialCells_Interaction/added_files_by_me/results/biolog_vs_metabolome.csv')
#FVA for three models

metab_data_big = pd.read_csv('Calbicans_LRhamnosus_EpithelialCells_Interaction/data/metabolomics_data.csv',index_col=0)

recon = cobra.io.load_matlab_model('Calbicans_LRhamnosus_EpithelialCells_Interaction/models/Recon3DModel_301.mat')
candida = cobra.io.read_sbml_model('Calbicans_LRhamnosus_EpithelialCells_Interaction/models/Suppl1_candida_albicans_GSMM.xml')
lacto = cobra.io.load_matlab_model('Calbicans_LRhamnosus_EpithelialCells_Interaction/models/Lactobacillus_rhamnosus_LMS2_1.mat')

candida_met = pd.read_csv('Calbicans_LRhamnosus_EpithelialCells_Interaction/data/models/Candida_albicans_CORECO_to_AGORA_map.csv',index_col=0)
agora_mets = pd.read_table('Calbicans_LRhamnosus_EpithelialCells_Interaction/added_files_by_me/data/recon-store-metabolites-1.tsv')
manual_match = pd.read_table('Calbicans_LRhamnosus_EpithelialCells_Interaction/data/models/manual_met_match.csv')
manual_match = list(manual_match['agora_met; name'])

ids = agora_mets[agora_mets['pubChemId'].isin(list(metab_data_big['PUBCHEM']))]

#lacto essentials for growth
lacto_essentials = ['EX_ca2(e)',
                    'EX_cl(e)',
                    'EX_cobalt2(e)',
                    'EX_cu2(e)',
                    'EX_fe2(e)',
                    'EX_fe3(e)',
                    'EX_k(e)',
                    'EX_mg2(e)',
                    'EX_mn2(e)',
                    'EX_zn2(e)',
                    'EX_h2o(e)',
                    'EX_na1(e)']

recon_essentials = ['EX_ca2[e]',
                    'EX_cl[e]',
                    'EX_cobalt2[e]',
                    'EX_cu2[e]',
                    'EX_fe2[e]',
                    'EX_fe3[e]',
                    'EX_k[e]',
                    'EX_mg2[e]',
                    'EX_mn2[e]',
                    'EX_zn2[e]',
                    'EX_h2o[e]',
                    'EX_na1[e]']

candida_essentials = ['Ex_CHEBI29033']

# met mapping
agora_mets = pd.read_table('Calbicans_LRhamnosus_EpithelialCells_Interaction/added_files_by_me/data/recon-store-metabolites-1.tsv')
manual_match = pd.read_table('Calbicans_LRhamnosus_EpithelialCells_Interaction/data/models/manual_met_match.csv')
manual_match = list(manual_match['agora_met; name'])

ids = agora_mets[agora_mets['keggId'].isin(list(metab_data_big['KEGG']))]

agora_ids = list()
id_by_kegg = list()
id_by_pub = list()
id_by_hmdb = list()

for i in range(metab_data_big.shape[0]):
    if not pd.isna(metab_data_big['KEGG'][i]) and agora_mets[agora_mets['keggId']==metab_data_big['KEGG'][i]].shape[0]!=0:
        id_by_kegg.append(agora_mets[agora_mets['keggId']==metab_data_big['KEGG'][i]]['abbreviation'].iloc[0])
    else:
        id_by_kegg.append('')

    if not pd.isna(metab_data_big['PUBCHEM'][i]) and agora_mets[agora_mets['pubChemId']==metab_data_big['PUBCHEM'][i]].shape[0]!=0:
        id_by_pub.append(agora_mets[agora_mets['pubChemId']==metab_data_big['PUBCHEM'][i]]['abbreviation'].iloc[0])
    else:
        id_by_pub.append('')

    #hmdb modification
    if pd.isna(metab_data_big['HMDB'][i]):
        id_by_hmdb.append('')
        continue
    else:    
        hmdb_id = 'HMDB00' + metab_data_big['HMDB'][i][4:]
    
    if not pd.isna(hmdb_id) and agora_mets[agora_mets['hmdb']==hmdb_id].shape[0]!=0:
        id_by_hmdb.append(agora_mets[agora_mets['hmdb']==hmdb_id]['abbreviation'].iloc[0])
    else:
        id_by_hmdb.append('')

#dataframe of mapped metabolites to agora ids by using three databases.
#In case columns are not identical in each row, it should interrogated
#to see which agora ID is correct for the corresponding metabolite
df_mapped_ids = pd.DataFrame({'id_by_hmdb':id_by_hmdb,'id_by_kegg':id_by_kegg,'id_by_pub':id_by_pub})

#consider the identical ones to be agora id
#print rest to investigate manually
for i in range(df_mapped_ids.shape[0]):
    if all(df_mapped_ids.iloc[i,:]==df_mapped_ids.iloc[i,0]) and all(df_mapped_ids.iloc[i,:]!=''):
        agora_ids.append(df_mapped_ids.iloc[i,0])
    else:
        agora_ids.append('')
        print(df_mapped_ids.iloc[i,:])
        print(metab_data_big.index[i])
        print()

#manual consideration for the inconsistent rows
agora_ids[6] = '1mncam'
agora_ids[17] = '23dhmb'
agora_ids[21] = 'M00653'
agora_ids[27] = 'bhb'
agora_ids[29] = '3hmp'
agora_ids[33] = '3mop'
agora_ids[34] = 'Lcyst'
agora_ids[37] = '4gubut'
agora_ids[52] = 'alltn'
agora_ids[57] = 'lipoate'
agora_ids[58] = 'xylt'
agora_ids[61] = 'arg_L'
agora_ids[69] = 'crn'
agora_ids[80] = '3sala'
agora_ids[89] = 'dhor_S'
agora_ids[93] = 'ethmalac'
agora_ids[95] = 'fru'
agora_ids[99] = 'gluala'
agora_ids[111] = 'glu_L'
agora_ids[115] = 'glyc_R'
agora_ids[117] = 'glyc3p'
agora_ids[118] = 'M02913'
agora_ids[125] = 'guln'
agora_ids[141] = 'lac_L'
agora_ids[144] = 'mal_L'
agora_ids[153] = 'mev_R'
agora_ids[154] = 'mvlac'
agora_ids[159] = 'CE1556'
agora_ids[161] = 'acgam'
agora_ids[162] = 'acglu'
agora_ids[165] = 'acile_L'
agora_ids[167] = 'C02712'
agora_ids[168] = 'acnam'
agora_ids[186] = 'orn'
agora_ids[188] = 'pcs'
agora_ids[191] = 'pnto_R'
agora_ids[197] = 'plac'
agora_ids[202] = 'Lpipecol'
agora_ids[208] = 'pydxn'
agora_ids[211] = 'rbt'
agora_ids[214] = '1pyr5c'
agora_ids[215] = 'M03165'
agora_ids[218] = 'so4'


#mapping using manual match file
for i in manual_match:
    s = i.split(';')
    if s[0] not in agora_ids and s[1][1:] in list(metab_data_big.index):
        agora_ids[list(metab_data_big.index).index(s[1][1:])] = s[0]
agora_ids[list(metab_data_big.index).index("mannonate*")] = "mana"
agora_ids[list(metab_data_big.index).index("aconitate [cis or trans]")] = "acon_C"
agora_ids[list(metab_data_big.index).index("arabonate/xylonate")] = "xylnt"
agora_ids[list(metab_data_big.index).index("diacetylspermidine*")] = "CE1059"
agora_ids[list(metab_data_big.index).index("glycerophosphoglycerol")] = "g3pg"
agora_ids[list(metab_data_big.index).index("palmitoyl sphingomyelin (d18:1/16:0)")] = "sphmyln18116_hs"
agora_ids[list(metab_data_big.index).index("aconitate [cis or trans]")] = "acon_C"
agora_ids[20] = '2hb'
agora_ids[55] = 'HC00591'



## met mapping for candida
candida_ids = list()
for i in agora_ids:
    if i != '' and i in list(candida_met['abbreviation']):
        candida_ids.append(candida_met[candida_met['abbreviation'] == i].iloc[0,2])
        if pd.isna(candida_ids[len(candida_ids)-1]):
            candida_ids[len(candida_ids)-1] = ''
    else:
        candida_ids.append('')
# check if some ids can be found by using kegg ids
kegg_ids = list(metab_data_big['KEGG'])
for i in range(len(candida_ids)):
    if candida_ids[i]=='' and kegg_ids[i] in candida.metabolites:
        print('yes')
        if 'Ex_' + kegg_ids[i] in candida.reactions:
            candida_ids[i] = kegg_ids[i]

candida_ex_ids = ['Ex_'+i if i!='' else '' for i in candida_ids]

####
#make list of exchange reaction ids based on agora_ids separately for each 
recon_ex_ids = list()
for i in agora_ids:
    if i!='':
        recon_ex_ids.append('EX_'+i+'[e]')
    else:
        recon_ex_ids.append('')

lacto_ex_ids = list()
for i in agora_ids:
    if i!='':
        lacto_ex_ids.append('EX_'+i+'(e)')
    else:
        lacto_ex_ids.append('')



#convert metab_data to influx values. each value divided by the maximum value in the matab_data multiplied to 1000
metab_data_big.iloc[:,3:] = metab_data_big.iloc[:,3:].fillna(0)
metab_data_big.iloc[:,3:] = metab_data_big.iloc[:,3:]/metab_data_big.iloc[:,3:].max().max()*-1000


#set up the diet with respect to metabolome data and make list of ex reactions for each model with respect to metabolome data (available compounds in models)
#take the related media and uptake rates needed for each model
blank_6 = metab_data_big.loc[:,metab_data_big.columns.isin(['Blank_6 h'])].values
blank_12 = metab_data_big.loc[:,metab_data_big.columns.isin(['Blank_12 h'])].values
ec_spent_6 = metab_data_big.loc[:,metab_data_big.columns.isin(['EC_6 h', 'EC_6 h.1', 'EC_6 h.2', 'EC_6 h.3', 'EC_6 h.4'])].values
ec_spent_12 = metab_data_big.loc[:,metab_data_big.columns.isin(['EC_12 h', 'EC_12 h.1', 'EC_12 h.2', 'EC_12 h.3', 'EC_12 h.4'])].values

def set_diet_ref_based(model,ex,upt):
    ex_rxns = list()
    for i in model.reactions:
        if i.id[:3]=='EX_' or i.id[:3]=='Ex_':
            #if i.id in ex and present[ex.index(i.id)]:
            if i.id in ex:
                i.bounds=(upt[ex.index(i.id)].mean(),1000)
                ex_rxns.append(i)
            else:
                i.bounds = (0,1000)

            #if the model is lacto, allow influx of the essentials
            if model.id=='model' and i.id in lacto_essentials:
                i.bounds=(-1000,1000)
    return model,ex_rxns

#blank6
recon,recon_ex_rxns_blank_6 = set_diet_ref_based(recon,recon_ex_ids,blank_6)
lacto,lacto_ex_rxns_blank_6 = set_diet_ref_based(lacto,lacto_ex_ids,blank_6)

r = recon.reactions.get_by_id('biomass_maintenance')
sol = recon.optimize()
r.bounds = (sol.objective_value,sol.objective_value)
for i in recon.reactions:
    if i.id[:2]=='EX':
        i.lower_bound = -1000
r = lacto.reactions.get_by_id('biomass525')
sol = lacto.optimize()
r.bounds = (sol.objective_value,sol.objective_value)
for i in lacto.reactions:
    if i.id[:2]=='EX':
        i.lower_bound = -1000
r.bounds = (0,0)

fva_recon_blank_6 = pd.DataFrame()
for i in range(len(recon_ex_rxns_blank_6)):
    fva_recon_blank_6 = pd.concat([flux_variability_analysis(recon,recon_ex_rxns_blank_6[i],fraction_of_optimum = 0.9),
                                   fva_recon_blank_6])
fva_lacto_blank_6 = pd.DataFrame()
for i in range(len(lacto_ex_rxns_blank_6)):
    fva_lacto_blank_6 = pd.concat([flux_variability_analysis(lacto,lacto_ex_rxns_blank_6[i],fraction_of_optimum = 0.9),
                                   fva_lacto_blank_6])

r = recon.reactions.get_by_id('biomass_maintenance')
r.bounds = (0,1000)
r = lacto.reactions.get_by_id('biomass525')
r.bounds = (0,1000)

#blank12
recon,recon_ex_rxns_blank_12 = set_diet_ref_based(recon,recon_ex_ids,blank_12)
lacto,lacto_ex_rxns_blank_12 = set_diet_ref_based(lacto,lacto_ex_ids,blank_12)

r = recon.reactions.get_by_id('biomass_maintenance')
sol = recon.optimize()
r.bounds = (sol.objective_value,sol.objective_value)
for i in recon.reactions:
    if i.id[:2]=='EX':
        i.lower_bound = -1000
r = lacto.reactions.get_by_id('biomass525')
sol = lacto.optimize()
r.bounds = (sol.objective_value,sol.objective_value)
for i in lacto.reactions:
    if i.id[:2]=='EX':
        i.lower_bound = -1000

fva_recon_blank_12 = pd.DataFrame()
for i in range(len(recon_ex_rxns_blank_12)):
    fva_recon_blank_12 = pd.concat([flux_variability_analysis(recon,recon_ex_rxns_blank_12[i],fraction_of_optimum = 0.9),
                                    fva_recon_blank_12])
fva_lacto_blank_12 = pd.DataFrame()
for i in range(len(lacto_ex_rxns_blank_12)):
    fva_lacto_blank_12 = pd.concat([flux_variability_analysis(lacto,lacto_ex_rxns_blank_12[i],fraction_of_optimum = 0.9),
                                    fva_lacto_blank_12])

r = recon.reactions.get_by_id('biomass_maintenance')
r.bounds = (0,1000)
r = lacto.reactions.get_by_id('biomass525')
r.bounds = (0,1000)

#ec_spent_6
recon,recon_ex_rxns_ec_6 = set_diet_ref_based(recon,recon_ex_ids,ec_spent_6)
lacto,lacto_ex_rxns_ec_6 = set_diet_ref_based(lacto,lacto_ex_ids,ec_spent_6)

r = recon.reactions.get_by_id('biomass_maintenance')
sol = recon.optimize()
r.bounds = (sol.objective_value,sol.objective_value)
for i in recon.reactions:
    if i.id[:2]=='EX':
        i.lower_bound = -1000
r = lacto.reactions.get_by_id('biomass525')
sol = lacto.optimize()
r.bounds = (sol.objective_value,sol.objective_value)
for i in lacto.reactions:
    if i.id[:2]=='EX':
        i.lower_bound = -1000

fva_recon_ec_6 = pd.DataFrame()
for i in range(len(recon_ex_rxns_ec_6)):
    fva_recon_ec_6 = pd.concat([flux_variability_analysis(recon,recon_ex_rxns_ec_6[i],fraction_of_optimum = 0.9),
                                fva_recon_ec_6])
fva_lacto_ec_6 = pd.DataFrame()
for i in range(len(lacto_ex_rxns_ec_6)):
    fva_lacto_ec_6 = pd.concat([flux_variability_analysis(lacto,lacto_ex_rxns_ec_6[i],fraction_of_optimum = 0.9),
                                   fva_lacto_ec_6])

r = recon.reactions.get_by_id('biomass_maintenance')
r.bounds = (0,1000)
r = lacto.reactions.get_by_id('biomass525')
r.bounds = (0,1000)

#ec_spent_12
recon,recon_ex_rxns_ec_12 = set_diet_ref_based(recon,recon_ex_ids,ec_spent_12)
lacto,lacto_ex_rxns_ec_12 = set_diet_ref_based(lacto,lacto_ex_ids,ec_spent_12)

r = recon.reactions.get_by_id('biomass_maintenance')
sol = recon.optimize()
r.bounds = (sol.objective_value,sol.objective_value)
for i in recon.reactions:
    if i.id[:2]=='EX':
        i.lower_bound = -1000
r = lacto.reactions.get_by_id('biomass525')
sol = lacto.optimize()
r.bounds = (sol.objective_value,sol.objective_value)
for i in lacto.reactions:
    if i.id[:2]=='EX':
        i.lower_bound = -1000

fva_recon_ec_12 = pd.DataFrame()
for i in range(len(recon_ex_rxns_ec_12)):
    fva_recon_ec_12 = pd.concat([flux_variability_analysis(recon,recon_ex_rxns_ec_12[i],fraction_of_optimum = 0.9),
                                 fva_recon_ec_12])
fva_lacto_ec_12 = pd.DataFrame()
for i in range(len(lacto_ex_rxns_ec_12)):
    fva_lacto_ec_12 = pd.concat([flux_variability_analysis(lacto,lacto_ex_rxns_ec_12[i],fraction_of_optimum = 0.9),
                                    fva_lacto_ec_12])

r = recon.reactions.get_by_id('biomass_maintenance')
r.bounds = (0,1000)
r = lacto.reactions.get_by_id('biomass525')
r.bounds = (0,1000)

fva_df = pd.DataFrame(index=unq_metabolites,data={'EC_blank6':['none']*len(unq_metabolites),
                                                  'Lr_blank6':['none']*len(unq_metabolites),
                                                  'EC_blank12':['none']*len(unq_metabolites),
                                                  'Lr_blank12':['none']*len(unq_metabolites),
                                                  'EC_ec6':['none']*len(unq_metabolites),
                                                  'Lr_ec6':['none']*len(unq_metabolites),
                                                  'EC_ec12':['none']*len(unq_metabolites),
                                                  'EC_ec12':['none']*len(unq_metabolites)})
                                                  

fva_recon_ec_6[abs(fva_recon_ec_6)<0.001] = 0
fva_recon_ec_12[abs(fva_recon_ec_12)<0.001] = 0
fva_recon_blank_6[abs(fva_recon_blank_6)<0.001] = 0
fva_recon_blank_12[abs(fva_recon_blank_12)<0.001] = 0
fva_lacto_ec_6[abs(fva_lacto_ec_6)<0.001] = 0
fva_lacto_ec_12[abs(fva_lacto_ec_12)<0.001] = 0
fva_lacto_blank_6[abs(fva_lacto_blank_6)<0.001] = 0
fva_lacto_blank_12[abs(fva_lacto_blank_12)<0.001] = 0

for i in unq_metabolites:
    
    ind_big = list(metab_data_big.index).index(i)
    
    if lacto_ex_ids[ind_big] != '' and lacto_ex_ids[ind_big] in fva_lacto_blank_6.index:
        if fva_lacto_blank_6.loc[lacto_ex_ids[ind_big],'maximum']> 0 and fva_lacto_blank_6.loc[lacto_ex_ids[ind_big],'minimum']<0:
            fva_df.loc[i,'Lr_blank6'] = 'both'
        elif fva_lacto_blank_6.loc[lacto_ex_ids[ind_big],'maximum']>0 and fva_lacto_blank_6.loc[lacto_ex_ids[ind_big],'minimum']>=0:
            fva_df.loc[i,'Lr_blank6'] = 'producible'
        elif abs(fva_lacto_blank_6.loc[lacto_ex_ids[ind_big],'maximum'])<=0 and fva_lacto_blank_6.loc[lacto_ex_ids[ind_big],'minimum']<0:
            fva_df.loc[i,'Lr_blank6'] = 'consumable'

    if recon_ex_ids[ind_big] != '' and recon_ex_ids[ind_big] in fva_recon_blank_6.index:
        if fva_recon_blank_6.loc[recon_ex_ids[ind_big],'maximum']>0 and fva_recon_blank_6.loc[recon_ex_ids[ind_big],'minimum']<0:
            fva_df.loc[i,'EC_blank6'] = 'both'
        elif fva_recon_blank_6.loc[recon_ex_ids[ind_big],'maximum']>0 and fva_recon_blank_6.loc[recon_ex_ids[ind_big],'minimum']>=0:
            fva_df.loc[i,'EC_blank6'] = 'producible'
        elif fva_recon_blank_6.loc[recon_ex_ids[ind_big],'maximum']<=0 and fva_recon_blank_6.loc[recon_ex_ids[ind_big],'minimum']<0:
            fva_df.loc[i,'EC_blank6'] = 'consumable'

    if lacto_ex_ids[ind_big] != '' and lacto_ex_ids[ind_big] in fva_lacto_ec_6.index:
        if fva_lacto_ec_6.loc[lacto_ex_ids[ind_big],'maximum']>0 and fva_lacto_ec_6.loc[lacto_ex_ids[ind_big],'minimum']<0:
            fva_df.loc[i,'Lr_ec6'] = 'both'
        elif fva_lacto_ec_6.loc[lacto_ex_ids[ind_big],'maximum']>0 and fva_lacto_ec_6.loc[lacto_ex_ids[ind_big],'minimum']>=0:
            fva_df.loc[i,'Lr_ec6'] = 'producible'
        elif fva_lacto_ec_6.loc[lacto_ex_ids[ind_big],'maximum']<=0 and fva_lacto_ec_6.loc[lacto_ex_ids[ind_big],'minimum']<0:
            fva_df.loc[i,'Lr_ec6'] = 'consumable'

    if recon_ex_ids[ind_big] != '' and recon_ex_ids[ind_big] in fva_recon_ec_6.index:
        if fva_recon_ec_6.loc[recon_ex_ids[ind_big],'maximum']>0 and fva_recon_ec_6.loc[recon_ex_ids[ind_big],'minimum']<0:
            fva_df.loc[i,'EC_ec6'] = 'both'
        elif fva_recon_ec_6.loc[recon_ex_ids[ind_big],'maximum']>0 and fva_recon_ec_6.loc[recon_ex_ids[ind_big],'minimum']>=0:
            fva_df.loc[i,'EC_ec6'] = 'producible'
        elif fva_recon_ec_6.loc[recon_ex_ids[ind_big],'maximum']<=0 and fva_recon_ec_6.loc[recon_ex_ids[ind_big],'minimum']<0:
            fva_df.loc[i,'EC_ec6'] = 'consumable'

    if lacto_ex_ids[ind_big] != '' and lacto_ex_ids[ind_big] in fva_lacto_blank_12.index:
        if fva_lacto_blank_12.loc[lacto_ex_ids[ind_big],'maximum']>0 and fva_lacto_blank_12.loc[lacto_ex_ids[ind_big],'minimum']<0:
            fva_df.loc[i,'Lr_blank12'] = 'both'
        elif fva_lacto_blank_12.loc[lacto_ex_ids[ind_big],'maximum']>0 and fva_lacto_blank_12.loc[lacto_ex_ids[ind_big],'minimum']>=0:
            fva_df.loc[i,'Lr_blank12'] = 'producible'
        elif fva_lacto_blank_12.loc[lacto_ex_ids[ind_big],'maximum']<=0 and fva_lacto_blank_12.loc[lacto_ex_ids[ind_big],'minimum']<0:
            fva_df.loc[i,'Lr_blank12'] = 'consumable'

    if recon_ex_ids[ind_big] != '' and recon_ex_ids[ind_big] in fva_recon_blank_12.index:
        if fva_recon_blank_12.loc[recon_ex_ids[ind_big],'maximum']>0 and fva_recon_blank_12.loc[recon_ex_ids[ind_big],'minimum']<0:
            fva_df.loc[i,'EC_blank12'] = 'both'
        elif fva_recon_blank_12.loc[recon_ex_ids[ind_big],'maximum']>0 and fva_recon_blank_12.loc[recon_ex_ids[ind_big],'minimum']>=0:
            fva_df.loc[i,'EC_blank12'] = 'producible'
        elif fva_recon_blank_12.loc[recon_ex_ids[ind_big],'maximum']<=0 and fva_recon_blank_12.loc[recon_ex_ids[ind_big],'minimum']<0:
            fva_df.loc[i,'EC_blank12'] = 'consumable'

    if lacto_ex_ids[ind_big] != '' and lacto_ex_ids[ind_big] in fva_lacto_ec_12.index:
        if fva_lacto_ec_12.loc[lacto_ex_ids[ind_big],'maximum']>0 and fva_lacto_ec_12.loc[lacto_ex_ids[ind_big],'minimum']<0:
            fva_df.loc[i,'Lr_ec12'] = 'both'
        elif fva_lacto_ec_12.loc[lacto_ex_ids[ind_big],'maximum']>0 and fva_lacto_ec_12.loc[lacto_ex_ids[ind_big],'minimum']>=0:
            fva_df.loc[i,'Lr_ec12'] = 'producible'
        elif fva_lacto_ec_12.loc[lacto_ex_ids[ind_big],'maximum']<=0 and fva_lacto_ec_12.loc[lacto_ex_ids[ind_big],'minimum']<0:
            fva_df.loc[i,'Lr_ec12'] = 'consumable'

    if recon_ex_ids[ind_big] != '' and recon_ex_ids[ind_big] in fva_recon_ec_12.index:
        if fva_recon_ec_12.loc[recon_ex_ids[ind_big],'maximum']>0 and fva_recon_ec_12.loc[recon_ex_ids[ind_big],'minimum']<0:
            fva_df.loc[i,'EC_ec12'] = 'both'
        elif fva_recon_ec_12.loc[recon_ex_ids[ind_big],'maximum']>0 and fva_recon_ec_12.loc[recon_ex_ids[ind_big],'minimum']>=0:
            fva_df.loc[i,'EC_ec12'] = 'producible'
        elif fva_recon_ec_12.loc[recon_ex_ids[ind_big],'maximum']<=0 and fva_recon_ec_12.loc[recon_ex_ids[ind_big],'minimum']<0:
            fva_df.loc[i,'EC_ec12'] = 'consumable'
    

fva_df[fva_df=='producible'] = 'secretion'
fva_df[fva_df=='consumable'] = 'uptake'


    
#########################################################
#saving the results
results1 = pd.DataFrame({'anova_adj_c':anova_adj_c,
                        'anova_adj_n':anova_adj_n,
                        'anova_adj_p':anova_adj_p,
                        'anova_adj_s':anova_adj_s,
                        'anova_c':anova_c,
                        'anova_n':anova_n,
                        'anova_p':anova_p,
                        'anova_s':anova_s,
                        'wilcox_c':wilcox_c,
                        'wilcox_n':wilcox_n,
                        'wilcox_p':wilcox_p,
                        'wilcox_s':wilcox_s,
                        'wilcox_adj_c':wilcox_adj_c,
                        'wilcox_adj_n':wilcox_adj_n,
                        'wilcox_adj_p':wilcox_adj_p,
                        'wilcox_adj_s':wilcox_adj_s,
                        'ttest_c':ttest_c,
                        'ttest_n':ttest_n,
                        'ttest_p':ttest_p,
                        'ttest_s':ttest_s,
                        'ttest_adj_c':ttest_adj_c,
                        'ttest_adj_n':ttest_adj_n,
                        'ttest_adj_p':ttest_adj_p,
                        'ttest_adj_s':ttest_adj_s,
                        'rel_growths_c':rel_growths_c,
                        'rel_growths_n':rel_growths_n,
                        'rel_growths_p':rel_growths_p,
                        'rel_growths_s':rel_growths_s},
                        index = unq_metabolites)
                   
results2 = pd.DataFrame({'lr_6_blank':mean_values_all[0],
                        'lr_12_blank':mean_values_all[1],
                        'ca_6_blank':mean_values_all[2],
                        'ca_12_blank':mean_values_all[3],
                        'ec_6_blank':mean_values_all[4],
                        'ec_12_blank':mean_values_all[5],
                        'lr_ec_6_blank':mean_values_all[6],
                        'lr_ec_12_blank':mean_values_all[7],
                        'ca_ec_6_blank':mean_values_all[8],
                        'ca_ec_12_blank':mean_values_all[9],
                        'lr_ec_ca_6_blank':mean_values_all[10],
                        'lr_ec_ca_12_blank':mean_values_all[11],
                        'lr_ec_6_vs_lr_ec_ca_6':mean_values_all[12],
                        'lr_ec_12_vs_lr_ec_ca_12':mean_values_all[13],
                        'ec_6_vs_lr_ec_6':mean_values_all[14],
                        'ec_12_vs_lr_ec_12':mean_values_all[15],
                        'lr_6_vs_lr_ec_6':mean_values_all[16],
                        'lr_12_vs_lr_ec_12':mean_values_all[17],
                        'ca_6_vs_ca_ec_6':mean_values_all[18],
                        'ca_12_vs_ca_ec_12':mean_values_all[19],
                        'ec_6+lr_ec_ca_6':mean_values_all[20],
                        'ec_12+lr_ec_ca_12':mean_values_all[21],
                        'ca_ec_6+lr_ec_ca_6':mean_values_all[22],
                        'ca_ec_12+lr_ec_ca_12':mean_values_all[23],
                        
                        'lr_6_blank_wilcox':wilcoxon_all[0],
                        'lr_12_blank_wilcox':wilcoxon_all[1],
                        'ca_6_blank_wilcox':wilcoxon_all[2],
                        'ca_12_blank_wilcox':wilcoxon_all[3],
                        'ec_6_blank_wilcox':wilcoxon_all[4],
                        'ec_12_blank_wilcox':wilcoxon_all[5],
                        'lr_ec_6_blank_wilcox':wilcoxon_all[6],
                        'lr_ec_12_blank_wilcox':wilcoxon_all[7],
                        'ca_ec_6_blank_wilcox':wilcoxon_all[8],
                        'ca_ec_12_blank_wilcox':wilcoxon_all[9],
                        'lr_ec_ca_6_blank_wilcox':wilcoxon_all[10],
                        'lr_ec_ca_12_blank_wilcox':wilcoxon_all[11],
                        'lr_ec_6_vs_lr_ec_ca_6_wilcox':wilcoxon_all[12],
                        'lr_ec_12_vs_lr_ec_ca_12_wilcox':wilcoxon_all[13],
                        'ec_6_vs_lr_ec_6_wilcox':wilcoxon_all[14],
                        'ec_12_vs_lr_ec_12_wilcox':wilcoxon_all[15],
                        'lr_6_vs_lr_ec_6_wilcox':wilcoxon_all[16],
                        'lr_12_vs_lr_ec_12_wilcox':wilcoxon_all[17],
                        'ca_6_vs_ca_ec_6_wilcox':wilcoxon_all[18],
                        'ca_12_vs_ca_ec_12_wilcox':wilcoxon_all[19],
                        'ec_6+lr_ec_ca_6_wilcox':wilcoxon_all[20],
                        'ec_12+lr_ec_ca_12_wilcox':wilcoxon_all[21],
                        'ca_ec_6+lr_ec_ca_6_wilcox':wilcoxon_all[22],
                        'ca_ec_12+lr_ec_ca_12_wilcox':wilcoxon_all[23],
                        
                        'lr_6_blank_ttest_adj':pvalues_all_corrected[:79],
                        'lr_12_blank_ttest_adj':pvalues_all_corrected[79:2*79],
                        'ca_6_blank_ttest_adj':pvalues_all_corrected[2*79:3*79],
                        'ca_12_blank_ttest_adj':pvalues_all_corrected[3*79:4*79],
                        'ec_6_blank_ttest_adj':pvalues_all_corrected[4*79:5*79],
                        'ec_12_blank_ttest_adj':pvalues_all_corrected[5*79:6*79],
                        'lr_ec_6_blank_ttest_adj':pvalues_all_corrected[6*79:7*79],
                        'lr_ec_12_blank_ttest_adj':pvalues_all_corrected[7*79:8*79],
                        'ca_ec_6_blank_ttest_adj':pvalues_all_corrected[8*79:9*79],
                        'ca_ec_12_blank_ttest_adj':pvalues_all_corrected[9*79:10*79],
                        'lr_ec_ca_6_blank_ttest_adj':pvalues_all_corrected[10*79:11*79],
                        'lr_ec_ca_12_blank_ttest_adj':pvalues_all_corrected[11*79:12*79],
                        'lr_ec_6_vs_lr_ec_ca_6_ttest_adj':pvalues_all_corrected[12*79:13*79],
                        'lr_ec_12_vs_lr_ec_ca_12_ttest_adj':pvalues_all_corrected[13*79:14*79],
                        'ec_6_vs_lr_ec_6_ttest_adj':pvalues_all_corrected[14*79:15*79],
                        'ec_12_vs_lr_ec_12_ttest_adj':pvalues_all_corrected[15*79:16*79],
                        'lr_6_vs_lr_ec_6_ttest_adj':pvalues_all_corrected[16*79:17*79],
                        'lr_12_vs_lr_ec_12_ttest_adj':pvalues_all_corrected[17*79:18*79],
                        'ca_6_vs_ca_ec_6_ttest_adj':pvalues_all_corrected[18*79:19*79],
                        'ca_12_vs_ca_ec_12_ttest_adj':pvalues_all_corrected[19*79:20*79],
                        'ec_6+lr_ec_ca_6_ttest_adj':pvalues_all_corrected[20*79:21*79],
                        'ec_12+lr_ec_ca_12_ttest_adj':pvalues_all_corrected[21*79:22*79],
                        'ca_ec_6+lr_ec_ca_6_ttest_adj':pvalues_all_corrected[22*79:23*79],
                        'ca_ec_12+lr_ec_ca_12_ttest_adj':pvalues_all_corrected[23*79:],

                        'lr_6_blank_ttest':pvalues_all[0],
                        'lr_12_blank_ttest':pvalues_all[1],
                        'ca_6_blank_ttest':pvalues_all[2],
                        'ca_12_blank_ttest':pvalues_all[3],
                        'ec_6_blank_ttest':pvalues_all[4],
                        'ec_12_blank_ttest':pvalues_all[5],
                        'lr_ec_6_blank_ttest':pvalues_all[6],
                        'lr_ec_12_blank_ttest':pvalues_all[7],
                        'ca_ec_6_blank_ttest':pvalues_all[8],
                        'ca_ec_12_blank_ttest':pvalues_all[9],
                        'lr_ec_ca_6_blank_ttest':pvalues_all[10],
                        'lr_ec_ca_12_blank_ttest':pvalues_all[11],
                        'lr_ec_6_vs_lr_ec_ca_6_ttest':pvalues_all[12],
                        'lr_ec_12_vs_lr_ec_ca_12_ttest':pvalues_all[13],
                        'ec_6_vs_lr_ec_6_ttest':pvalues_all[14],
                        'ec_12_vs_lr_ec_12_ttest':pvalues_all[15],
                        'lr_6_vs_lr_ec_6_ttest':pvalues_all[16],
                        'lr_12_vs_lr_ec_12_ttest':pvalues_all[17],
                        'ca_6_vs_ca_ec_6_ttest':pvalues_all[18],
                        'ca_12_vs_ca_ec_12_ttest':pvalues_all[19],
                        'ec_6+lr_ec_ca_6_ttest':pvalues_all[20],
                        'ec_12+lr_ec_ca_12_ttest':pvalues_all[21],
                        'ca_ec_6+lr_ec_ca_6_ttest':pvalues_all[22],
                        'ca_ec_12+lr_ec_ca_12_ttest':pvalues_all[23], 

                        'lr_6_blank_wilcox_adj':wilcoxon_all_corrected[:79],
                        'lr_12_blank_wilcox_adj':wilcoxon_all_corrected[79:2*79],
                        'ca_6_blank_wilcox_adj':wilcoxon_all_corrected[2*79:3*79],
                        'ca_12_blank_wilcox_adj':wilcoxon_all_corrected[3*79:4*79],
                        'ec_6_blank_wilcox_adj':wilcoxon_all_corrected[4*79:5*79],
                        'ec_12_blank_wilcox_adj':wilcoxon_all_corrected[5*79:6*79],
                        'lr_ec_6_blank_wilcox_adj':wilcoxon_all_corrected[6*79:7*79],
                        'lr_ec_12_blank_wilcox_adj':wilcoxon_all_corrected[7*79:8*79],
                        'ca_ec_6_blank_wilcox_adj':wilcoxon_all_corrected[8*79:9*79],
                        'ca_ec_12_blank_wilcox_adj':wilcoxon_all_corrected[9*79:10*79],
                        'lr_ec_ca_6_blank_wilcox_adj':wilcoxon_all_corrected[10*79:11*79],
                        'lr_ec_ca_12_blank_wilcox_adj':wilcoxon_all_corrected[11*79:12*79],
                        'lr_ec_6_vs_lr_ec_ca_6_wilcox_adj':wilcoxon_all_corrected[12*79:13*79],
                        'lr_ec_12_vs_lr_ec_ca_12_wilcox_adj':wilcoxon_all_corrected[13*79:14*79],
                        'ec_6_vs_lr_ec_6_wilcox_adj':wilcoxon_all_corrected[14*79:15*79],
                        'ec_12_vs_lr_ec_12_wilcox_adj':wilcoxon_all_corrected[15*79:16*79],
                        'lr_6_vs_lr_ec_6_wilcox_adj':wilcoxon_all_corrected[16*79:17*79],
                        'lr_12_vs_lr_ec_12_wilcox_adj':wilcoxon_all_corrected[17*79:18*79],
                        'ca_6_vs_ca_ec_6_wilcox_adj':wilcoxon_all_corrected[18*79:19*79],
                        'ca_12_vs_ca_ec_12_wilcox_adj':wilcoxon_all_corrected[19*79:20*79],
                        'ec_6+lr_ec_ca_6_wilcox_adj':wilcoxon_all_corrected[20*79:21*79],
                        'ec_12+lr_ec_ca_12_wilcox_adj':wilcoxon_all_corrected[21*79:22*79],
                        'ca_ec_6+lr_ec_ca_6_wilcox_adj':wilcoxon_all_corrected[22*79:23*79],
                        'ca_ec_12+lr_ec_ca_12_wilcox_adj':wilcoxon_all_corrected[23*79:]},
                        index = metab_data.index)

results = pd.merge(results1,results2,how='outer',left_index=True,right_index=True)

results.to_csv('Calbicans_LRhamnosus_EpithelialCells_Interaction/results/biolog_vs_metabolome.csv')

fva_df.to_csv('Calbicans_LRhamnosus_EpithelialCells_Interaction/results/fva_lr_iec.csv')
