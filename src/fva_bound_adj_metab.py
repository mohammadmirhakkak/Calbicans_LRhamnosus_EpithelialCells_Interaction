import cobra
import pandas as pd
from cobra import Reaction
from cobra.flux_analysis import flux_variability_analysis
import numpy as np
import copy


#importing the models and mapping files
recon = cobra.io.load_matlab_model('Calbicans_LRhamnosus_EpithelialCells_Interaction/models/Recon3DModel_301.mat')
candida = cobra.io.read_sbml_model('Calbicans_LRhamnosus_EpithelialCells_Interaction/models/Suppl1_candida_albicans_GSMM.xml')
lacto = cobra.io.load_matlab_model('Calbicans_LRhamnosus_EpithelialCells_Interaction/models/Lactobacillus_rhamnosus_LMS2_1.mat')
vmh_rxns = pd.read_table('Calbicans_LRhamnosus_EpithelialCells_Interaction/data/vmh_reactions.tsv')
candida_rxns = pd.read_excel('Calbicans_LRhamnosus_EpithelialCells_Interaction/data/Suppl2_model_info.xlsx')

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
                      

#import mapping file for candida model
candida_met = pd.read_csv('Calbicans_LRhamnosus_EpithelialCells_Interaction/data/Candida_albicans_CORECO_to_AGORA_map.csv',index_col=0)

#metabolome data. calculate ratios for
#EC vs blank (for EC)
#Lr vs blank (for Lr)
#Lr_EC vs EC (for EC)
#Lr_EC vs Lr (for Lr)
#Ca_Lr_EC vs Lr_LC (for Lr and EC)
#Ca vs blank (for Ca)
#Ca_EC vs Ca (for Ca)

#import metabolomics data
metab_data = pd.read_csv('Calbicans_LRhamnosus_EpithelialCells_Interaction/data/metabolomics_data.csv',index_col=0)


#collect headers of metabolome data
blank = ['Blank_6 h','Blank_12 h']
lr_6 = ['Blank_6 h','Lr_6 h', 'Lr_6 h.1', 'Lr_6 h.2', 'Lr_6 h.3', 'Lr_6 h.4']
lr_12 = ['Blank_12 h','Lr_12 h', 'Lr_12 h.1', 'Lr_12 h.2', 'Lr_12 h.3', 'Lr_12 h.4']
ec_6 = ['Blank_6 h','EC_6 h','EC_6 h.1', 'EC_6 h.2', 'EC_6 h.3', 'EC_6 h.4']
ec_12 = ['Blank_12 h','EC_12 h','EC_12 h.1', 'EC_12 h.2', 'EC_12 h.3', 'EC_12 h.4']
ca_ec_6 = ['Blank_6 h','Ca_EC_6 h','Ca_EC_6 h.1', 'Ca_EC_6 h.2', 'Ca_EC_6 h.3', 'Ca_EC_6 h.4']
ca_ec_12 = ['Blank_12 h','Ca_EC_12 h','Ca_EC_12 h.1', 'Ca_EC_12 h.2', 'Ca_EC_12 h.3', 'Ca_EC_12 h.4']
ca_6 = ['Blank_6 h','Ca_6 h', 'Ca_6 h.1', 'Ca_6 h.2', 'Ca_6 h.3', 'Ca_6 h.4']
ca_12 = ['Blank_12 h','Ca_12 h', 'Ca_12 h.1', 'Ca_12 h.2', 'Ca_12 h.3', 'Ca_12 h.4']
lr_ec_6 = ['Blank_6 h','Lr_EC_6 h','Lr_EC_6 h.1','Lr_EC_6 h.2','Lr_EC_6 h.3','Lr_EC_6 h.4']
lr_ec_12 = ['Blank_12 h','Lr_EC_12 h','Lr_EC_12 h.1','Lr_EC_12 h.2','Lr_EC_12 h.3','Lr_EC_12 h.4']
lr_ec_ca_6 = ['Ca_LR_EC_6 h','Ca_LR_EC_6 h.1','Ca_LR_EC_6 h.2','Ca_LR_EC_6 h.3','Ca_LR_EC_6 h.4']
lr_ec_ca_12 = ['Ca_LR_EC_12 h','Ca_LR_EC_12 h.1','Ca_LR_EC_12 h.2','Ca_LR_EC_12 h.3','Ca_LR_EC_12 h.4']


#convert metab_data to influx values. each value divided by the maximum value in the matab_data multiplied to 1000
metab_data.iloc[:,3:] = metab_data.iloc[:,3:].fillna(0)
metab_data.iloc[:,3:] = metab_data.iloc[:,3:]/metab_data.iloc[:,3:].max().max()*-1000


#comparison to blank
ratios_all = list()
row_sum_all = list()
nominators_all = list()
denominators_all = list()
for c in [blank,lr_6,lr_12,ec_6,ec_12,ca_6,ca_12,lr_ec_6,lr_ec_12,ca_ec_6,ca_ec_12]:
    sub_df = metab_data.iloc[:,metab_data.columns.isin(c)]
    sub_df = sub_df.fillna(0)#replace NA values with 0
    row_sum = list(sub_df.sum(axis=1))
    
    nominator = list()
    denominator = list()
    for i in range(sub_df.shape[0]):
        nominator.append(sub_df.iloc[i,1:].mean())
        denominator.append(sub_df.iloc[i,0])
        
    nominators_all.append(nominator)
    denominators_all.append(denominator)
    row_sum_all.append(row_sum)


#prepare vectors for non-blank comparisons
del lr_6[0]
del lr_12[0]
del ca_6[0]
del ca_12[0]
del lr_6[0]
del ec_6[0]
del ec_12[0]
del ca_ec_6[0]
del ca_ec_12[0]
del lr_ec_6[0]
del lr_ec_12[0]

#non-blank comparisons
for c in [ec_6+lr_ec_6, ec_12+lr_ec_12, ca_6+ca_ec_6, ca_12+ca_ec_12, ec_6+ca_ec_6, ec_12+ca_ec_12, lr_6+lr_ec_6, lr_12+lr_ec_12]:
    sub_df = metab_data.iloc[:,metab_data.columns.isin(c)]
    sub_df = sub_df.fillna(0)#replace NA values with 0
    row_sum = list(sub_df.sum(axis=1))
    
    nominator = list()
    denominator = list()
    for i in range(sub_df.shape[0]):
        nominator.append(sub_df.iloc[i,5:].mean())
        denominator.append(sub_df.iloc[i,:5].mean())
        
    nominators_all.append(nominator)
    denominators_all.append(denominator)
    row_sum_all.append(row_sum)




# met mapping
agora_mets = pd.read_table('Calbicans_LRhamnosus_EpithelialCells_Interaction/data/recon-store-metabolites-1.tsv')
manual_match = pd.read_table('Calbicans_LRhamnosus_EpithelialCells_Interaction/data/manual_met_match.csv')
manual_match = list(manual_match['agora_met; name'])

ids = agora_mets[agora_mets['keggId'].isin(list(metab_data['KEGG']))]

agora_ids = list()
id_by_kegg = list()
id_by_pub = list()
id_by_hmdb = list()

for i in range(metab_data.shape[0]):
    if not pd.isna(metab_data['KEGG'][i]) and agora_mets[agora_mets['keggId']==metab_data['KEGG'][i]].shape[0]!=0:
        id_by_kegg.append(agora_mets[agora_mets['keggId']==metab_data['KEGG'][i]]['abbreviation'].iloc[0])
    else:
        id_by_kegg.append('')

    if not pd.isna(metab_data['PUBCHEM'][i]) and agora_mets[agora_mets['pubChemId']==metab_data['PUBCHEM'][i]].shape[0]!=0:
        id_by_pub.append(agora_mets[agora_mets['pubChemId']==metab_data['PUBCHEM'][i]]['abbreviation'].iloc[0])
    else:
        id_by_pub.append('')

    #hmdb modification
    if pd.isna(metab_data['HMDB'][i]):
        id_by_hmdb.append('')
        continue
    else:    
        hmdb_id = 'HMDB00' + metab_data['HMDB'][i][4:]
    
    if not pd.isna(hmdb_id) and agora_mets[agora_mets['hmdb']==hmdb_id].shape[0]!=0:
        id_by_hmdb.append(agora_mets[agora_mets['hmdb']==hmdb_id]['abbreviation'].iloc[0])
    else:
        id_by_hmdb.append('')

#dataframe of mapped metabolites to agora ids by using three databases.
#In case columns are not identical in each row, it should be interrogated
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
        print(metab_data.index[i])
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


#mapping using manually prepared mapping data
for i in manual_match:
    s = i.split(';')
    if s[0] not in agora_ids and s[1][1:] in list(metab_data.index):
        agora_ids[list(metab_data.index).index(s[1][1:])] = s[0]

#continue manual mapping
agora_ids[list(metab_data.index).index("mannonate*")] = "mana"
agora_ids[list(metab_data.index).index("aconitate [cis or trans]")] = "acon_C"
agora_ids[list(metab_data.index).index("arabonate/xylonate")] = "xylnt"
agora_ids[list(metab_data.index).index("diacetylspermidine*")] = "CE1059"
agora_ids[list(metab_data.index).index("glycerophosphoglycerol")] = "g3pg"
agora_ids[list(metab_data.index).index("palmitoyl sphingomyelin (d18:1/16:0)")] = "sphmyln18116_hs"
agora_ids[list(metab_data.index).index("aconitate [cis or trans]")] = "acon_C"

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
kegg_ids = list(metab_data['KEGG'])
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


####
#set up the diet with respect to metabolome data and make list of ex reactions for each model with respect to metabolome data (available compounds in models)


#recon
recon_ex_rxns = list()
for i in recon.reactions:
    if i.id[:2]=='EX':
        if i.id in recon_ex_ids:
            i.bounds=(-1000,1000)
            recon_ex_rxns.append(i)
        else:
            i.bounds = (0,1000)

#lacto
lacto_ex_rxns = list()
for i in lacto.reactions:
    if i.id[:2]=='EX':
        if i.id in lacto_ex_ids:
            i.bounds=(-1000,1000)
            lacto_ex_rxns.append(i)
        else:
            i.bounds = (0,1000)

#candida
candida_ex_rxns = list()
for i in candida.reactions:
    if i.id[:2]=='Ex':
        if i.id in candida_ex_ids:
            i.bounds=(-1000,1000)
            candida_ex_rxns.append(i)
        else:
            i.bounds = (0,1000)


#this function allows influx of metabolites in the model which are present in the ref media like Blank
#in order to perform FBA in pfba_impl function. So, using FBA-derived values and the ratios of respective
#conditions (division of metabolites), the boundaries will be adjusted and FVA will be performed.
def set_diet_ref_based(model,ex,upt):
    for i in model.reactions:
        if i.id[:3]=='EX_' or i.id[:3]=='Ex_':
            if i.id in ex:
                i.bounds=(upt[ex.index(i.id)].mean(),1000)
            else:
                i.bounds = (0,1000)
            

            #if the model is lacto, allow influx of the essentials
            if model.id=='model' and i.id in lacto_essentials:
                i.bounds=(-1,1000)
            if model.id=='Recon3DModel' and i.id in recon_essentials:
                i.bounds=(-1,1000)
            if model.id=='MODEL1604280052' and i.id in candida_essentials:
                i.bounds=(-1,1000)    
    return model

def fva_impl(model,ex_ids,rats,cond_refs,metab_data,rxns_ids_for_fva):

    adp_fva = list()
    ref_fva = list()
    ref_fba = list()

    rxns_for_fva = list()
    for i in rxns_ids_for_fva:
        rxns_for_fva.append(model.reactions.get_by_id(i))
    
    for i in range(len(rats)):

        ####set the model with respect to reference metabolome data to do FBA
        #present = np.array([False]*len(ex_ids))
        #present[metab_data.loc[:,metab_data.columns.isin(cond_refs[i])].sum(axis=1)>0] = True
        present = metab_data.loc[:,metab_data.columns.isin(cond_refs[i])].values
        model = set_diet_ref_based(model,ex_ids,present)
        sol = model.optimize()
        ref_fba.append(sol.objective_value)
        fva = flux_variability_analysis(model,rxns_for_fva)
        ref_fva.append(fva)
        ####
        
        for j in range(len(ex_ids)):
    
            if ex_ids[j] in model.reactions:
                    
                r = model.reactions.get_by_id(ex_ids[j])

                flx = sol.fluxes[ex_ids[j]]
                r.bounds = (rats[i][j],1000)
                
            
        fva = flux_variability_analysis(model,rxns_for_fva)
        adp_fva.append(fva)
    return adp_fva,sol.status,ref_fva,ref_fba


#collect reference media for each cell

refs_for_lacto = [['Blank_6 h'],['Blank_6 h'],['Blank_12 h'],
                  ['EC_6 h', 'EC_6 h.1', 'EC_6 h.2', 'EC_6 h.3', 'EC_6 h.4'],
                  ['EC_12 h', 'EC_12 h.1', 'EC_12 h.2', 'EC_12 h.3', 'EC_12 h.4']]

refs_for_recon = [['Blank_6 h'],['Blank_12 h'],['Blank_6 h'],['Blank_12 h'],
                  ['EC_6 h', 'EC_6 h.1', 'EC_6 h.2', 'EC_6 h.3', 'EC_6 h.4'],
                  ['EC_12 h', 'EC_12 h.1', 'EC_12 h.2', 'EC_12 h.3', 'EC_12 h.4'],
                  ['EC_6 h', 'EC_6 h.1', 'EC_6 h.2', 'EC_6 h.3', 'EC_6 h.4'],
                  ['EC_12 h', 'EC_12 h.1', 'EC_12 h.2', 'EC_12 h.3', 'EC_12 h.4']]

refs_for_candida = [['Blank_6 h'],['Blank_12 h'],
                    ['EC_6 h', 'EC_6 h.1', 'EC_6 h.2', 'EC_6 h.3', 'EC_6 h.4'],
                    ['EC_12 h', 'EC_12 h.1', 'EC_12 h.2', 'EC_12 h.3', 'EC_12 h.4']]


rxns_ids_for_fva = [[i.id for i in lacto.reactions if i.id[:2]!='EX']]
adp_fva_lacto_all_conds = []
ref_fva_lacto_all_conds = []
for condition in rxns_ids_for_fva:
    lacto_copied = lacto.copy()
    adp_fva_lacto,status,ref_fva_lacto,ref_fba_lacto = fva_impl(lacto_copied,lacto_ex_ids,[nominators_all[i] for i in [0,3,4,15,16]],refs_for_lacto,metab_data,condition)
    adp_fva_lacto_all_conds.append(adp_fva_lacto)
    ref_fva_lacto_all_conds.append(ref_fva_lacto)


rxns_ids_for_fva = [[i.id for i in recon.reactions if i.id[:2]!='EX']]
adp_fva_recon_all_conds = []
ref_fva_recon_all_conds = []
for condition in rxns_ids_for_fva:
    recon_copied = recon.copy()
    adp_fva_recon,status,ref_fva_recon,ref_fba_recon = fva_impl(recon_copied,recon_ex_ids,[nominators_all[i] for i in [7,8,9,10,11,12,15,16]],refs_for_recon,metab_data,condition)
    adp_fva_recon_all_conds.append(adp_fva_recon)
    ref_fva_recon_all_conds.append(ref_fva_recon)


rxns_ids_for_fva = [[i.id for i in candida.reactions if i.id[:2]!='Ex']]
adp_fva_candida_all_conds = []
ref_fva_candida_all_conds = []
for condition in rxns_ids_for_fva:
    candida_copied = candida.copy()
    adp_fva_candida,status,ref_fva_candida,ref_fba_candida = fva_impl(candida_copied,candida_ex_ids,[nominators_all[i] for i in [3,4,11,12]],refs_for_candida,metab_data,condition)
    adp_fva_candida_all_conds.append(adp_fva_candida)
    ref_fva_candida_all_conds.append(ref_fva_candida)
    

######
#analyzing the data

#take pathways associated to recon and lacto reactions

def fva2path(model,fva_list,fva_ref,info,key_col,key_rxn):
    
    for i in range(len(fva_list[0])):

        #a bound changes if deviation is more than 10% when both bounds are greater than 0.001
        #or when one is below 0.001 and the other greater than 0.001

        #so first get rid of close to zero values in the dataframes and replace them with zero
        fva_list[0][i][abs(fva_list[0][i])<0.001] = 0
        fva_ref[0][i][abs(fva_ref[0][i])<0.001] = 0
        
        ratio_min = fva_list[0][i]['minimum'] / fva_ref[0][i]['minimum']    
        ratio_max = fva_list[0][i]['maximum'] / fva_ref[0][i]['maximum']

        #drop the NaN values in ratio_min and ratio_max. it removes unchanged bounds (both zero)
        ratio_min = ratio_min.dropna()
        ratio_max = ratio_max.dropna()

        #filter ratios for greater than 10% deviation
        ratio_min = ratio_min[(abs(ratio_min)>1.1).values | (abs(ratio_min)<0.9).values]
        ratio_max = ratio_max[(abs(ratio_max)>1.1).values | (abs(ratio_max)<0.9).values]

        ratio_min = set(ratio_min.index)
        ratio_max = set(ratio_max.index)

        # if lower_bound changed OR upper_bound changed
        diffs = ratio_min | ratio_max
        diffs = list(diffs)

        #compute per pathway

        #associated pathways to all diffs
        paths = list(info[info[key_rxn].isin(diffs)][key_col])

        #unique of pathway
        paths_unq = list(set(paths))

        #count diff reactions per pathway
        num_diff_paths = [paths.count(j) for j in paths_unq]

        #to count size of each pathway in the model
        all_model_rxns = [j.id for j in model.reactions]

        #take all pathways in the model
        paths = list(info[info[key_rxn].isin(all_model_rxns)][key_col])

        num_all_paths = [paths.count(j) for j in paths_unq]

        num_all_paths = np.array(num_all_paths)
        num_diff_paths = np.array(num_diff_paths)
        relative = list(num_diff_paths/num_all_paths)

        if i==0: 
            df = pd.DataFrame({'path_change'+str(i):relative,'path_count':list(num_all_paths)},index=paths_unq)
        else:
            df_to_join = pd.DataFrame({'path_change'+str(i):relative,'path_count':list(num_all_paths)},index=paths_unq)
            df = pd.merge(df,df_to_join,how='outer',left_index=True,right_index=True,on='path_count')
        
    return df

df_recon = fva2path(recon,adp_fva_recon_all_conds,ref_fva_recon_all_conds,vmh_rxns,'subsystem','abbreviation')
df_lacto = fva2path(lacto,adp_fva_lacto_all_conds,ref_fva_lacto_all_conds,vmh_rxns,'subsystem','abbreviation')
df_candida = fva2path(candida,adp_fva_candida_all_conds,ref_fva_candida_all_conds,candida_rxns,'Subsystem','ID')

df_lacto = df_lacto.fillna(0)
df_recon = df_recon.fillna(0)
df_candida = df_candida.fillna(0)


#remove pathways having less than 4 reactions
df_lacto = df_lacto[df_lacto.path_count>3]
df_recon = df_recon[df_recon.path_count>3]
df_candida = df_candida[df_candida.path_count>3]


df_lacto.columns = ['blank 12h vs blank 6h',
                    'ptw_size',
                    'IEC vs blank (6h)',
                    'IEC vs blank (12h)',
                    'Ca+IEC vs IEC (6h)',
                    'Ca+IEC vs IEC (12h)']

df_recon.columns = ['Lr+IEC vs blank (6h)',
                    'ptw_size',
                    'Lr+IEC vs blank (12h)',
                    'Ca+IEC vs blank (6h)',
                    'Ca+IEC vs blank (12h)',
                    'Lr+IEC vs IEC (6h)',
                    'Lr+IEC vs IEC (12h)',
                    'Ca+IEC vs IEC (6h)',
                    'Ca+IEC vs IEC (12h)']

df_candida.columns = ['IEC vs blank (6h)',
                      'ptw_size',
                      'IEC vs blank (12h)',
                      'Lr+IEC vs IEC (6h)',
                      'Lr+IEC vs IEC (12h)']


df_lacto.to_csv('Calbicans_LRhamnosus_EpithelialCells_Interaction/results/lacto_fva_pro_adj_bound.csv')
df_recon.to_csv('Calbicans_LRhamnosus_EpithelialCells_Interaction/results/recon_fva_pro_adj_bound.csv')
df_candida.to_csv('Calbicans_LRhamnosus_EpithelialCells_Interaction/results/candida_fva_adj_bound.csv')

#make a table containing information for all reaction changes of Candida model in EC+Lr vs EC at 12hr
table = pd.concat([ref_fva_candida_all_conds[0][0],adp_fva_candida_all_conds[0][0]],axis=1)
table = pd.concat([table,ref_fva_candida_all_conds[0][1]],axis=1)
table = pd.concat([table,adp_fva_candida_all_conds[0][1]],axis=1)
table = pd.concat([table,adp_fva_candida_all_conds[0][2]],axis=1)
table = pd.concat([table,adp_fva_candida_all_conds[0][3]],axis=1)
candida_rxns = candida_rxns.drop('Unnamed: 0',axis=1)
candida_rxns.index = candida_rxns['ID']
candida_rxns = candida_rxns.drop('ID',axis=1)
table = pd.merge(table,candida_rxns,left_index=True,right_index=True)
table.columns = ['blank6_lb','blank6_ub','EC6_lb','EC6_ub',
                 'blank12_lb','blank12_ub','EC12_lb','EC12_ub',
                 'Lr+EC6_lb','Lr+EC6_ub','Lr+EC12_lb','Lr+EC12_ub',
                 'Name','EC','Equation','Subsystem','Gene association']
table = table.sort_values(by='Subsystem')
table.to_csv('Documents/candilactoepi_hube/added_files_by_me/results/updated_results_19_NOV/all_candida_rxns_metabolome_adj_bounds_complete.csv')
