#%%
import pandas as pd
import os
from Bio import AlignIO
from Bio import SeqIO
import subprocess
from HIV_indel_func import *

bp = os.path.join(os.getcwd(),'data')
path_to_seq = os.path.join(bp, 'eu_env_aa_one_seq_pP.fasta')
path_to_bioNJ = os.path.join(bp, 'bioNJ_env_aa_one_seq_pP_subtypeB.nwk')
path_to_initial_msa = os.path.join(bp,'initial_msa_env_aa_one_seq_pP_subtypeB.fas')
path_to_date_database = os.path.join(bp,'env_date_database.txt')
path_to_msa = os.path.join(bp,'msa_env_aa_one_seq_pP_subtypeB.fas')
path_to_ancestral = os.path.join(bp,'ancestral_env_aa_one_seq_pP_subtypeB')
path_to_nexus_dated = os.path.join(bp, 'initial_msa_env_aa_one_seq_pP_subtypeB.fas.timetree.nex')


path_to_iqtree = os.path.join('/Applications/iqtree-2.0.6-MacOSX/bin/iqtree2')
path_to_indelMaP = os.path.join('/Users/iglh/Desktop/indelMaP-main/src/indelMaP_MSA.py')
path_to_indelMaP_ASR = os.path.join('/Users/iglh/Desktop/indelMaP-main/src/indelMaP_ASR.py')

#%% extract subtypeB and write accession numbers
path_to_seq_B = extract_subtypeB(bp,path_to_seq)
path_to_seq_B_clean = delete_duplicated(bp, path_to_seq_B)
write_accession_no(bp, path_to_seq_B_clean)
#%% delete ambigous characters for ProPIP
path_to_seq_B_not_amb = delete_ambiguous(bp, path_to_seq_B_clean, 'Protein')

#%% run indelMaP to get initial MSA
subprocess.run('python '+path_to_indelMaP+' -s '+path_to_seq_B+' -t '+path_to_bioNJ+' -a Protein -o '+path_to_initial_msa.split('.')[0]+' -q HIVb', shell=True)

#%% extract start and end sites for env regions
start_sites, end_sites = extract_regions(path_to_initial_msa)       

#%% write the date file for IQTree                         
date_file_path=write_date_file(bp, path_to_initial_msa, path_to_date_database)    
#%% estimate tip dated phylogenetic tree for the whole sequence
subprocess.run(path_to_iqtree+ ' -s '+path_to_initial_msa+' -t '+path_to_bioNJ+' -m HIVb '+
               '-o B.FR.1983.HXB2-LAI-IIIB-BRU.K03455.LAI'+
               ' --date '+date_file_path+' --date-options "-l 0.0"'+' --redo', shell=True)

#%% run indelMaP to get a MSA
path_to_tip_dated_tree = os.path.join(path_to_initial_msa+'.timetree.nwk')
subprocess.run('python '+path_to_indelMaP+' -s '+path_to_seq_B+' -t '+path_to_tip_dated_tree+' -a Protein -o '+path_to_msa.split('.')[0]+' -q HIVb', shell=True)

#%% run indelMaP_ASR to estimate ancestral reconstruction and infer evolutionary events
subprocess.run('python '
               +path_to_indelMaP_ASR+' -m '+path_to_msa+' -t '+path_to_tip_dated_tree+' -a Protein -o '+path_to_ancestral+' -q HIVb', shell=True)

# %%
events_df = write_evolutionary_events(bp, path_to_ancestral, path_to_nexus_dated)
# %%
path_to_events_all = os.path.join(bp, 'evolutionary_events.csv')
path_to_msa_events_all = os.path.join(bp, 'ancestral_env_aa_one_seq_pP_subtypeB_all_evolutionary_events.fas')
path_to_newick_dated_w_ROOOT = os.path.join(bp, 'tree_env_aa_root.timetree.nwk')
calculate_indel_rates(bp, path_to_events_all, path_to_newick_dated_w_ROOOT, path_to_msa_events_all)
        
#%%
rates_df = pd.read_csv(os.path.join(bp,'evolutionary_rates.csv'))
plot_indel_rates_gp120(bp, rates_df)
# %%

rates_df = pd.read_csv(os.path.join(bp,'evolutionary_rates.csv'))
plot_indel_rates_gp120_variable(bp, rates_df)

# %%
path_to_events_all = os.path.join(bp, 'evolutionary_events.csv')
path_to_newick_dated_w_ROOT = os.path.join(bp, 'tree_env_aa_root.timetree.nwk')
plot_indel_lengths(bp,path_to_events_all,path_to_newick_dated_w_ROOT)

# %%
path_to_newick_dated_w_ROOT = os.path.join(bp, 'tree_env_aa_root.timetree.nwk')
path_to_msa_events_all = os.path.join(bp, 'ancestral_env_aa_one_seq_pP_subtypeB_all_evolutionary_events.fas')
path_to_events_all = os.path.join(bp, 'evolutionary_events.csv')
calc_plot_rates_indel_event(bp, path_to_newick_dated_w_ROOT,path_to_events_all,path_to_msa_events_all)

# %%
