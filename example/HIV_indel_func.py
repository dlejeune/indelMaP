import os
import re
import pandas as pd
from Bio import SeqIO, AlignIO
from ete3 import PhyloNode
from Bio.Phylo import convert
from researchpy import summary_cont
import matplotlib.colors as mc
import colorsys
import seaborn_image as isns
from matplotlib.patches import PathPatch
import seaborn as sns
import matplotlib.pyplot as plt


def extract_subtypeB(bp,path_to_seq):
    seq = SeqIO.parse(open(path_to_seq),'fasta')
    with open(os.path.join(bp, 'seq_env_subtypeB.fas'), 'w') as f:
        for record in seq:
            label = record.id.split('.')
            if record.id.split('.')[0] == 'B' and record.id.split('.')[2]!='-':
                print('>',record.id, file=f)
                print(record.seq, file=f)
    return os.path.join(bp, 'seq_env_subtypeB.fas')

def delete_duplicated(bp, path_to_seq):
    seq = SeqIO.parse(open(path_to_seq),'fasta')
    seq_comp = SeqIO.parse(open(path_to_seq),'fasta')
    out = os.path.join(bp, 'seq_env_subtypeB_clean.fas') 
    with open(out, 'w') as f:
        for record in seq:
            dup = False
            for record_comp in seq_comp:
                if str(record.name) != str(record_comp.name) and str(record.seq) == str(record_comp.seq):
                    print("duplicated")
                    dup = True
            if not dup:
                print(">"+record.name, file=f)
                print(record.seq, file=f)
    return out


def write_accession_no(bp,path_to_seq_B):
    seq = SeqIO.parse(open(path_to_seq_B),'fasta')
    with open(os.path.join(bp,'accession_file.txt'),'w') as f:
        for record in seq:
            print(record.id.split('.')[-2], file=f)

def delete_ambiguous(bp, path_to_seq_B,alphabet):
    if alphabet == 'Protein':
        characters = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M',
                      'F', 'P', 'S', 'T', 'W', 'Y', 'V']
    elif alphabet == 'DNA':
        characters = ['A', 'C', 'G', 'T']
    seq_B = SeqIO.parse(open(path_to_seq_B),'fasta')
    with open(os.path.join(bp,'seq_env_subtypeB_not_ambiguous.fas'), 'w') as f:
        for record in seq_B:
            new_seq=''
            for site in range(len(record.seq)):
                if record.seq[site] in characters:
                    new_seq+=record.seq[site]
            record.seq = new_seq
            print('>'+record.id, file=f)
            print(record.seq, file=f)
    return os.path.join(bp,'seq_env_subtypeB_not_ambiguous.fas')

def extract_regions(path_to_msa):
    msa = AlignIO.read(path_to_msa, format='fasta')
    for record in msa._records:
        if re.search(r'K03455', record.id):
            count = 0
            for i in range(len(record.seq)):
                if record.seq[i] != '-' and record.seq[i] != '*':
                    if count == 0:
                        start_signal_peptide=i
                    if count == 29:
                        end_signal_peptide=i
                        print('End of env signal peptide', i)
                    elif count == 30:
                        start_c1=i
                    elif count == 128:
                        end_c1=i
                        print('C1 from', start_c1, 'to', end_c1)
                    elif count == 129:
                        start_v1=i
                    elif count == 157:
                        end_v1=i
                        print('V1 from', start_v1, 'to', end_v1)
                    elif count == 156:
                        start_v2=i
                    elif count == 196:
                        end_v2=i
                        print('V2 from', start_v2, 'to', end_v2)
                    elif count == 197:
                        start_c2=i
                    elif count == 293:
                        end_c2=i
                        print('C2 from', start_c2, 'to', end_c2)
                    elif count == 294:
                        start_v3=i
                    elif count == 331:
                        end_v3=i
                        print('V3 from', start_v3, 'to', end_v3)
                    elif count == 332:
                        start_c3=i
                    elif count == 382:
                        end_c3=i
                        print('C3 from', start_c3, 'to', end_c3)
                    elif count == 383:
                        start_v4=i
                    elif count == 418:
                        end_v4=i
                        print('V4 from', start_v4, 'to', end_v4)
                    elif count == 419:
                        start_c4=i
                    elif count == 458:
                        end_c4=i
                        print('C4 from', start_c4, 'to', end_c4)
                    elif count == 459:
                        start_v5=i
                    elif count == 469:
                        end_v5=i
                        print('V5 from', start_v5, 'to', end_v5)
                    elif count == 471:
                        start_c5=i
                    elif count == 510:
                        end_c5=i
                        print('C5 from', start_c5, 'to', end_c5) 
                    elif count == 510:
                        print('End of gp120', i)
                    elif count == 511:
                        start_gp41=i    
                    elif count == 855:
                        end_gp41=i
                        print('gp41 from', start_gp41, 'to', end_gp41,
                        'or end of alignment', len(record.seq))
                    count += 1
    start_signal_peptide_c1 = end_signal_peptide+1
    end_signal_peptide_c1 = start_c1-1
    start_c1_v1=end_c1+1
    end_c1_v1=start_v1-1
    start_v1_v2=end_v1+1
    end_v1_v2=start_v2-1
    start_v2_c2=end_v2+1
    end_v2_c2=start_c2-1
    start_c2_v3=end_c2+1
    end_c2_v3=start_v3-1
    start_v3_c3=end_v3+1
    end_v3_c3=start_c3-1
    start_c3_v4=end_c3+1
    end_c3_v4=start_v4-1
    start_v4_c4=end_v4+1
    end_v4_c4=start_c4-1
    start_c4_v5=end_c4+1
    end_c4_v5=start_v5-1
    start_v5_c5=end_v5+1
    end_v5_c5=start_c5-1
    start_c5_gp41=end_c5+1
    end_c5_gp41=start_gp41-1
    start_gp41_end=end_gp41+1
    end_gp41_end=len(record.seq)-1

    start_sites = [start_signal_peptide, start_signal_peptide_c1, start_c1, start_c1_v1, start_v1, start_v1_v2, 
                start_v2, start_v2_c2, start_c2, start_c2_v3, start_v3, start_v3_c3, start_c3, 
                start_c3_v4, start_v4, start_v4_c4, start_c4, start_c4_v5, start_v5, start_v5_c5,
                start_c5, start_c5_gp41, start_gp41, start_gp41_end]

    end_sites =  [end_signal_peptide, end_signal_peptide_c1, end_c1, end_c1_v1, end_v1, end_v1_v2, 
                end_v2, end_v2_c2, end_c2, end_c2_v3, end_v3, end_v3_c3, end_c3, 
                end_c3_v4, end_v4, end_v4_c4, end_c4, end_c4_v5, end_v5, end_v5_c5, end_c5, end_c5_gp41, end_gp41, end_gp41_end]  
                 
    return start_sites, end_sites

def define_region_env(start_sites, end_sites, site):
    regions = ['env signal peptide', 'between signal peptide and C1', 'gp120 - C1', 'gp120 - between C1 and V1', 
               'gp120 - V1', 'gp120 - between V1 and V2', 'gp120 - V2', 'gp120 - between V2 and C2', 'gp120 - C2',
               'gp120 - between C2 and V3', 'gp120 - V3', 'gp120 - between V3 and C3', 'gp120 - C3', 'gp120 - between C3 and V4', 
               'gp120 - V4', 'gp120 - between V4 and C4', 
               'gp120 - C4', 'gp120 - between C4 and V5', 'gp120 - V5', 'gp120 - between V5 and C5', 'gp120 - C5', 
               'between C5 and gp41', 'gp41', 'between gp41 and end']
    
    for i in range(len(regions)):
        if site >= start_sites[i] and site <= end_sites[i]:
            region = regions[i]
    
    return region

def write_date_file(bp, path_to_msa, date_database):

    months = {'Jan':'01','Feb':'02','Mar':'03','Apr':'04','May':'05','Jun':'06',
              'Jul':'07','Aug':'08','Sep':'09','Oct':'10','Nov':'11','Dec':'12'}
    with open(date_database,'r') as data:
        acc_date = {}
        d_tf = True

        for line in data:
            if "ACCESSION" in line:
                if d_tf==False:
                    acc_date[accession]='no date'
                accession = line.split(' ')[3].strip("\n")
                d_tf=False
            elif '/collection_date' in line:
                d = line.split('=')[-1].strip('\n').strip('"')
                if len(d.split('-')) == 3:
                    d = d.split('-')[2]+'-'+months[d.split('-')[1]]+'-'+d.split('-')[0]
                elif len(d.split('-')) == 2:
                    d = d.split('-')[1]+'-'+months[d.split('-')[0]]
                else:
                    d='no date'

                acc_date[accession]=d
                d_tf=True
    msa = AlignIO.read(path_to_msa, format='fasta')
    with open(os.path.join(bp, 'date_file_env.txt'),'w') as f:
        for record in msa._records:
            if record.id.split('.')[-2] not in acc_date.keys():
                if record.id.split('.')[2] != '-':
                    print(record.id, record.id.split('.')[2], file=f)
            elif acc_date[record.id.split('.')[-2]] == 'no date':
                if record.id.split('.')[2] != '-':
                    print(record.id, record.id.split('.')[2], file=f)
            else:
                print(record.id, acc_date[record.id.split('.')[-2]], file=f )

    return os.path.join(bp, 'date_file_env.txt')

def write_evolutionary_events(bp, path_to_ancestral, path_to_nexus_dated):
    
    path_to_events_leaves = os.path.join(bp,  path_to_ancestral+'_leaves_evolutionary_events.fas')
    path_to_events_internal = os.path.join(bp,  path_to_ancestral+'_internal_evolutionary_events.fas')
    path_to_events_all = os.path.join(bp, path_to_ancestral+'_all_evolutionary_events.fas')
    path_to_newick_dated_with_years = os.path.join(bp, 'msa_env.fas.timetree.nwk')
    path_to_newick_dated = os.path.join(bp, 'tree_env_aa.timetree.nwk')
    path_to_newick_dated_w_ROOT = os.path.join(bp, 'tree_env_aa_root.timetree.nwk')

    events_leaves = AlignIO.read(path_to_events_leaves, format='fasta')
    events_internal = AlignIO.read(path_to_events_internal, format='fasta')

    with open(path_to_events_all, 'w') as f:
        for record in events_internal:
            print('>'+record.id, file=f)
            print(record.seq, file=f)
        for record in events_leaves:
            print('>'+record.id, file=f)
            print(record.seq, file=f)

    convert(path_to_nexus_dated, "nexus", path_to_newick_dated_with_years, "newick")

    f = open(path_to_newick_dated_with_years,'r')
    if os.path.exists(path_to_newick_dated):
        os.remove(path_to_newick_dated)
    f1 = open(path_to_newick_dated, 'w')
    for line in f.readlines():
        regex = re.compile(r'\[(.*?)\]\]')
        match = regex.search(line)
        if match:
            line = re.sub(r'\[(.*?)\]\]', '', line)
        print(line, file=f1)       
    f.close()
    f1.close()

    tree_dated =PhyloNode(newick=path_to_newick_dated, format=1)
    no_internal = len(tree_dated) 
    for node in tree_dated.traverse('preorder'):
        if node.name == '':
            if not node.is_leaf():
                node.name = 'N' + str(no_internal)
                no_internal += 1

    tree_dated.write(format=1, outfile=path_to_newick_dated)
    f = open(path_to_newick_dated, 'r')
    f1 = open(path_to_newick_dated_w_ROOT, 'w')
    for line in f.readlines():
        regex = re.compile(r';')
        match = regex.search(line)
        if match:
            line = re.sub(r';', 'ROOT;', line)
        print(line, file=f1)       

    tree_dated = PhyloNode(newick=path_to_newick_dated_w_ROOT, alignment=path_to_events_all, format=1)

    start_sites, end_sites = extract_regions(path_to_events_all)
    events_df = pd.DataFrame(columns=['parent', 'child', 'time in years', 
                                 'event', 'character', 'region', 
                                 'site number in alignment', 
                                 'position in alignment',
                                 'position w/o placeholders'])

    for node in tree_dated.traverse('preorder'):
        if not node.is_root():
            print(node.name)
            pos_wo_placeholders = 0
            for site in range(len(node.sequence)):
                char = node.sequence[site]
                parent = node.up.name
                char_up = node.up.sequence[site]
                time = node.dist
                region = define_region_env(start_sites, end_sites, site)
                site_no = site+1
                if char != '*':
                    pos_wo_placeholders += 1
                if char.islower() and char_up == '*':
                    event = 'insertion'
                    events_df.loc[len(events_df)] = [parent, node.name, time, 
                                                 event, char, region, 
                                                 site_no, site, pos_wo_placeholders]
                elif char == '-' and char_up != '-':
                    event = 'deletion' 
                    events_df.loc[len(events_df)] = [parent, node.name, time, 
                                                 event, char, region, 
                                                 site_no, site, pos_wo_placeholders]
                elif char not in ['*', '-'] and char.upper() != char_up.upper() and char_up not in ['*', '-']:
                    event = 'substitution'
                    char = char_up+' -> '+char
                    events_df.loc[len(events_df)] = [parent, node.name, time, 
                                                 event, char, region, 
                                                 site_no, site, pos_wo_placeholders]

    events_df.to_csv(os.path.join(bp,'evolutionary_events.csv'))
    return events_df

def calculate_indel_rates(bp,path_to_events_all, path_to_newick_dated_w_ROOT, path_to_msa_events_all):
    tree_string = [line for line in open(path_to_newick_dated_w_ROOT).readlines()]
    bl = re.findall(r":[-+]?(?:\d*\.*\d+)",tree_string[0])
    bl_float = [float(item.split(':')[1]) for item in bl]
    min_bl = min(bl_float)
    events_df = pd.read_csv(path_to_events_all)
    tree_dated = PhyloNode(newick=path_to_newick_dated_w_ROOT, alignment=path_to_msa_events_all, format=1)
    start_sites, end_sites = extract_regions(path_to_msa_events_all)

    rates_df = pd.DataFrame(columns=['parent', 'child', 'tMRCA', 'env region',
                                    'number of AA in region',
                                    'number of subsitutions', 'substitution/AA', 'substitution/AA/year',
                                    'number of deletions', 'deletion/AA', 'deletion/AA/year',
                                    'number of insertions', 'insertion/AA', 'insertion/AA/year',
                                    'number of insertions - deletions', '(insertions-deletions)/AA', 
                                    '(insertions-deletions)/AA/year'])

    regs = ['signal peptide','gp120 - C1','gp120 - V1','gp120 - V2','gp120 - C2','gp120 - V3','gp120 - C3','gp120 - V4','gp120 - C4','gp120 - V5','gp120 - C5','gp41']
    lengths_ref_regs = [0]*len(regs)
    for node in events_df['child'].unique():
        sub_df = events_df[events_df['child'] == node]
        parent = sub_df['parent'].unique()[0]
        tMRCA =  sub_df['time in years'].unique()[0]
        if min_bl >= 1.0:
            min_bl = 0.0
        if tMRCA>min_bl:
            seq = tree_dated.search_nodes(name=node)[0].sequence
            for i in range(len(regs)):
                start_reg = start_sites[i*2]
                end_reg = end_sites[i*2]
                length_reg = 0
                coverage = 0
                for site in range(start_reg,end_reg+1):
                    if seq[site] !='*':
                        length_reg+=1
                    if seq[site] != '-' and seq[site] != '*':
                        coverage+=1
                if re.search(r'K03455', node):
                    lengths_ref_regs[i]=length_reg
                if coverage >= lengths_ref_regs[i]*0.1 and not re.search(r'K03455', node):
                    subs = 0
                    ins = 0
                    dels = 0
                    for row in range(len(sub_df)):
                        if regs[i] in sub_df['region'].iloc[row]:
                            if sub_df['event'].iloc[row] == 'substitution':
                                subs += 1
                            if sub_df['event'].iloc[row] == 'insertion':
                                ins += 1
                            if sub_df['event'].iloc[row] == 'deletion':
                                dels += 1
                    
                    rates_df.loc[len(rates_df)] = [parent, node, tMRCA, regs[i], length_reg, 
                                                subs, subs/length_reg, (subs/length_reg)/tMRCA,
                                                dels, dels/length_reg, (dels/length_reg)/tMRCA,
                                                ins, ins/length_reg, (ins/length_reg)/tMRCA,
                                                (ins-dels), (ins-dels)/length_reg, ((ins-dels)/length_reg)/tMRCA]

    rates_df.to_csv(os.path.join(bp,'evolutionary_rates.csv'))

def lighten_color(color, amount=0.8):  
    # --------------------- SOURCE: @IanHincks ---------------------
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])


def plot_indel_rates_gp120(bp, rates_df):
    print(summary_cont(rates_df.groupby('env region')['substitution/AA/year']))
    print(summary_cont(rates_df.groupby('env region')['insertion/AA/year']))
    print(summary_cont(rates_df.groupby('env region')['deletion/AA/year']))
    # rates_df = pd.read_csv(os.path.join(bp, 'evolutionary_rates.csv')))
    rates_per_year = rates_df[['env region','substitution/AA/year','insertion/AA/year', 'deletion/AA/year']]

    rpy_melt = pd.melt(rates_per_year, id_vars=['env region'], value_vars=['insertion/AA/year', 'deletion/AA/year'],
                    var_name='indel event', value_name='events/AA/year')
    #, 'substitution/AA/year'
    def split(x):
        return x.split('/')[0]
    def split_(x):
        return x.split('-')[1].lstrip()

    rpy_melt['indel event'] = rpy_melt['indel event'].apply(split) 
    rpy_gp120 = rpy_melt.loc[rpy_melt['env region'].isin(['gp120 - V1', 'gp120 - V2', 'gp120 - V3', 'gp120 - V4', 'gp120 - V5',
                                                        'gp120 - C1', 'gp120 - C2', 'gp120 - C3', 'gp120 - C4', 'gp120 - C5'])]
    rpy_gp120['env region'] = ['gp120']*len(rpy_gp120)
    rpy_melt = rpy_melt.loc[rpy_melt['env region'].isin(['gp120 - V1', 'gp120 - V2', 'gp120 - V3', 'gp120 - V4', 'gp120 - V5',
                                                        'gp120 - C1', 'gp120 - C2', 'gp120 - C3', 'gp120 - C4', 'gp120 - C5'])]
    rpy_melt['env region'] = rpy_melt['env region'].apply(split_)
    rpy_melt = pd.concat([rpy_melt,rpy_gp120], axis=0)

    isns.set_context('paper',fontfamily='Verdana')
    isns.set_context("paper", rc={"font.size":18,"axes.titlesize":18,"axes.labelsize":18})
    g = sns.FacetGrid(rpy_melt, height=6*(1.065), aspect=1.7, despine=False)
    g.map(sns.barplot, x="events/AA/year", y="env region", hue='indel event', order=(['gp120', 'C1','V1','V2','C2','V3','C3','V4','C4','V5','C5']),
        #['gp120 - C1','gp120 - V1','gp120 - V2','gp120 - C2','gp120 - V3','gp120 - C3','gp120 - V4','gp120 - C4','gp120 - V5','gp120 - C5']
        errorbar=('ci', 95), orient='h', linewidth=0.4,
                capsize=.1, palette={'insertion':lighten_color('#00a0d1'), 'deletion':lighten_color('#fa5a0f'), 'substitution':'red'},
                data=rpy_melt, width=0.9)

    for ax in g.axes.flat:
        
        ax.set_axisbelow(True)
        ax.grid(True, which='both', axis='x', zorder=-1, linestyle='dashed',
                linewidth=0.55)
        plt.setp(ax.get_yticklabels(),  ha='center', fontsize=16)
        # rotation=90,
        plt.setp(ax.get_xticklabels(), fontsize=16)
        ax.tick_params(axis='y', direction='out', pad=60)
        ax.tick_params(axis=u'both', which=u'both',length=0)
    plt.legend()

    g.savefig(os.path.join(bp,'env_protein_indel_rates_HIV.png'), dpi=600)

def plot_indel_rates_gp120_variable(bp, rates_df):
   
    rates_per_year = rates_df[['env region','substitution/AA/year','insertion/AA/year', 'deletion/AA/year']]
    rpy_melt = pd.melt(rates_per_year, id_vars=['env region'], value_vars=['insertion/AA/year', 'deletion/AA/year'],
                    var_name='indel event', value_name='events/AA/year')
    #, 'substitution/AA/year'
    def split(x):
        return x.split('/')[0]
    def split_(x):
        return x.split('-')[1].lstrip()

    rpy_melt['indel event'] = rpy_melt['indel event'].apply(split) 
    rpy_gp120 = rpy_melt.loc[rpy_melt['env region'].isin(['gp120 - V1', 'gp120 - V2', 'gp120 - V3', 'gp120 - V4', 'gp120 - V5',
                                                        'gp120 - C1', 'gp120 - C2', 'gp120 - C3', 'gp120 - C4', 'gp120 - C5'])]
    rpy_gp120['env region'] = ['gp120']*len(rpy_gp120)
    rpy_melt = rpy_melt.loc[rpy_melt['env region'].isin(['gp120 - V1', 'gp120 - V2', 'gp120 - V3', 'gp120 - V4', 'gp120 - V5',
                                                        'gp120 - C1', 'gp120 - C2', 'gp120 - C3', 'gp120 - C4', 'gp120 - C5'])]
    rpy_melt['env region'] = rpy_melt['env region'].apply(split_)
    rpy_melt = pd.concat([rpy_melt,rpy_gp120], axis=0)
    rpy_melt['protein region'] = rpy_melt['env region']
    rpy_melt['sites/AA/year'] = rpy_melt['events/AA/year']
    isns.set_context('paper',fontfamily='Times New Roman')
    sns.set_context("paper", rc={"font.size":10,"axes.titlesize":10,"axes.labelsize":10})  
    g = sns.FacetGrid(rpy_melt, height=1.9, aspect=1.7, despine=False)
    g.map(sns.barplot, x="sites/AA/year", y="protein region", hue='indel event', order=(['gp120', 'V1','V2','V3','V4','V5']),
        #['gp120 - C1','gp120 - V1','gp120 - V2','gp120 - C2','gp120 - V3','gp120 - C3','gp120 - V4','gp120 - C4','gp120 - V5','gp120 - C5']
        errorbar=('ci', 95), orient='h', linewidth=0.2, errwidth=0.5, 
                capsize=.1, palette={'insertion':'#00539C', 'deletion':lighten_color('#EEA47F',1.5), 'substitution':'red'},
                data=rpy_melt, width=0.8)

    for ax in g.axes.flat:
        sns.set_context("paper", rc={"font.size":10,"axes.titlesize":10,"axes.labelsize":10})
        ax.set_axisbelow(True)
        ax.grid(True, which='both', axis='x', zorder=-1, linestyle='dashed',
                linewidth=0.55)
        plt.setp(ax.get_yticklabels(),  ha='center', fontsize=8)
        # rotation=90,
        plt.setp(ax.get_xticklabels(), fontsize=8)
        ax.tick_params(axis='y', direction='out', pad=13)
        ax.tick_params(axis=u'both', which=u'both',length=0)
    # plt.legend(fontsize=8, bbox_to_anchor=(-0.5,3.7))
    g.savefig(os.path.join(bp,'env_protein_indel_rates_HIV_gp120_variable.png'), dpi=400)

def plot_indel_lengths(bp,path_to_events_all, path_to_newick_dated_w_ROOT):

    events_df = pd.read_csv(path_to_events_all)
    tree_string = [line for line in open(path_to_newick_dated_w_ROOT).readlines()]
    bl = re.findall(r":[-+]?(?:\d*\.*\d+)",tree_string[0])
    bl_float = [float(item.split(':')[1]) for item in bl]
    min_bl = min(bl_float)
    length_df = pd.DataFrame(columns=['parent', 'child', 'tMRCA', 'env region','event', 'length'])
    regs = ['env signal peptide','gp120 - C1','gp120 - V1','gp120 - V2','gp120 - C2','gp120 - V3','gp120 - C3','gp120 - V4','gp120 - C4','gp120 - V5','gp120 - C5','gp41']
    length = 1

    for row in range(1,len(events_df)):
        node = events_df['child'].iloc[row]
        prev_node = events_df['child'].iloc[row-1]
        region = events_df['region'].iloc[row]
        prev_region =  events_df['region'].iloc[row-1]
        prev_parent = events_df['parent'].iloc[row-1]
        prev_tMRCA = events_df['time in years'].iloc[row-1]
        event = events_df['event'].iloc[row]
        pos = events_df['position w/o placeholders'].iloc[row]
        prev_event = events_df['event'].iloc[row-1]
        prev_pos = events_df['position w/o placeholders'].iloc[row-1]
        if event == prev_event and pos == prev_pos+1 and node == prev_node and region == prev_region:
            length += 1
            pass
        else:
            if prev_tMRCA > min_bl:
                length_df.loc[len(length_df)] = [prev_parent, prev_node, prev_tMRCA, prev_region, prev_event, length]
                length = 1

    length_df.to_csv(os.path.join(bp, 'indel_lengths.csv'))

    only_length = length_df[['env region','event','length']]
    print(summary_cont(only_length.groupby(['env region', 'event'])['length']))
    def split(x):
        return x.split('/')[0]
    def split_(x):
        return x.split('-')[1].lstrip()

    rpy_gp120 = only_length.loc[only_length['env region'].isin(['gp120 - V1', 'gp120 - V2', 'gp120 - V3', 'gp120 - V4', 'gp120 - V5',
                                                            'gp120 - C1', 'gp120 - C2', 'gp120 - C3', 'gp120 - C4', 'gp120 - C5'])]
    rpy_gp120['env region'] = ['gp120']*len(rpy_gp120)
    only_length = only_length.loc[only_length['env region'].isin(['gp120 - V1', 'gp120 - V2', 'gp120 - V3', 'gp120 - V4', 'gp120 - V5'])]
    only_length['env region'] = only_length['env region'].apply(split_)
    only_length = pd.concat([only_length,rpy_gp120], axis=0, ignore_index=True)
    mean_length = only_length.groupby(['env region', 'event']).length.agg(['mean', 'median'])
    isns.set_context('paper',fontfamily='Times New Roman')
    isns.set_context("paper", rc={"font.size":8,"axes.titlesize":8,"axes.labelsize":8})
    only_length_indel = only_length[only_length['event'].isin(['insertion','deletion'])]

    g = sns.displot(data=only_length_indel, x='length', hue='event', multiple='stack', col="env region", 
                    discrete=True, facet_kws={'sharey': False, 'sharex': True, 'despine': False, 'margin_titles':False,'subplot_kws':{'title':None}}, 
                    height=1.7, aspect=0.9, col_order=(['gp120', 'V1','V2','V3','V4','V5']), kind='hist',
                    palette={'insertion':'#00539C', 'deletion':lighten_color('#EEA47F',1.5)}, col_wrap=3)

    g.tight_layout()
    # plt.subplots_adjust(hspace=0.03, wspace=0.02)
    for ax in g.axes.flat:
        isns.set_context('paper',fontfamily='Times New Roman')
        sns.set_context("paper", rc={"font.size":8,"axes.titlesize":8,"axes.labelsize":8})
        ax.set_xlim(0,12)
        # ax.set_ylim(0,1800)
        ax.set_xticks(range(2,12,2))
        # ax.set_yticks([0, 500,1500])
        # ax.set_yticks([0,1000],minor=True)
        spec = ax.get_title().split(' = ')[1]
        # select the data for the species
        data = mean_length.loc[spec, :]
        # print data as needed or comment out
        # print(data)
        # plot the lines
        ax.axvline(x=data.loc['insertion']['median'], c='#00539C', ls='-', lw=.5)
        ax.axvline(x=data.loc['deletion']['median'], c=lighten_color('#EEA47F',1.7), ls='--', lw=.5)

        ax.set_axisbelow(True)
        # ax.grid(True, which='both', axis='y', zorder=-1, linestyle='dashed',
        #         linewidth=0.55)
        plt.setp(ax.get_yticklabels(), rotation=90, ha='center', fontsize=8)
        plt.setp(ax.get_xticklabels(), fontsize=8)
        ax.tick_params(axis='y', direction='out', pad=5)
        ax.tick_params(axis=u'both', which=u'both',length=0)
        ax.xaxis.set_label_coords(0.5, -0.2)
        ax.yaxis.set_label_coords(-0.18,0.5)
            
    g.savefig(os.path.join(bp,'Indel_length_HIV.png'), dpi=400)



def calc_plot_rates_indel_event(bp, path_to_newick_dated_w_ROOT,path_to_events_all,path_to_msa_events_all):
    tree_string = [line for line in open(path_to_newick_dated_w_ROOT).readlines()]
    bl = re.findall(r":[-+]?(?:\d*\.*\d+)",tree_string[0])
    bl_float = [float(item.split(':')[1]) for item in bl]
    min_bl = min(bl_float)
    events_df = pd.read_csv(path_to_events_all)
    tree_dated = PhyloNode(newick=path_to_newick_dated_w_ROOT, alignment=path_to_msa_events_all, format=1)
    start_sites, end_sites = extract_regions(path_to_msa_events_all)

    rates_df = pd.DataFrame(columns=['parent', 'child', 'tMRCA', 'env region',
                                    'number of AA in region',
                                    'number of subsitutions', 'substitution/AA', 'substitution/AA/year',
                                    'number of deletions', 'deletion/AA', 'deletion/AA/year',
                                    'number of insertions', 'insertion/AA', 'insertion/AA/year',
                                    'number of insertions - deletions', '(insertions-deletions)/AA', 
                                    '(insertions-deletions)/AA/year'])

    regs = ['signal peptide','gp120 - C1','gp120 - V1','gp120 - V2','gp120 - C2','gp120 - V3','gp120 - C3','gp120 - V4','gp120 - C4','gp120 - V5','gp120 - C5','gp41']
    lengths_ref_regs = [0]*len(regs)
    for node in events_df['child'].unique():
        sub_df = events_df[events_df['child'] == node]
        parent = sub_df['parent'].unique()[0]
        tMRCA =  sub_df['time in years'].unique()[0]
        if min_bl >= 1.0:
            min_bl = 0.0
        if tMRCA>min_bl:
            seq = tree_dated.search_nodes(name=node)[0].sequence
            for i in range(len(regs)):
                start_reg = start_sites[i*2]
                end_reg = end_sites[i*2]
                length_reg = 0
                coverage = 0
                for site in range(start_reg,end_reg+1):
                    if seq[site] !='*':
                        length_reg+=1
                    if seq[site] != '-' and seq[site] != '*':
                        coverage+=1
                if re.search(r'K03455', node):
                    lengths_ref_regs[i]=length_reg
                if coverage >= lengths_ref_regs[i]*0.1 and not re.search(r'K03455', node):
                    subs = 0
                    ins = 0
                    dels = 0
                    for row in range(len(sub_df)):
                        if regs[i] in sub_df['region'].iloc[row]:
                            if sub_df['event'].iloc[row] == 'substitution':
                                subs += 1
                        
                            if sub_df['event'].iloc[row] == 'insertion' and (sub_df['event'].iloc[row-1] != 'insertion' or sub_df['position w/o placeholders'].iloc[row]-1!=sub_df['position w/o placeholders'].iloc[row-1] or sub_df['region'].iloc[row] !=sub_df['region'].iloc[row-1]):
                                ins += 1
                            
                            if sub_df['event'].iloc[row] == 'deletion' and (sub_df['event'].iloc[row-1] != 'deletion' or sub_df['position w/o placeholders'].iloc[row]-1!=sub_df['position w/o placeholders'].iloc[row-1] or sub_df['region'].iloc[row] != sub_df['region'].iloc[row-1]):
                                dels += 1
                    
                    rates_df.loc[len(rates_df)] = [parent, node, tMRCA, regs[i], length_reg, 
                                                subs, subs/length_reg, (subs/length_reg)/tMRCA,
                                                dels, dels/length_reg, (dels/length_reg)/tMRCA,
                                                ins, ins/length_reg, (ins/length_reg)/tMRCA,
                                                (ins-dels), (ins-dels)/length_reg, ((ins-dels)/length_reg)/tMRCA]

    # rates_df.to_csv(os.path.join(bp,'evolutionary_rates.csv'))
    print(summary_cont(rates_df.groupby('env region')['insertion/AA/year']))
    print(summary_cont(rates_df.groupby('env region')['deletion/AA/year']))
    rates_per_year = rates_df[['env region','substitution/AA/year','insertion/AA/year', 'deletion/AA/year']]
    rpy_melt = pd.melt(rates_per_year, id_vars=['env region'], value_vars=['insertion/AA/year', 'deletion/AA/year'],
                    var_name='indel event', value_name='events/AA/year')
    #, 'substitution/AA/year'
    def split(x):
        return x.split('/')[0]
    def split_(x):
        return x.split('-')[1].lstrip()

    rpy_melt['indel event'] = rpy_melt['indel event'].apply(split) 
    rpy_gp120 = rpy_melt.loc[rpy_melt['env region'].isin(['gp120 - V1', 'gp120 - V2', 'gp120 - V3', 'gp120 - V4', 'gp120 - V5',
                                                        'gp120 - C1', 'gp120 - C2', 'gp120 - C3', 'gp120 - C4', 'gp120 - C5'])]
    rpy_gp120['env region'] = ['gp120']*len(rpy_gp120)
    rpy_melt = rpy_melt.loc[rpy_melt['env region'].isin(['gp120 - V1', 'gp120 - V2', 'gp120 - V3', 'gp120 - V4', 'gp120 - V5',
                                                        'gp120 - C1', 'gp120 - C2', 'gp120 - C3', 'gp120 - C4', 'gp120 - C5'])]
    rpy_melt['env region'] = rpy_melt['env region'].apply(split_)
    rpy_melt = pd.concat([rpy_melt,rpy_gp120], axis=0)
    rpy_melt['protein region'] = rpy_melt['env region']

    isns.set_context('paper',fontfamily='Times New Roman')
    sns.set_context("paper", rc={"font.size":10,"axes.titlesize":10,"axes.labelsize":10})  
    g = sns.FacetGrid(rpy_melt, height=1.9, aspect=1.7, despine=False)
    g.map(sns.barplot, x="events/AA/year", y="protein region", hue='indel event', order=(['gp120', 'V1','V2','V3','V4','V5']),
        #['gp120 - C1','gp120 - V1','gp120 - V2','gp120 - C2','gp120 - V3','gp120 - C3','gp120 - V4','gp120 - C4','gp120 - V5','gp120 - C5']
        errorbar=('ci', 95), orient='h', linewidth=0.2, errwidth=0.5, 
                capsize=.1, palette={'insertion':'#00539C', 'deletion':lighten_color('#EEA47F',1.5), 'substitution':'red'},
                data=rpy_melt, width=0.8)

    for ax in g.axes.flat:
        sns.set_context("paper", rc={"font.size":10,"axes.titlesize":10,"axes.labelsize":10})
        ax.set_axisbelow(True)
        ax.grid(True, which='both', axis='x', zorder=-1, linestyle='dashed',
                linewidth=0.55)
        plt.setp(ax.get_yticklabels(),  ha='center', fontsize=8)
        # rotation=90,
        plt.setp(ax.get_xticklabels(), fontsize=8)
        ax.tick_params(axis='y', direction='out', pad=13)
        ax.tick_params(axis=u'both', which=u'both',length=0)
    # plt.legend(fontsize=8, bbox_to_anchor=(-0.5,3.7))
    g.savefig(os.path.join(bp,'env_protein_indel_rates_HIV_gp120_variable_event.png'), dpi=400)


def plot_hotspots(bp, path_to_events_all):
    events_df = pd.read_csv(path_to_events_all)
    only_position = events_df[['region','event','position w/o placeholders']]
    
    def split_(x):
        return x.split('-')[1].lstrip()

    rpy_gp120 = only_position.loc[only_position['region'].isin(['gp120 - V1', 'gp120 - V2', 'gp120 - V3', 'gp120 - V4', 'gp120 - V5',
                                                            'gp120 - C1', 'gp120 - C2', 'gp120 - C3', 'gp120 - C4', 'gp120 - C5'])]
    
    regs_list = ['gp120 - V1', 'gp120 - V2', 'gp120 - V3', 'gp120 - V4', 'gp120 - V5',
                                                            'gp120 - C1', 'gp120 - C2', 'gp120 - C3', 'gp120 - C4', 'gp120 - C5']
    
    rpy_gp120['region'] = ['gp120']*len(rpy_gp120)
    only_position = only_position.loc[only_position['region'].isin(regs_list)]
    only_position['region'] = only_position['region'].apply(split_)
    only_position = pd.concat([only_position,rpy_gp120], axis=0, ignore_index=True)
    isns.set_context('paper',fontfamily='Times New Roman')
    isns.set_context("paper", rc={"font.size":8,"axes.titlesize":8,"axes.labelsize":8})
    only_position_indel = only_position[only_position['event'].isin(['insertion','deletion'])]
    for reg in only_position_indel.region.unique(): 
        isns.set_context('paper',fontfamily='Times New Roman')
        isns.set_context("paper", rc={"font.size":8,"axes.titlesize":10,"axes.labelsize":10})
        g = sns.displot(data=only_position_indel[only_position_indel['region']==reg], x='position w/o placeholders', hue='event', multiple='stack', 
                        discrete=True, 
                        height=1.7, aspect=2,binwidth=2,
                        palette={'insertion':'#00539C', 'deletion':lighten_color('#EEA47F',1.5)})

        g.tight_layout()
        for ax in g.axes.flat:
            isns.set_context('paper',fontfamily='Times New Roman')
            sns.set_context("paper", rc={"font.size":8,"axes.titlesize":8,"axes.labelsize":8})
        
            ax.set_axisbelow(True)
            plt.setp(ax.get_yticklabels(), rotation=90, ha='center', fontsize=8)
            plt.setp(ax.get_xticklabels(), fontsize=8)
            ax.tick_params(axis='y', direction='out', pad=5)
            ax.tick_params(axis=u'both', which=u'both',length=0)
            ax.set_xlabel('env region = '+reg)
            ax.xaxis.set_label_coords(0.5, -0.2)
            ax.yaxis.set_label_coords(-0.1,0.5)
        g.savefig(os.path.join(bp,'indel_hotspots.png'), dpi=400)

