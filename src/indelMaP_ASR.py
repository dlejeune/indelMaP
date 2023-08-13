#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 29 22:23:09 2021

@author: claraiglhaut
"""

import numpy as np
from ete3 import PhyloNode
from random import sample
from calculateC import calculateC
from HelperFunctions import determineC
from RateMatrix import WAG, blosum, HIVb, JC69, K80, GTR
import os
import argparse

def indelMaP_Leaf(leaf, alphabet, indel_aware=True):
    '''
    Initializes the parsimony sets and scores at the leaf nodes.

    Parameters
    ----------
    leaf : PhlyoNode or PhyloTree
        tree leaves with the aligned sequences.

    Returns
    -------
    None.

    '''
    pars_sets = []
    characters_protein = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P',
                          'S', 'T', 'W', 'Y', 'V']
    characters_dna = ['T', 'C', 'A', 'G']
    for i in range(len(leaf.sequence)):
        in_character = leaf.sequence[i].upper()
        if alphabet == 'Protein':
            if in_character == 'X':
                pars_sets.append(set(characters_protein))
            elif in_character == 'B':
                pars_sets.append(set(['D', 'N']))
            elif in_character == 'Z':
                pars_sets.append(set(['E', 'Q']))
            elif in_character == 'J':
                pars_sets.append(set(['I', 'L']))
            else:
                pars_sets.append(set(in_character))
        if alphabet == 'DNA':
            if in_character == 'X' or in_character == 'N':
                pars_sets.append(set(characters_dna))
            elif in_character == 'V':
                pars_sets.append(set(['A', 'C', 'G']))
            elif in_character == 'H':
                pars_sets.append(set(['A', 'C', 'T']))
            elif in_character == 'D':
                pars_sets.append(set(['A', 'G', 'T']))
            elif in_character == 'B':
                pars_sets.append(set(['C', 'G', 'T']))
            elif in_character == 'M':
                pars_sets.append(set(['A', 'C']))
            elif in_character == 'R':
                pars_sets.append(set(['A', 'G']))
            elif in_character == 'W':
                pars_sets.append(set(['A', 'T']))
            elif in_character == 'S':
                pars_sets.append(set(['C', 'G']))
            elif in_character == 'Y':
                pars_sets.append(set(['C', 'T']))
            elif in_character == 'K':
                pars_sets.append(set(['G', 'T']))
            else:
                pars_sets.append(set(in_character))    
    pars_scores = [0]*len(leaf.sequence)
    ins_gaps = [False] *len(leaf.sequence)
    if indel_aware:
        ins_flags = [False]*len(leaf.sequence)
        leaf.add_features(insertion_flags = ins_flags)
    leaf.add_features(parsimony_sets = pars_sets)
    leaf.add_features(parsimony_scores = pars_scores)
    leaf.add_features(insertion_gaps = ins_gaps)

  
def indelMaP_Internal(tree, i, C_all, gi_f, ge_f, indel_aware, branch_length, bl_percentiles):
    '''
    Creates the parsimony sets and scores for the internal nodes. 

    Parameters
    ----------
    tree : PhyloNode or PhlyoTree
        Internal nodes a tree structue.
    i : int
        Index for the sequence.
    cost_matrix: double dictionary
        Gives the scores for matching, mismatching characters
        Give None to get the unweighted version.
    go : float
        Gap opening penalty.
    ge : float
        Gap extension penalty.
    indel_aware: boolean
        Default is set to True. If set to false, insertions are penalized
        multiple times.

    Returns
    -------
    None.

    '''
    left_set = tree.children[0].parsimony_sets
    left_score = tree.children[0].parsimony_scores[i]
    right_set = tree.children[1].parsimony_sets
    right_score  = tree.children[1].parsimony_scores[i]
    left_dist = tree.children[0].dist
    right_dist = tree.children[1].dist
    C_left, gi_left, ge_left = determineC(C_all, left_dist, gi_f, ge_f, 
                                          branch_length, bl_percentiles)
    C_right, gi_right, ge_right = determineC(C_all, right_dist, gi_f, ge_f, 
                                             branch_length, bl_percentiles)
    length_MSA = len(left_set)
    if i == 0:
        pars_sets = [set()] * length_MSA
        pars_scores = [0] * length_MSA
        ins_gaps = [False] * length_MSA
        tree.add_features(parsimony_sets = pars_sets)
        tree.add_features(parsimony_scores = pars_scores)
        tree.add_features(insertion_gaps = ins_gaps)
        if indel_aware and not hasattr(tree, 'insertion_flags'):
            ins_flags = [False] * length_MSA
            tree.add_features(insertion_flags = ins_flags)
    #matching two gap characters
    if left_set[i] == set('-') and right_set[i] == set('-'):
        tree.parsimony_sets[i] = set('-')
        tree.parsimony_scores[i] = left_score + right_score
        tree.insertion_gaps[i] = True 
    elif (left_set[i] == set('-') and right_set[i] != set('-')):
        min_score_r = np.inf
        for right_character in right_set[i]:
            for int_character in right_set[i]:
                score_r = C_right[int_character][right_character]
                if score_r < min_score_r:
                    min_score_r = score_r
        if i > 0:
            #if we expand an insertion site
            if indel_aware and tree.insertion_flags[i]:
                tmp = i-1
                while (tree.insertion_gaps[tmp] and (not tree.insertion_flags[tmp]) and tmp > 0):
                    tmp -= 1
                if left_set[tmp] == set('-') and tree.insertion_flags[tmp]:
                    tree.parsimony_scores[i] = left_score + right_score + ge_left + min_score_r
                else:
                    tree.parsimony_scores[i] = left_score + right_score + gi_left + min_score_r
            elif indel_aware and not tree.insertion_flags[i]:
                tmp = i-1
                while (tree.insertion_gaps[tmp] or tree.insertion_flags[tmp]) and tmp > 0:
                    tmp -= 1
                if left_set[tmp] == set('-') and (not tree.insertion_gaps[tmp]) and (not tree.insertion_flags[tmp]):
                    tree.parsimony_scores[i] = left_score + right_score + ge_left + min_score_r
                else:
                    tree.parsimony_scores[i] = left_score + right_score + gi_left + min_score_r
            else:
                tree.parsimony_scores[i] = left_score + right_score + gi_left + min_score_r
        else:
            tree.parsimony_scores[i] = left_score + right_score + gi_left + min_score_r
        if indel_aware and tree.insertion_flags[i]:
            tree.parsimony_sets[i] = set('-')
            for node in tree.iter_ancestors():
                if hasattr(node, 'insertion_gaps'):
                    node.insertion_gaps[i] = True 
                else:
                    ins_gaps = [False] * length_MSA
                    node.add_features(insertion_gaps = ins_gaps)
                    node.insertion_gaps[i] = True 
        else:
            tree.parsimony_sets[i] = right_set[i]
    elif (left_set[i] != set('-') and right_set[i] == set('-')):
        min_score_l = np.inf
        for left_character in left_set[i]:
            for int_character in left_set[i]:
                score_l = C_left[int_character][left_character]
                if score_l < min_score_l:
                    min_score_l = score_l
        if i > 0:
            #if we expand an insertion site
            if indel_aware and tree.insertion_flags[i]:
                tmp = i-1
                while tree.insertion_gaps[tmp] and (not tree.insertion_flags[tmp]) and tmp > 0:
                    tmp -= 1
                if right_set[tmp] == set('-') and tree.insertion_flags[tmp]:
                    tree.parsimony_scores[i] = left_score + right_score + ge_right + min_score_l
                else:
                    tree.parsimony_scores[i] = left_score + right_score + gi_right + min_score_l
            # if we expand a deletion site
            elif indel_aware and not tree.insertion_flags[i]:
                tmp = i-1
                while (tree.insertion_gaps[tmp] or tree.insertion_flags[tmp]) and tmp > 0:
                    tmp -= 1
                if right_set[tmp] == set('-') and (not tree.insertion_gaps[tmp]) and (not tree.insertion_flags[tmp]):
                    tree.parsimony_scores[i] = left_score + right_score + ge_right + min_score_l
                else:
                    tree.parsimony_scores[i] = left_score + right_score + gi_right + min_score_l
            else:
                tree.parsimony_scores[i] = left_score + right_score + gi_right + min_score_l
        else:
            tree.parsimony_scores[i] = left_score + right_score + gi_right + min_score_l
        if indel_aware and tree.insertion_flags[i]:
            tree.parsimony_sets[i] = set('-')
            for node in tree.iter_ancestors():
                if hasattr(node, 'insertion_gaps'):
                    node.insertion_gaps[i] = True 
                else:
                    ins_gaps = [False] * length_MSA
                    node.add_features(insertion_gaps = ins_gaps)
                    node.insertion_gaps[i] = True 
        else:
            tree.parsimony_sets[i] = left_set[i]
    elif left_set[i].intersection(right_set[i]):
        tree.parsimony_sets[i] = left_set[i].intersection(right_set[i])
        min_score = np.inf
        for left_character in left_set[i]:
            for int_character in tree.parsimony_sets[i]:
                score_l = C_left[int_character][left_character]
                for right_character in right_set[i]:
                    score_r = C_right[int_character][right_character]
                    score = score_l + score_r
                    if score < min_score:
                        min_score = score  
        tree.parsimony_scores[i] = left_score + right_score + min_score
    elif not left_set[i].intersection(right_set[i]):
        tree.parsimony_sets[i] = left_set[i].union(right_set[i])
        min_score = np.inf
        for left_character in left_set[i]:
            for int_character in tree.parsimony_sets[i]:
                score_l = C_left[int_character][left_character]
                for right_character in right_set[i]:
                    score_r = C_right[int_character][right_character]
                    score = score_l + score_r
                    if score < min_score:
                        min_score = score        
        tree.parsimony_scores[i] = left_score + right_score + min_score


def indelMaP_Ancestral(tree, C_all, gi_f, ge_f, indel_aware, branch_length, bl_percentiles, alphabet):
    seq = ''
    events = ''
    for i in range(len(tree.parsimony_sets)):
        if tree.is_root():
            character = best_reconstruction(tree, C_all, gi_f, ge_f, indel_aware, branch_length, bl_percentiles, i, alphabet)     
        else:
            if tree.parsimony_sets[i] == set('-'):
                if indel_aware and (tree.insertion_gaps[i] or tree.insertion_flags[i]) and tree.up.evolutionary_events[i]=='*':
                    character = '*'
                else:
                    character = '-'            
            else:
                if len(tree.parsimony_sets[i]) > 1 and tree.up.sequence[i] in tree.parsimony_sets[i]:
                        if indel_aware and tree.up.insertion_flags[i]:
                            character = tree.up.sequence[i].lower()
                        else:
                            character = tree.up.sequence[i]
                else:
                    if indel_aware and tree.up.insertion_flags[i]:
                        character = best_reconstruction(tree, C_all, gi_f, ge_f, indel_aware, branch_length, bl_percentiles, i, alphabet).lower()
                    else:
                        character = best_reconstruction(tree, C_all, gi_f, ge_f, indel_aware, branch_length, bl_percentiles, i, alphabet)       
        events += character
        if character == '*':
            seq += '-'
        else:
            seq += character.upper()
    tree.add_features(evolutionary_events = events)
    tree.add_features(sequence = seq)

def best_reconstruction(tree, C_all, gi_f, ge_f, indel_aware, branch_length, bl_percentiles, i, alphabet):
    if len(tree.parsimony_sets[i]) == 1:
        character = list(tree.parsimony_sets[i])[0]
        if character == '-' and (indel_aware and (tree.insertion_gaps[i] or tree.insertion_flags[i])):
                character = '*'
    elif tree.is_leaf():
        # Dealing with ambigous characters
        if alphabet == 'Protein':
            characters_protein = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P',
                          'S', 'T', 'W', 'Y', 'V']
            if tree.parsimony_sets[i] == set(characters_protein):
                character = 'X'
            elif tree.parsimony_sets[i] == set(['D', 'N']):
                character = 'B'
            elif tree.parsimony_sets[i] == set(['E', 'Q']):
                character = 'Z'
            elif tree.parsimony_sets[i] == set(['I', 'L']):
                character = 'J'
        elif alphabet == 'DNA':
            characters_dna = ['T', 'C', 'A', 'G']
            if tree.parsimony_sets[i] == set(characters_dna):
                character = 'X'
            elif tree.parsimony_sets[i] == set(['A', 'C', 'G']):
                character = 'V'
            elif tree.parsimony_sets[i]== set(['A', 'C', 'T']):
                character = 'H'
            elif tree.parsimony_sets[i] == set(['A', 'G', 'T']):
                charachter = 'D'
            elif tree.parsimony_sets[i] == set(['C', 'G', 'T']):
                character = 'B'
            elif tree.parsimony_sets[i] == set(['A', 'C']):
                character = 'M'
            elif tree.parsimony_sets[i] == set(['A', 'G']):
                character = 'R'
            elif tree.parsimony_sets[i] == set(['A', 'T']):
                character = 'W'
            elif tree.parsimony_sets[i] == set(['C', 'G']):
                character = 'S'
            elif tree.parsimony_sets[i] == set(['C', 'T']):
                character = 'Y'
            elif tree.parsimony_sets[i] == set(['G', 'T']):
                character = 'K'
        
    else:
        left_dist = tree.children[0].dist
        right_dist = tree.children[1].dist
        C_left, gi_left, ge_left = determineC(C_all, left_dist, gi_f, ge_f, 
                                                    branch_length, bl_percentiles)
        C_right, gi_right, ge_right = determineC(C_all, right_dist, gi_f, ge_f, 
                                                        branch_length, bl_percentiles)
        left_set = tree.children[0].parsimony_sets
        right_set = tree.children[1].parsimony_sets
        min_char = []
        if left_set[i] == set('-'):
            min_score_r = np.inf
            for right_character in right_set[i]:
                for int_character in tree.parsimony_sets[i]:
                    score_r = C_right[int_character][right_character]
                    if score_r == min_score_r:
                        min_char += [int_character]
                    elif score_r < min_score_r:
                        min_score_r = score_r
                        min_char = [int_character]
        elif right_set[i] == set('-'):
            min_score_l= np.inf
            for left_character in left_set[i]:
                for int_character in tree.parsimony_sets[i]:
                    score_l = C_left[int_character][left_character]
                    if score_l == min_score_l:
                        min_char += [int_character]
                    elif score_l < min_score_l:
                        min_score_l = score_l
                        min_char = [int_character]
        else:
            min_score = np.inf
            for left_character in left_set[i]:
                for int_character in tree.parsimony_sets[i]:
                    score_l = C_left[int_character][left_character]
                    for right_character in right_set[i]:
                        score_r = C_right[int_character][right_character]
                        score = score_l + score_r
                        if score == min_score:
                            min_char += [int_character]
                        elif score < min_score:
                            min_score = score
                            min_char = [int_character]
        if len(min_char)>1:
            character = sample(sorted(min_char),1)[0]   
        else:
            character = min_char[0]
    return character
                    
                
def indelMaP_ASR(tree_file, msa_file, alphabet, out_file=os.path.abspath(os.getcwd())+'/msa', Q=None, gi_f=2.5, ge_f=0.5, indel_aware=True,
              branch_length = True, ancestor_reconstruction=True):
    '''
    Calculates the parsimony score for the whole tree while accounting for 
    insertions and deletions.

    Parameters
    ----------
    tree_file : str
        Path to file containing a guide tree in newick format.
    msa_file : str
        Path to alignemnt file in fasta format.
    alphabet : str
        Choose between 'Protein' or 'DNA'.
    output_file : str
        Path and filename for the outputfiles without suffix;
        if left empty the files will save to the current working directory. Files will only be printed if ancestral_reconstruction is True.
    Q : numpy.ndarray
        Choose the substitution model; - Protein: WAG, blosum - DNA: JC69, K80(alpha,beta);
        if no substitution model is specified each substitution is associated with the same cost.
        You can specify your own model by giving a symmetric transition rate matrix with an average substitution rate of 1,
        as two dimensional np.array enclosed in a tuple (np.array,), substitution rates are given as float
    gi_f : float, default = 2.5
        Gives factor for gap opening panelty.
    ge_f : float, default = 0.5
        Gives factor for gap extension penalty.
    indel_aware : boolean, default = True
        If set to False the algorithm does not distinguish between insertion and deletion events.
    branch_length : boolean, default = True
        If set to False the cost matrix is calculated based on the transition probability matrix for branch length 0.5 for each branch;
        if True the cost matrix is calculated based on one of the four distances [0.1,0.3,0.5,0.7], the distance closest to the branch length is chosen.
    ancestor_reconstruction : boolean, default = True
        If set to true ancestral sites and evolutionary events are reconstructed.

    Returns
    -------
    tree_score : numpy.float64
        parsimony score for the whole tree.

    ''' 
    tree = PhyloNode(newick = tree_file, alignment = msa_file, format=1, quoted_node_names=True)
    length_MSA = len(tree.get_leaves()[0].sequence)
    C_all, bl_percentiles = calculateC(alphabet, Q, branch_length, tree_file, tree)
    # find insertion points and mark them with an insertion flag set to True
    if indel_aware:
        for i in range(length_MSA):
            leaf_res = []
            no_leaves = 0
            for leaf in tree.iter_leaves():
                no_leaves += 1
                if i == 0:
                    #sets and scores for leaves
                    indelMaP_Leaf(leaf, alphabet, indel_aware)
                if leaf.sequence[i] != '-':
                    leaf_res.append(leaf)
            if len(leaf_res) == 1:
                ancestor = leaf_res[0]
            else:
                ancestor = tree.get_common_ancestor(leaf_res)
            if not ancestor.is_root():
                if hasattr(ancestor.up, 'insertion_flags'):
                    ancestor.up.insertion_flags[i] = True
                else:
                    ins_flags = [False] * length_MSA
                    ancestor.up.add_features(insertion_flags = ins_flags)
                    ancestor.up.insertion_flags[i] = True
            for node in tree.iter_leaves():
                if not node in ancestor.iter_leaves():
                    node.insertion_gaps[i] = True
    else:
        for leaf in tree.iter_leaves():
            indelMaP_Leaf(leaf, alphabet, indel_aware)
    #find internal sets and scores
    for node in tree.traverse('postorder'):
        if not node.is_leaf():
            for i in range(length_MSA):
                indelMaP_Internal(node, i, C_all, gi_f, ge_f, 
                                     indel_aware, branch_length, bl_percentiles)
    # sum the parsimony scores at the root over the whole sequence
    tree_score = sum(tree.parsimony_scores)
    if ancestor_reconstruction:
        with open(out_file + '_internal_evolutionary_events.fas', 'w') as f1, open(out_file + '_internal_ancestral_reconstruction.fas', 'w') as f3, open(out_file + '_leaves_evolutionary_events.fas', 'w') as f2:
            no_internal = len(tree) + 1
            for node in tree.traverse('preorder'):
                if node.name == '':
                    if node.is_root():
                        node.name = 'ROOT'
                    elif not node.is_leaf():
                        node.name = 'N' + str(no_internal)
                        no_internal += 1
                        
                indelMaP_Ancestral(node, C_all, gi_f, ge_f, indel_aware, branch_length, bl_percentiles)
                
                if not node.is_leaf():
                    print('>'+ node.name+'\n'+node.evolutionary_events, file=f1)
                    print('>'+ node.name+ '\n'+ node.sequence, file=f3)
                else:
                    print('>'+ node.name+'\n'+node.evolutionary_events, file=f2)
        f1.close()
        f2.close()
        f3.close()
        tree.write(format=1, outfile=out_file+'_tree.nwk')
    return tree_score


def get_arguments_from_CLI():
    parser = argparse.ArgumentParser(prog = 'Indel-aware parsimony ancestral reconstruction', description='The program reconstructs ancestral sequences for protein or nucleotide sequences along a given guide tree under indel-aware parsimony.')
    parser.add_argument('--msa_file', '-m', help='Path to alignemnt file in fasta format.', required = True)
    parser.add_argument('--tree_file', '-t',  help='Path to file containing a guide tree in newick format.', required = True)
    parser.add_argument('--output_file', '-o', default=os.path.abspath(os.getcwd())+'/msa', help = 'Path and filename for the outputfiles without suffix; if left empty the files will save to the current working directory. Files will only be printed if ancestral_reconstruction is True.', required = False)
    parser.add_argument('--alphabet', '-a',  help='Specify sequence alphabet, choose between DNA or Protein', type=str, required = True)
    parser.add_argument('--RateMatrix','-q', default = None, type=str, nargs='+', help = 'Choose the substitution model; - Protein: WAG, blosum - DNA: JC69, K80{alpha,beta}; if no substitution model is specified the default for DNA sequences is K80 with transition to transversion ratio of 2 and for Protein sequences WAG; if each substitution should be associated with the same cost you can type None. You can specify your own model by giving a symmetric transition rate matrix with an average substitution rate of 1. Columns have to be separated by a comma and rows with a colon, and transition rates are given as float. Example: -q=-1,0.3333333333333333,0.3333333333333333,0.3333333333333333:0.3333333333333333,-1,0.3333333333333333,0.3333333333333333:0.3333333333333333,0.3333333333333333,-1,0.3333333333333333:0.3333333333333333,0.3333333333333333,0.3333333333333333,-1')
    parser.add_argument('--gap_opening_factor','-go', default = 2.5, help = 'The gap opening cost is given by the gap_opening_factor*average_substitution_cost; default = 2.5', type=float, required = False)
    parser.add_argument('--gap_extension_factor', '-ge', default = 0.5, help = 'The gap extension cost is given by the gap_extension_factor*average_substitution_cost; default = 0.5', type=float, required = False)
    parser.add_argument('--indel-aware', '-ia', default=True, help = 'If set to False the algorithm does not distinguish between insertion and deletion events; default = True', type=bool, required = False)
    parser.add_argument('--branch_length', '-bl', default=True, help = 'If set to False the cost matrix is calculated based on the transition probability matrix for branch length 0.5 for each branch; if True the cost matrix is calculated based on one of the four distances [0.1,0.3,0.5,0.7], the distance closest to the branch length is chosen.', type=bool, required=False)
    parser.add_argument('--ancestral_reconstruction', '-asr', default=True, help = 'If set to true ancestral sites and evolutionary events are reconstructed.', type=bool, required=False)
    args = parser.parse_args()
    return args

def main():
    args = get_arguments_from_CLI()
    tree_file = args.tree_file
    sequence_file = args.msa_file
    output_file = args.output_file
    alphabet = args.alphabet
    Q = args.RateMatrix
    if Q == None:
        if alphabet == 'Protein':
            q = WAG
        elif alphabet == 'DNA':
            q = K80(2,1)
    elif Q[0] == 'WAG':
        q = WAG
    elif Q[0] == 'blosum':
        q = blosum
    elif Q[0] == 'HIVb':
        q = HIVb
    elif 'K80' in Q[0]:
        q = K80(float(Q[0].split('0')[1]),float(Q[1].split('0')[1]))
    elif Q[0] == 'JC69':
        q = JC69
    elif Q[0] == 'None':
        q = None
    elif 'GTR' in Q[0]:
        q = GTR(float(Q[0].split('GTR')[1]),float(Q[1].split('R')[1]),
                float(Q[2].split('R')[1]),float(Q[3].split('R')[1]),
                float(Q[4].split('R')[1]),float(Q[5].split('R')[1]),
                float(Q[6].split('R')[1]),float(Q[7].split('R')[1]),
                float(Q[8].split('R')[1]),float(Q[9].split('R')[1]))    
    else:
        print('User defined rate matrix')
        Q = Q[0].split(':')
        q=np.empty((len(Q),len(Q)))
        for i in range(len(Q)):
            row = Q[i].split(',')
            for j in range(len(row)):
                q[i][j] = float(row[j])
        q = (q,1)
    gi_f = args.gap_opening_factor
    ge_f = args.gap_extension_factor
    indel_aware = args.indel_aware
    branch_length = args.branch_length
    ancestral_reconstruction = args.ancestral_reconstruction
    parsimony_score = indelMaP_ASR(tree_file, sequence_file, alphabet, output_file,
                          q, gi_f, ge_f, indel_aware,
                          branch_length, ancestral_reconstruction)
    print('\nThe files including the ancestral sequence reconstruction and evoultionary events were saved to:\n'+ output_file +  '_internal_evolutionary_events.fas\n'+  output_file + '_internal_ancestral_reconstruction.fas\n'  + output_file + '_leaves_evolutionary_events.fas\n' + '\nIndel-aware parsimony score: '+str(parsimony_score)+'\n')

if __name__ == '__main__':
    main()
