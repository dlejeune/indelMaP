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
from RateMatrix import WAG, blosum, JC69, K80
import os
import argparse

def ParsLeaf(leaf, alphabet, indel_aware=True):
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
    for i in range(len(leaf.sequence)):
        if alphabet == 'Protein':
    
            if leaf.sequence[i] == 'X':
                pars_sets.append(set(alphabet))
            
            elif leaf.sequence[i] == 'B':
                pars_sets.append(set(['D', 'N']))
            
            elif leaf.sequence[i] == 'Z':
                pars_sets.append(set(['E', 'Q']))
            
            elif leaf.sequence[i] == 'J':
                pars_sets.append(set(['I', 'L']))
                
            else:
                pars_sets.append(set(leaf.sequence[i]))
        
        if alphabet == 'DNA':
            if leaf.sequence[i] == 'X' or leaf.sequence[i] == 'N':
                pars_sets.append(set(alphabet))
                
            elif leaf.sequence[i] == 'V':
                pars_sets.append(set(['A', 'C', 'G']))
            
            elif leaf.sequence[i] == 'H':
                pars_sets.append(set(['A', 'C', 'T']))
                    
            elif leaf.sequence[i] == 'D':
                pars_sets.append(set(['A', 'G', 'T']))
                
            elif leaf.sequence[i] == 'B':
                pars_sets.append(set(['C', 'G', 'T']))
            
            elif leaf.sequence[i] == 'M':
                pars_sets.append(set(['A', 'C']))
            
            elif leaf.sequence[i] == 'R':
                pars_sets.append(set(['A', 'G']))
            
            elif leaf.sequence[i] == 'W':
                pars_sets.append(set(['A', 'T']))
            
            elif leaf.sequence[i] == 'S':
                pars_sets.append(set(['C', 'G']))
                
            elif leaf.sequence[i] == 'Y':
                pars_sets.append(set(['C', 'T']))
            
            elif leaf.sequence[i] == 'K':
                pars_sets.append(set(['G', 'T']))
            
            else:
                pars_sets.append(set(leaf.sequence[i]))
    
    pars_scores = [0]*len(leaf.sequence)
    ins_gaps = [False] *len(leaf.sequence)

    if indel_aware:
        ins_flags = [False]*len(leaf.sequence)
        leaf.add_features(insertion_flags = ins_flags)
    
    leaf.add_features(parsimony_sets = pars_sets)
    leaf.add_features(parsimony_scores = pars_scores)
    leaf.add_features(insertion_gaps = ins_gaps)

  
def ParsInternal(tree, i, C_all, gi_f, ge_f, indel_aware, branch_length):
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
                                          branch_length)
    C_right, gi_right, ge_right = determineC(C_all, right_dist, gi_f, ge_f, 
                                             branch_length)
    
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
        if i > 0:
            
                
            #if we expand an insertion site
            if indel_aware and tree.insertion_flags[i]:
                tmp = i-1
                while (tree.insertion_gaps[tmp] and (not tree.insertion_flags[tmp]) and tmp > 0):
                    tmp -= 1
                #print(left_set[tmp])
                if left_set[tmp] == set('-') and tree.insertion_flags[tmp]:
                    tree.parsimony_scores[i] = left_score + right_score + ge_left
                else:
                    tree.parsimony_scores[i] = left_score + right_score + gi_left
            
            elif indel_aware and not tree.insertion_flags[i]:
                tmp = i-1
                while (tree.insertion_gaps[tmp] or tree.insertion_flags[tmp]) and tmp > 0:
                    tmp -= 1
                
                if left_set[tmp] == set('-') and (not tree.insertion_gaps[tmp]) and (not tree.insertion_flags[tmp]):
                    tree.parsimony_scores[i] = left_score + right_score + ge_left
                else:
                    tree.parsimony_scores[i] = left_score + right_score + gi_left
            else:
                tree.parsimony_scores[i] = left_score + right_score + gi_left
        else:
            tree.parsimony_scores[i] = left_score + right_score + gi_left
                
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
        if i > 0:
                
            #if we expand an insertion site
            if indel_aware and tree.insertion_flags[i]:
                tmp = i-1
                while tree.insertion_gaps[tmp] and (not tree.insertion_flags[tmp]) and tmp > 0:
                    tmp -= 1
                
                if right_set[tmp] == set('-') and tree.insertion_flags[tmp]:
                    tree.parsimony_scores[i] = left_score + right_score + ge_right
                else:
                    tree.parsimony_scores[i] = left_score + right_score + gi_right
                    
            elif indel_aware and not tree.insertion_flags[i]:
                tmp = i-1
            
                while (tree.insertion_gaps[tmp] or tree.insertion_flags[tmp]) and tmp > 0:
                    tmp -= 1
                
                if right_set[tmp] == set('-') and (not tree.insertion_gaps[tmp]) and (not tree.insertion_flags[tmp]):
                    tree.parsimony_scores[i] = left_score + right_score + ge_right
                else:
                    tree.parsimony_scores[i] = left_score + right_score + gi_right
            else:
                tree.parsimony_scores[i] = left_score + right_score + gi_right
        else:
            tree.parsimony_scores[i] = left_score + right_score + gi_right
        
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
                score_l = C_left[left_character][int_character]
                
                for right_character in right_set[i]:
                    score_r = C_right[right_character][int_character]
                    
                    score = score_l + score_r
                    
                    if score < min_score:
                        min_score = score  
                        
        tree.parsimony_scores[i] = left_score + right_score + min_score
    
    elif not left_set[i].intersection(right_set[i]):
        tree.parsimony_sets[i] = left_set[i].union(right_set[i])
        
        min_score = np.inf
        for left_character in left_set[i]:
            for int_character in tree.parsimony_sets[i]:
                score_l = C_left[left_character][int_character]
                
                for right_character in right_set[i]:
                    score_r = C_right[right_character][int_character]
                    
                    score = score_l + score_r
                    
                    if score < min_score:
                        min_score = score        
                
        tree.parsimony_scores[i] = left_score + right_score + min_score
    
        
def ParsAncestral(tree, indel_aware):
    seq = ''
    events = ''
    
    for i in range(len(tree.parsimony_sets)):
        if tree.is_root():
            character = sample(tree.parsimony_sets[i],1)[0]
            
            if character == '-' and (indel_aware and (tree.insertion_gaps[i] or tree.insertion_flags[i])):
                character = '*'
        
        else:
            if tree.parsimony_sets[i] == set('-'):
                if indel_aware and (tree.insertion_gaps[i] or tree.insertion_flags[i]):
                    character = '*'
                else:
                    character = '-'
            
            else:
                if len(tree.parsimony_sets[i]) > 1 and tree.up.sequence[i] in tree.parsimony_sets[i]:
                        character = tree.up.sequence[i]
                else:
                    if indel_aware and tree.up.insertion_flags[i]:
                        character = sample(tree.parsimony_sets[i],1)[0].lower()
                    else:
                        if tree.up.evolutionary_events[i].islower():
                            character = sample(tree.parsimony_sets[i],1)[0].lower()
                        else:
                            character = sample(tree.parsimony_sets[i],1)[0]
            
        events += character
        if character == '*':
            seq += '-'
        else:
            seq += character.upper()
        
    tree.add_features(evolutionary_events = events)
    tree.add_features(sequence = seq)
                    
                
def ParsASR(tree_file, msa_file, alphabet, out_file=os.path.abspath(os.getcwd())+'/msa', Q=None, gi_f=2.5, ge_f=0.5, indel_aware=True,
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
    
    
    C_all = calculateC(alphabet, Q, branch_length)
    
    # find insertion points and mark them with an insertion flag set to True
    if indel_aware:
        for i in range(length_MSA):
            leaf_res = []
            no_leaves = 0
            for leaf in tree.iter_leaves():
                no_leaves += 1
                if i == 0:
                    #sets and scores for leaves
                    ParsLeaf(leaf, alphabet, indel_aware)
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
            ParsLeaf(leaf, alphabet, indel_aware)

                    
    #find internal sets and scores
    for node in tree.traverse('postorder'):
        if not node.is_leaf():
            for i in range(length_MSA):
                    
                ParsInternal(node, i, C_all, gi_f, ge_f, 
                                     indel_aware, branch_length)
    # sum the parsimony scores at the root over the whole sequence
            # print(node.insertion_gaps)
            # print(node.insertion_flags)
            # print(node.parsimony_scores)
    tree_score = sum(tree.parsimony_scores)
    
    if ancestor_reconstruction:
        with open(out_file + '_internal_evolutionary_events.fasta', 'w') as f1, open(out_file + '_internal_ancestral_reconstruction.fasta', 'w') as f3, open(out_file + '_leaves_evolutionary_events.fasta', 'w') as f2:
            
            no_internal = len(tree) + 1
            for node in tree.traverse('preorder'):
                
                if node.is_root():
                    node.name = 'ROOT'
                else:
                    node.name = 'N' + str(no_internal)
                    no_internal += 1
                    
                ParsAncestral(node, indel_aware)
                
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
    elif 'K80' in Q[0]:
        q = K80(float(Q[0].split('0')[1]),float(Q[1].split('0')[1]))
    elif Q[0] == 'JC69':
        q = JC69
    elif Q[0] == 'None':
        q = None
        
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
    parsimony_score = ParsASR(tree_file, sequence_file, alphabet, output_file,
                          q, gi_f, ge_f, indel_aware,
                          branch_length, ancestral_reconstruction)

    print('\nThe files including the ancestral sequence reconstruction and evoultionary events were saved to:\n'+ output_file +  '_internal_evolutionary_events.fasta\n'+  output_file + '_internal_ancestral_reconstruction.fasta\n'  + output_file + '_leaves_evolutionary_events.fasta\n' + '\nIndel-aware parsimony score: '+str(parsimony_score)+'\n')



if __name__ == '__main__':
    main()
