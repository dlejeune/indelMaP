#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 29 22:23:09 2021

@author: claraiglhaut
"""

import numpy as np
from ete3 import PhyloNode
from random import sample
from src.calculateC import calculateC
from src.HelperFunctions import determineC

def ParsLeaf(leaf, alphabet, phylogeny_aware=True):
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

    if phylogeny_aware:
        ins_flags = [False]*len(leaf.sequence)
        leaf.add_features(insertion_flags = ins_flags)
    
    leaf.add_features(parsimony_sets = pars_sets)
    leaf.add_features(parsimony_scores = pars_scores)
    leaf.add_features(insertion_gaps = ins_gaps)

  
def ParsInternal(tree, i, C_all, gi_f, ge_f, phylogeny_aware, branch_length):
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
    phylogeny_aware: boolean
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
        
        if phylogeny_aware and not hasattr(tree, 'insertion_flags'):
            ins_flags = [False] * length_MSA
            tree.add_features(insertion_flags = ins_flags)
    
    #matching two gap characters
    if left_set[i] == set('-') and right_set[i] == set('-'):
        tree.parsimony_sets[i] = set('-')
        tree.parsimony_scores[i] = left_score + right_score
        tree.insertion_gaps[i] = True
           
    
    elif (left_set[i] == set('-') and right_set[i] != set('-')):
        if i > 0:
            tmp = i-1
            while tree.insertion_gaps[tmp] == True and tmp > 0:
                tmp -= 1
                
            #if we expand an insertion site
            if left_set[tmp] == set('-') and not tree.insertion_gaps[tmp] == True:
                tree.parsimony_scores[i] = left_score + right_score + ge_left
            else:
                tree.parsimony_scores[i] = left_score + right_score + gi_left
        else:
            tree.parsimony_scores[i] = left_score + right_score + gi_left
                
        if phylogeny_aware and tree.insertion_flags[i]:
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
            tmp = i-1
            while tree.insertion_gaps[tmp] == True and tmp > 0:
                tmp -= 1
                
            #if we expand an insertion site
            if right_set[tmp] == set('-') and not tree.insertion_gaps[tmp] == True:
                tree.parsimony_scores[i] = left_score + right_score + ge_right
            else:
                tree.parsimony_scores[i] = left_score + right_score + gi_right
        else:
            tree.parsimony_scores[i] = left_score + right_score + gi_right
        
        if phylogeny_aware and tree.insertion_flags[i]:
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
        
def ParsAncestral(tree):
    seq = ''
    events = ''
    
    for i in range(len(tree.parsimony_sets)):
        if tree.is_root():
            character = sample(tree.parsimony_sets[i],1)[0]
            
            if character == '-' and (tree.insertion_gaps[i] or tree.insertion_flags[i]):
                character = '*'
        
        else:
            if tree.parsimony_sets[i] == set('-'):
                if (tree.insertion_gaps[i] or tree.insertion_flags[i]):
                    character = '*'
                else:
                    character = '-'
            
            else:
                if len(tree.parsimony_sets[i]) > 1 and tree.up.sequence[i] in tree.parsimony_sets[i]:
                        character = tree.up.sequence[i]
                else:
                    if tree.up.insertion_flags[i]:
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
                    
                
def ParsScore(tree_file, msa_file, out_file, alphabet, Q=None, gi_f=1, ge_f=1, phylogeny_aware=True,
              branch_length = True, ancestor_reconstruction=True):
    '''
    Calculates the parsimony score for the whole tree while accounting for 
    insertions and deletions.

    Parameters
    ----------
    tree_file : str
        path to newick tree file.
    msa_file : str
        path to alignemnt file in fasta format.
    cost_matrix: double dictionary
        Gives the scores for matching, mismatching characters
        Give None to get the unweighted version
    go : float
        Gap opening panelty
    ge : float
        Gap extension penalty
    phylogeny_aware: boolean
        Default is set to True. If set to false, insertions are penalized 
        multiple times
    branch_length : boolean, default = True
        If set to true the branch length is used to calculate a branch length 
        specific cost matrix, if set to False a standard distance of 0.62 is used
    ancestor_reconstruction : boolean, default = True
        If set to true ancestral sites and evolutionary events are reconstructed.

    Returns
    -------
    tree_score : int
        parsimony score for the whole tree.

    '''
    
    tree = PhyloNode(newick = tree_file, alignment = msa_file, format=1, quoted_node_names=True)
    length_MSA = len(tree.get_leaves()[0].sequence)
    
    
    C_all = calculateC(alphabet, Q, branch_length)
    
    # find insertion points and mark them with an insertion flag set to True
    if phylogeny_aware:
        for i in range(length_MSA):
            leaf_res = []
            no_leaves = 0
            for leaf in tree.iter_leaves():
                no_leaves += 1
                if i == 0:
                    #sets and scores for leaves
                    ParsLeaf(leaf, alphabet, phylogeny_aware)
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
            ParsLeaf(leaf, alphabet, phylogeny_aware)

                    
    #find internal sets and scores
    for node in tree.traverse('postorder'):
        if not node.is_leaf():
            for i in range(length_MSA):
                    
                ParsInternal(node, i, C_all, gi_f, ge_f, 
                                     phylogeny_aware, branch_length)

        
    # sum the parsimony scores at the root over the whole sequence
    tree_score = sum(tree.parsimony_scores)
    
    if ancestor_reconstruction:
        with (open(out_file + '_internal_evolutionary_events.fasta', 'w') as f1, 
              open(out_file + '_internal_ancestral_reconstruction.fasta', 'w') as f3,
              open(out_file + '_leaves_evolutionary_events.fasta', 'w') as f2):
            
                no_internal = no_leaves + 1
                for node in tree.traverse('preorder'):
                    if node.name == '':  
                        if node.is_root():
                            node.name = 'ROOT'
                        else:
                            node.name = 'N' + str(no_internal)
                            no_internal += 1
                        
                    ParsAncestral(node)
                    if not node.is_leaf():
                        print('>', node.name,'\n',node.evolutionary_events, file=f1)
                        print('>', node.name, '\n', node.sequence, file=f3)
                    else:
                        print('>', node.name,'\n',node.evolutionary_events, file=f2)
        f1.close()
        f2.close()
        f3.close()
        
        tree.write(format=1, outfile=out_file+'_tree.nwk')
            
    
    return tree_score



