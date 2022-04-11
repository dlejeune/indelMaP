#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 29 22:23:09 2021

@author: claraiglhaut
"""

import numpy as np
from ete3 import PhyloNode
from source.calculateC import calculateC

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

  
def ParsInternal(tree, i, cost_matrix, go, ge, phylogeny_aware):
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
                tree.parsimony_scores[i] = left_score + right_score + ge  
            else:
                tree.parsimony_scores[i] = left_score + right_score + go
        else:
            tree.parsimony_scores[i] = left_score + right_score + go
                
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
                tree.parsimony_scores[i] = left_score + right_score + ge  
            else:
                tree.parsimony_scores[i] = left_score + right_score + go
        else:
            tree.parsimony_scores[i] = left_score + right_score + go
        
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
            
                
    elif not left_set[i].intersection(right_set[i]):
        min_score = np.inf
        for left_character in left_set[i]:
            for right_character in right_set[i]:
                score = cost_matrix[left_character][right_character]
                
                if score < min_score:
                    min_score = score
                
        tree.parsimony_sets[i] = left_set[i].union(right_set[i])
        tree.parsimony_scores[i] = left_score + right_score + min_score
        
    else:
        min_score = np.inf
        for left_character in left_set[i].intersection(right_set[i]):
            for right_character in left_set[i].intersection(right_set[i]):
                
                score = cost_matrix[left_character][right_character]
                
                if score < min_score:
                    min_score = score
                        
        tree.parsimony_sets[i] = left_set[i].intersection(right_set[i])
        tree.parsimony_scores[i] = left_score + right_score + min_score
       

def ParsScore(tree_file, msa_file, alphabet, Q=np.array(None), gi_f=1, ge_f=1, phylogeny_aware=True,
              branch_length = True):
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

    Returns
    -------
    tree_score : int
        parsimony score for the whole tree.

    '''
    
    tree = PhyloNode(newick = tree_file, alignment = msa_file, format=1, quoted_node_names = True)
    length_MSA = len(tree.get_leaves()[0].sequence)
    
    
    if not branch_length:
        C, av_cost = calculateC(alphabet, Q, 0.62, branch_length)
        go = gi_f*av_cost
        ge = ge_f*av_cost
    
    # find insertion points and mark them with an insertion flag set to True
    if phylogeny_aware:
        for i in range(length_MSA):
            leaf_res = []
            for leaf in tree.iter_leaves():
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
                if branch_length:
                    C, av_cost = calculateC(alphabet, Q, node.dist, branch_length)
                    go = gi_f*av_cost
                    ge = ge_f*av_cost
                    
                ParsInternal(node, i, C, go, ge, 
                                     phylogeny_aware)
            #print(node.parsimony_scores)
        
    # sum the parsimony scores at the root over the whole sequence
    tree_score = sum(tree.parsimony_scores)
    
    return tree_score



