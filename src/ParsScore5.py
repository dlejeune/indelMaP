#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 13 15:24:28 2022

@author: claraiglhaut
"""


import numpy as np
from ete3 import PhyloNode
from random import sample
from src.calculateC import calculateC
from src.HelperFunctions import determineC

def ParsLeaf(leaf, alphabet):
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

    leaf.add_features(parsimony_sets = pars_sets)
    leaf.add_features(parsimony_scores = pars_scores)

  
def ParsInternal(tree, i):
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
    
        tree.add_features(parsimony_sets = pars_sets)
        tree.add_features(parsimony_scores = pars_scores)
            
    if left_set[i].intersection(right_set[i]):
        tree.parsimony_sets[i] = left_set[i].intersection(right_set[i])
        tree.parsimony_scores[i] = left_score + right_score
        
    elif not left_set[i].intersection(right_set[i]):
        tree.parsimony_sets[i] = left_set[i].union(right_set[i])
        tree.parsimony_scores[i] = left_score + right_score + 1
        
def ParsAncestral(tree):
    seq = ''
    
    for i in range(len(tree.parsimony_sets)):
        if not tree.is_root() and tree.up.sequence[i] in tree.parsimony_sets[i]:
            character = tree.up.sequence[i]
        else:
            character = sample(tree.parsimony_sets[i],1)[0]
    
        seq += character
        
    tree.add_features(sequence = seq)
                    
                
def ParsScore5(tree_file, msa_file, out_file, alphabet, ancestor_reconstruction=True):
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
    
        
    no_leaves = 0
    for leaf in tree.iter_leaves():
            ParsLeaf(leaf, alphabet)
            no_leaves += 1

                    
    #find internal sets and scores
    for node in tree.traverse('postorder'):
        if not node.is_leaf():
            for i in range(length_MSA):
                    
                ParsInternal(node, i)

        
    # sum the parsimony scores at the root over the whole sequence
    tree_score = sum(tree.parsimony_scores)
    
    if ancestor_reconstruction:
        with open(out_file + '_internal_ancestral_reconstruction.fasta', 'w') as f3:
            no_internal = no_leaves +1
            for node in tree.traverse('preorder'):
                if node.name == '':  
                    if node.is_root():
                        node.name = 'ROOT'
                    else:
                        node.name = 'N' + str(no_internal)
                        no_internal += 1
                    
                ParsAncestral(node)
                if not node.is_leaf():
                    print('>', node.name, '\n', node.sequence, file=f3)
        
        f3.close()
        tree.write(format=1, outfile=out_file+'_tree.nwk')
            
    
    return tree_score



