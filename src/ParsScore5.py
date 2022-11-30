#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 13 15:24:28 2022

@author: claraiglhaut
"""


import numpy as np
from ete3 import PhyloNode
from random import sample


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
    
    #initialize list for parsimony sets
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
    
    # initialize parsimony scores for each site
    pars_scores = [0]*len(leaf.sequence)
    
    # add parsimony sets and parsimony scores to each leaf
    leaf.add_features(parsimony_sets = pars_sets)
    leaf.add_features(parsimony_scores = pars_scores)

  
def ParsInternal(tree, i):
    '''
    Creates the parsimony sets and scores for the internal nodes. 

    Parameters
    ----------
    tree : PhyloNode or PhlyoTree
        Internal nodes a tree structure.
    i : int
        Index for the sequence.
    
    Returns
    -------
    None.

    '''
    
    # tree.children[0] accesses the left child for the subtree
    left_set = tree.children[0].parsimony_sets
    left_score = tree.children[0].parsimony_scores[i]
    
    # tree.children[1] accesses the right child for the subtree
    right_set = tree.children[1].parsimony_sets
    right_score  = tree.children[1].parsimony_scores[i]
    
    length_MSA = len(left_set)
    
    if i == 0:
        # initialize sets and parsimony scores for internal node
        # at the beginning of each sequence
        pars_sets = [set()] * length_MSA
        pars_scores = [0] * length_MSA
    
        tree.add_features(parsimony_sets = pars_sets)
        tree.add_features(parsimony_scores = pars_scores)
    
    # if the sets intersect we do not have an evolutionary event
    if left_set[i].intersection(right_set[i]):
        tree.parsimony_sets[i] = left_set[i].intersection(right_set[i])
        tree.parsimony_scores[i] = left_score + right_score
    
    # if the sets do not intersect we count 1 for an evolutionary event
    elif not left_set[i].intersection(right_set[i]):
        tree.parsimony_sets[i] = left_set[i].union(right_set[i])
        tree.parsimony_scores[i] = left_score + right_score + 1
        
def ParsAncestral(tree):
    '''
    Ancestral sequence reconstruction.
    
    Parameters
    ----------
    tree : PhyloNode or PhlyoTree
        Internal nodes a tree structure.
        
    Returns
    -------
    None.
    '''
    
    # initialize sequence
    seq = ''
    
    for i in range(len(tree.parsimony_sets)):
        
        # if we are not at the root and there is a character in the current parsimony set which was reconstructed at the site in the ancestor, we choose that character
        if not tree.is_root() and tree.up.sequence[i] in tree.parsimony_sets[i]:
            character = tree.up.sequence[i]
        # else we make a random choice between the characters in the currents set
        else:
            character = sample(tree.parsimony_sets[i],1)[0]
        
        # the character is added to the sequence
        seq += character
    
    # the sequence is added to the current node
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
    ancestor_reconstruction : boolean, default = True
        If set to true ancestral sites and evolutionary events are reconstructed.

    Returns
    -------
    tree_score : int
        parsimony score for the whole tree.

    '''
    # load tree structure and multiple sequence alignment
    tree = PhyloNode(newick = tree_file, alignment = msa_file, format=1, quoted_node_names=True)
    length_MSA = len(tree.get_leaves()[0].sequence)
    
    
    # initialize sets and scores at the leaves
    no_leaves = 0
    for leaf in tree.iter_leaves():
            ParsLeaf(leaf, alphabet)
            # count leaves to number internal nodes in ancestral reconstruction
            no_leaves += 1

                    
    # find internal sets and scores
    for node in tree.traverse('postorder'):
        if not node.is_leaf():
            for i in range(length_MSA):
                ParsInternal(node, i)

        
    # sum the parsimony scores at the root over the whole MSA length
    # to get the parsimony score for the whole tree
    tree_score = sum(tree.parsimony_scores)
    
    #
    if ancestor_reconstruction:
        # print ancestral reconstruction to out_file_internal_ancestral_reconstruction.fasta
        with open(out_file + '_internal_ancestral_reconstruction.fasta', 'w') as f3:
            no_internal = no_leaves +1
            for node in tree.traverse('preorder'):
                # name internal nodes if they are not already named
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
        # print tree with internal node names to out_file_tree.nwk
        tree.write(format=1, outfile=out_file+'_tree.nwk')
            
    
    return tree_score



