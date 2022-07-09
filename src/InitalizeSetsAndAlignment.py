#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 12 12:38:00 2022

@author: claraiglhaut
"""
import numpy as np

def InitalizeSetsAndAlignment(leaf, alphabet, phylogeny_aware):
    '''
    Initializes the parsimony sets and the alignments at the leaf nodes
    Parameters
    ----------
    leaf : PhlyoNode or PhyloTree
        Tree leaves with ungapped sequences
    phylogeny_aware : boolean, default = True
        If set to true possible insertions and deletions are treated 
        differently 
        
    Returns
    -------
    None.
    '''
    pars_sets = []
    align = np.empty((1, 0), dtype=str)
    if phylogeny_aware:
        ins_flags = []
        ins_perm = []
        
    for i in range(len(leaf.sequence)):
        if leaf.sequence[i] != '-':
            character = np.array(leaf.sequence[i]).reshape((1,1))
            align = np.concatenate((align, character), axis=1)
            
            if alphabet == 'Protein':
                characters = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 
                      'F', 'P', 'S', 'T', 'W', 'Y', 'V']
                if leaf.sequence[i] == 'X':
                    pars_sets.append(set(characters))
                
                elif leaf.sequence[i] == 'B':
                    pars_sets.append(set(['D', 'N']))
                
                elif leaf.sequence[i] == 'Z':
                    pars_sets.append(set(['E', 'Q']))
                
                elif leaf.sequence[i] == 'J':
                    pars_sets.append(set(['I', 'L']))
                    
                else:
                    pars_sets.append(set(leaf.sequence[i]))
            
            if alphabet == 'DNA':
                characters = ['T', 'C', 'A', 'G']
                if leaf.sequence[i] == 'X' or leaf.sequence[i] == 'N':
                    pars_sets.append(set(characters))
                    
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
            
            if phylogeny_aware:
                ins_flags.append(False)
                ins_perm.append(False)
    
    leaf.add_features(parsimony_sets = pars_sets)
    leaf.add_features(alignment = align)
    
    if phylogeny_aware:
        leaf.add_features(insertion_flags = ins_flags)
        leaf.add_features(insertion_permanent = ins_perm)