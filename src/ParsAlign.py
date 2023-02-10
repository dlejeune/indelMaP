#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  1 08:55:19 2022

@author: claraiglhaut
"""

import numpy as np
from ete3 import PhyloNode
from calculateC import calculateC
from HelperFunctions import determineT, setcurrentT, determineC
from RateMatrix import WAG, blosum, JC69, K80
import argparse
import os



def InitalizeSetsAndAlignment(leaf, alphabet, indel_aware):
    '''
    Initializes the parsimony sets and the alignments at the leaf nodes
    Parameters
    ----------
    leaf : PhlyoNode or PhyloTree
        Tree leaves with ungapped sequences
    indel_aware : boolean, default = True
        If set to true possible insertions and deletions are treated 
        differently 
        
    Returns
    -------
    None.
    '''
    pars_sets = []
    align = np.empty((1, 0), dtype=str)
    if indel_aware:
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
            
            if indel_aware:
                ins_flags.append(False)
                ins_perm.append(False)
    
    leaf.add_features(parsimony_sets = pars_sets)
    leaf.add_features(alignment = align)
    
    if indel_aware:
        leaf.add_features(insertion_flags = ins_flags)
        leaf.add_features(insertion_permanent = ins_perm)

def GenerateMatrices(tree, C_all, gi_f, ge_f, indel_aware, branch_length):
    '''
    Forward phase of the progressive algorithm. Generates the matrix S with
    the score and the trace back matrix T, which records the moves. 
    Returns the a weighted parsimony score and the trace back matrix.
    
    Parameters
    ----------
    tree : PhlyoTree or PhyloNode
        Current (sub-)tree
    cost_matrix : double dictionary
        Gives the scores for matching, mismatching characters
    gi : float
        Gives gap opening panelty
    ge : float
        Gives gap extension penalty
    indel_aware : boolean, default = True
        If set to true possible insertions and deletions are treated 
        differently 
        
    Returns
    -------
    parsimony_score : numpy.float64
        parsimony score of the alignment for the given tree
    T : numpy.ndarray
        Trace back matrix 
    '''
   
    
    left_sets = tree.children[0].parsimony_sets
    right_sets = tree.children[1].parsimony_sets
    
    if indel_aware:
        left_flags = tree.children[0].insertion_flags
        right_flags = tree.children[1].insertion_flags
        
        left_insertions = tree.children[0].insertion_permanent
        right_insertions = tree.children[1].insertion_permanent
    
    left_dist = tree.children[0].dist
    right_dist = tree.children[1].dist
    
    C_left, gi_left, ge_left = determineC(C_all, left_dist, gi_f, ge_f, 
                                          branch_length)
    C_right, gi_right, ge_right = determineC(C_all, right_dist, gi_f, ge_f, 
                                             branch_length)
        
    S = np.zeros((len(left_sets)+1, len(right_sets)+1), dtype=float)
    T = np.zeros((len(left_sets)+1, len(right_sets)+1), dtype=float)
    
    #fill the first column of the trace back matrix 
    for i in range(1,len(left_sets)+1):
        
        if indel_aware:
            if left_insertions[i-1]:
                S[i][0] = S[i-1][0]                 
            elif left_flags[i-1]:
                S[i][0] = S[i-1][0]             
            else:
                tmp_i = i-1
                while left_flags[tmp_i-1] and tmp_i > 0:
                    tmp_i -= 1
                    
                if tmp_i > 0 and S[tmp_i][0] != 0 and not left_flags[tmp_i-1]:
                    S[i][0] = S[i-1][0] + ge_left               
                else:
                    S[i][0] = S[i-1][0] + gi_left
        else:
            S[i][0] = gi_left + (i-1)*ge_left
        
        T[i][0] = 3
            
        
        
    #fill the first row of the trace back matrix
    #only match gaps with the left alignment - move horizontal
    for j in range(1,len(right_sets)+1):
        if indel_aware:
            if right_insertions[j-1]:
                S[0][j] = S[0][j-1] 
                
            elif right_flags[j-1]:
                S[0][j] = S[0][j-1] 

            else:
                tmp_j = j-1
                while right_flags[tmp_j-1] and tmp_j > 0:
                    tmp_j -= 1
            
                if tmp_j > 0 and S[0][tmp_j] != 0 and not right_flags[tmp_j-1]:
                    S[0][j] = S[0][j-1]  + ge_right
                else:
                   S[0][j] = S[0][j-1]  + gi_right
        else:
            S[0][j] = gi_right + (j-1)*ge_right 
        
        T[0][j] = 2
            
        
    
    #fill the score matrix S
    for i in range(1, len(left_sets)+1):
        for j in range(1, len(right_sets)+1):
            
            if indel_aware and left_insertions[i-1] and right_insertions[j-1]:
                S[i][j] = S[i-1][j-1]
                T[i][j] = 0 #skip 
            
            elif indel_aware and right_insertions[j-1]:
                S[i][j] = S[i][j-1]
                T[i][j] = 0
                
            elif indel_aware and left_insertions[i-1]:
                S[i][j] = S[i-1][j]
                T[i][j] = 0
            
            
            else:
                #matching sets with a non empty intersection
                if left_sets[i-1].intersection(right_sets[j-1]):
                    new_set = left_sets[i-1].intersection(right_sets[j-1])
                elif not left_sets[i-1].intersection(right_sets[j-1]):
                    new_set = left_sets[i-1].union(right_sets[j-1])
                    
                min_score_left = np.inf
                for left_character in left_sets[i-1]:
                    for int_character in new_set:
                        score = C_left[left_character][int_character]
                        
                        if score < min_score_left:
                            min_score_left = score
                
                min_score_right = np.inf
                for right_character in right_sets[j-1]:
                    for int_character in new_set:
                        score = C_right[right_character][int_character]
                        
                        if score < min_score_right:
                            min_score_right = score
                
                score_intersection = S[i-1][j-1] + min(min_score_left, min_score_right)
                                    
                    
                #if the algorithm is phylogeny aware gap placements in columns 
                #flagged as deletions are free
                if indel_aware and right_flags[j-1]:
                    score_gap_left = S[i][j-1]
                
                #otherwise gaps are penalized with gi for newly openend gaps 
                #and ge for gap extensions
                else:
                    tmp_j = j-1
                    while ((T[i][tmp_j] == 0 and tmp_j > 0)
                           or (indel_aware and T[i][tmp_j] == 2 
                               and right_flags[tmp_j-1] 
                               and tmp_j > 0)):
                        tmp_j -= 1                    
                        
                    if T[i][tmp_j] == 2:
                        if indel_aware:
                            if not right_flags[tmp_j-1]:
                                score_gap_left = S[i][j-1] + ge_right
                            else:
                                score_gap_left = S[i][j-1] + gi_right
                        else:
                            score_gap_left = S[i][j-1] + ge_right
                    else:
                        score_gap_left = S[i][j-1] + gi_right
                
                if indel_aware and left_flags[i-1]:
                    score_gap_right = S[i-1][j]
                
                else:
                    tmp_i = i-1
                      
                    while ((T[tmp_i][j] == 0 and tmp_i > 0)
                           or (indel_aware and T[tmp_i][j] == 3 and
                               left_flags[tmp_i-1] and tmp_i > 0)):
                        tmp_i -= 1
     
                    if T[tmp_i][j] == 3:
                        if indel_aware:
                            if not left_flags[tmp_i-1]:
                                score_gap_right = S[i-1][j] + ge_left
                            else:
                                score_gap_right = S[i-1][j] + gi_left
                        else:
                            score_gap_right = S[i-1][j] + ge_left
                    else:
                        score_gap_right = S[i-1][j] + gi_left
                    
                S[i][j] = min(score_intersection, score_gap_left, 
                              score_gap_right)
                
                T[i][j] = determineT(S[i][j], score_intersection, score_gap_right,
                                     score_gap_left)

    
    pars_score = S[len(left_sets)][len(right_sets)]
    #print(S)
    return pars_score, T
        
def GenerateMatricesAffine(tree, alphabet, C_all, gi_f, ge_f, indel_aware, branch_length):
    '''
    Forward phase of the progressive algorithm. Generates the matrix S with
    the score and the trace back matrix T, which records the moves. 
    Returns the a weighted parsimony score and the trace back matrix.
    
    Parameters
    ----------
    tree : PhlyoTree or PhyloNode
        Current (sub-)tree
    cost_matrix : double dictionary
        Gives the scores for matching, mismatching characters
    gi : float
        Gives gap opening penalty
    ge : float
        Gives gap extension penalty
    indel_aware : boolean, default = True
        If set to true possible insertions and deletions are treated 
        differently 
        
    Returns
    -------
    parsimony_score : numpy.float64
        parsimony score of the alignment for the given tree
    T : numpy.ndarray
        Trace back matrix 
    '''
   
    
    left_sets = tree.children[0].parsimony_sets
    right_sets = tree.children[1].parsimony_sets
    
    if indel_aware:
        left_flags = tree.children[0].insertion_flags
        right_flags = tree.children[1].insertion_flags
        
        left_insertions = tree.children[0].insertion_permanent
        right_insertions = tree.children[1].insertion_permanent
    
    
    left_dist = tree.children[0].dist
    right_dist = tree.children[1].dist
    
    C_left, gi_left, ge_left = determineC(C_all, left_dist, gi_f, ge_f,
                                          branch_length)
    C_right, gi_right, ge_right = determineC(C_all, right_dist, gi_f, ge_f,
                                             branch_length)
    
    S_M = np.zeros((len(left_sets)+1, len(right_sets)+1), dtype=float)
    S_Y = np.zeros((len(left_sets)+1, len(right_sets)+1), dtype=float)
    S_X = np.zeros((len(left_sets)+1, len(right_sets)+1), dtype=float)
    T_M = np.zeros((len(left_sets)+1, len(right_sets)+1), dtype=float)
    T_Y = np.zeros((len(left_sets)+1, len(right_sets)+1), dtype=float)
    T_X = np.zeros((len(left_sets)+1, len(right_sets)+1), dtype=float)
    
    #fill the first column of the trace back matrix 
    for i in range(1,len(left_sets)+1):
        S_M[i][0] = np.inf    
        S_Y[i][0] = np.inf
        if indel_aware:
            if left_insertions[i-1] or left_flags[i-1]:
                S_X[i][0] = S_X[i-1][0]        
            else:
                if S_X[i-1][0] != 0:
                    S_X[i][0] = S_X[i-1][0] + ge_right                
                else:
                    S_X[i][0] = S_X[i-1][0] + gi_right
        else:
            S_X[i][0] = gi_right + (i-1)*ge_right
        
        T_M[i][0] = T_X[i][0] = T_Y[i][0] = 3
            
        
        
    #fill the first row of the trace back matrix
    #only match gaps with the left alignment - move horizontal
    for j in range(1,len(right_sets)+1):
        S_M[0][j] = S_X[0][j] = np.inf
        if indel_aware:
            if right_insertions[j-1] or right_flags[j-1]:
                S_Y[0][j] = S_Y[0][j-1] 
            
            else:
                if S_Y[0][j-1] != 0:
                    S_Y[0][j] = S_Y[0][j-1] + ge_left
                else:
                   S_Y[0][j] = S_Y[0][j-1] + gi_left
        else:
            S_Y[0][j] = gi_left + (j-1)*ge_left
        
        T_M[0][j] = T_X[0][j] = T_Y[0][j] = 2
            
           
    #fill the score matrix S
    for i in range(1, len(left_sets)+1):
        for j in range(1, len(right_sets)+1):
            
            if indel_aware and left_insertions[i-1] and right_insertions[j-1]:
                S_M[i][j] = S_M[i-1][j-1]
                S_X[i][j] = S_X[i-1][j-1]
                S_Y[i][j] = S_Y[i-1][j-1]
                
                T_M[i][j] = T_X[i][j] = T_Y[i][j] = 0 #skip 
            
            elif indel_aware and right_insertions[j-1]:
                S_M[i][j] = S_M[i][j-1]
                S_X[i][j] = S_X[i][j-1]
                S_Y[i][j] = S_Y[i][j-1]
                
                T_M[i][j] = T_X[i][j] = T_Y[i][j] = 0 #skip 
                
            elif indel_aware and left_insertions[i-1]:
                S_M[i][j] = S_M[i-1][j]
                S_X[i][j] = S_X[i-1][j]
                S_Y[i][j] = S_Y[i-1][j]
                
                T_M[i][j] = T_X[i][j] = T_Y[i][j] = 0 #skip 
            
            
            else:
                min_S_M = min(S_M[i-1][j-1], S_X[i-1][j-1], S_Y[i-1][j-1])
                
                
                #matching sets with a non empty intersection
                if left_sets[i-1].intersection(right_sets[j-1]):
                    new_set = left_sets[i-1].intersection(right_sets[j-1])
                
                elif not left_sets[i-1].intersection(right_sets[j-1]):
                    new_set = left_sets[i-1].union(right_sets[j-1])
                    
                min_score = np.inf
                for left_character in left_sets[i-1]:
                    for int_character in new_set:
                        score_l = C_left[left_character][int_character]
                        
                        for right_character in right_sets[j-1]:
                            score_r = C_right[right_character][int_character]
                            
                            score = score_l + score_r
                            
                            if score < min_score:
                                min_score = score
                
                S_M[i][j] = min_S_M + min_score
                
                T_M[i][j] = determineT(min_S_M, S_M[i-1][j-1], S_X[i-1][j-1], S_Y[i-1][j-1])
                    
                #if the algorithm is phylogeny aware gap placements in columns 
                #flagged as deletions are free
                if indel_aware and right_flags[j-1]:
                    S_Y[i][j] = min(S_M[i][j-1], S_X[i][j-1], S_Y[i][j-1])
                    T_Y[i][j] = determineT(S_Y[i][j], S_M[i][j-1], S_X[i][j-1], S_Y[i][j-1])
                #otherwise gaps are penalized with gi for newly openend gaps 
                #and ge for gap extensions
                else:
                    tmp = j-1
                
                    while tmp > 0 and right_flags[tmp-1] and (T_Y[i][tmp] == 2 or T_Y[i][tmp] == 0) and indel_aware:
                        tmp -= 1
                    
                    if (T_Y[i][tmp] == 1 or T_Y[i][tmp] == 3) and right_flags[tmp-1]:
                        s_y = S_Y[i][tmp] + gi_left
                    else: 
                        s_y = S_Y[i][tmp] + ge_left
                    s_x = S_X[i][j-1] + gi_left
                    s_m = S_M[i][j-1] + gi_left
                    
                    S_Y[i][j] = min(s_m, s_x, s_y)
                    T_Y[i][j] = determineT(S_Y[i][j], s_m, s_x, s_y)
                
                if indel_aware and left_flags[i-1]:
                    
                    S_X[i][j] = min(S_M[i-1][j], S_Y[i-1][j], S_X[i-1][j])
                    T_X[i][j] = determineT(S_X[i][j], S_M[i-1][j], S_X[i-1][j], S_Y[i-1][j])
                else:
                    tmp = i-1
                    while tmp > 0  and left_flags[tmp-1] and (T_X[tmp][j] == 3 or T_X[tmp][j] == 0) and indel_aware:
                        tmp -= 1
                    
                    if (T_X[tmp][j] == 1 or T_X[tmp][j] == 2) and left_flags[tmp-1]:
                        s_x = S_X[tmp][j] + gi_right
                    else:
                        s_x = S_X[tmp][j] + ge_right
                    s_m = S_M[i-1][j] + gi_right
                    s_y = S_Y[i-1][j] + gi_right
                    
                    S_X[i][j] = min(s_m, s_x, s_y)
                    T_X[i][j] = determineT(S_X[i][j], s_m, s_x, s_y)
                
    
    pars_score = min(S_M[len(left_sets)][len(right_sets)], 
                     S_X[len(left_sets)][len(right_sets)],
                     S_Y[len(left_sets)][len(right_sets)])
    
    start_T = determineT(pars_score, S_M[len(left_sets)][len(right_sets)], 
                     S_X[len(left_sets)][len(right_sets)],
                     S_Y[len(left_sets)][len(right_sets)])
    
    
    
    
    #print(pars_score, '\n', S_M, '\n', T_M, '\n', S_X, '\n', T_X, '\n', S_Y, '\n', T_Y)
    
    return pars_score, T_M, T_X, T_Y, start_T

def TraceBack(T, tree, indel_aware):
    '''
    Finds the alignment for the (sub-)tree and adds it to the (sub-)tree root. 
    Adds the parsimony sets to the (sub-)tree root.
    Parameters
    ----------
    T : numpy.ndarray
        Trace back matrix for the alognment
    tree : PhyloTree or PhyloNode
        Current (sub-)tree
    indel_aware : boolean, default = True
        If set to true possible insertions and deletions are treated 
        differently 
        
    Returns
    -------
    None.
    '''
    
    #get the alignemts from the left and right child
    left_alignment = tree.children[0].alignment
    right_alignment = tree.children[1].alignment
    
    #get the parsimony sets from the left and right child
    left_sets = tree.children[0].parsimony_sets
    right_sets = tree.children[1].parsimony_sets
    
    if indel_aware:
        left_flags = tree.children[0].insertion_flags
        right_flags = tree.children[1].insertion_flags
        
        left_insertions = tree.children[0].insertion_permanent
        right_insertions = tree.children[1].insertion_permanent

    number_of_rows = len(left_alignment) + len(right_alignment)
    
    align = np.empty((number_of_rows, 0), dtype=str)
    
    if indel_aware:
        ins_flags = []
        ins_perm = []
    
    i = left_alignment.shape[1]
    j = right_alignment.shape[1]
    
    pars_sets = []
    while i > 0 or j > 0:
        #permanent insertions 
        if T[i][j] == 0 and indel_aware:
            if left_insertions[i-1] and right_insertions[j-1]:
                
                new_col0 = np.empty((number_of_rows,1), dtype=str)
                for n in range(len(left_alignment)):
                    new_col0[n] = left_alignment[n][i-1]
                    
                for m in range(len(right_alignment)):
                    new_col0[len(left_alignment)+m] = '-' 
          
                pars_sets.insert(0, set('-'))
                ins_flags.insert(0,True)
                ins_perm.insert(0,True)
                    
                i = i-1 
                
                new_col1 = np.empty((number_of_rows,1), dtype=str)
                for n in range(len(left_alignment)):
                    new_col1[n] = '-'
                for m in range(len(right_alignment)):
                    new_col1[len(left_alignment)+m] = right_alignment[m][j-1]
               
                pars_sets.insert(0, set('-'))
                ins_flags.insert(0,True)
                ins_perm.insert(0,True)
            
                j = j-1  
                
                
                new_col = np.concatenate((new_col1, new_col0), axis=1)
            
            elif left_insertions[i-1]:
                new_col = np.empty((number_of_rows,1), dtype=str)
                for n in range(len(left_alignment)):
                    new_col[n] = left_alignment[n][i-1]
                    
                for m in range(len(right_alignment)):
                    new_col[len(left_alignment)+m] = '-' 
          
                pars_sets.insert(0, set('-'))
                ins_flags.insert(0,True)
                ins_perm.insert(0,True)
     
                i = i-1  
            
            elif right_insertions[j-1]:
                new_col = np.empty((number_of_rows,1), dtype=str)
                for n in range(len(left_alignment)):
                    new_col[n] = '-'
                for m in range(len(right_alignment)):
                    new_col[len(left_alignment)+m] = right_alignment[m][j-1]
               
                pars_sets.insert(0, set('-'))
                ins_flags.insert(0,True)
                ins_perm.insert(0,True)
                
                j = j-1 
                
        
        #move diagonal - match two colums with residues 
        elif T[i][j] == 1:
            new_col = np.empty((number_of_rows,1), dtype=str)
            for n in range(len(left_alignment)):
                new_col[n] = left_alignment[n][i-1]
            for m in range(len(right_alignment)):
                new_col[len(left_alignment)+m] = right_alignment[m][j-1]
                
            if not left_sets[i-1].intersection(right_sets[j-1]):
                pars_sets.insert(0, left_sets[i-1].union(right_sets[j-1]))
            else:
                pars_sets.insert(0, left_sets[i-1].intersection(right_sets[j-1]))
            
            if indel_aware:
                ins_flags.insert(0, False)
                ins_perm.insert(0,False)
                
            i = i-1
            j = j-1 
        
        #move vertical - put a gap column in the left alignment and match 
        #with the column of the right alignment
        elif T[i][j] == 2:
            new_col = np.empty((number_of_rows,1), dtype=str)
            for n in range(len(left_alignment)):
                new_col[n] = '-'
            for m in range(len(right_alignment)):
                new_col[len(left_alignment)+m] = right_alignment[m][j-1]
           
            
            if indel_aware:
                ins_flags.insert(0, True)
                
                if right_flags[j-1]:
                    ins_perm.insert(0,True)
                    pars_sets.insert(0, set('-'))
                else:
                    ins_perm.insert(0,False)
                    pars_sets.insert(0, right_sets[j-1])
            else:
                pars_sets.insert(0, right_sets[j-1])
            
            j = j-1    
        
        #move horizontal - put a gap column in the right alignment and match 
        #with the column of the left alignment
        elif T[i][j] == 3:
            new_col = np.empty((number_of_rows,1), dtype=str)
            for n in range(len(left_alignment)):
                new_col[n] = left_alignment[n][i-1]
                
            for m in range(len(right_alignment)):
                new_col[len(left_alignment)+m] = '-' 
      
            
            if indel_aware:
                ins_flags.insert(0, True)
                
                if left_flags[i-1]:
                    ins_perm.insert(0,True)
                    pars_sets.insert(0, set('-'))
                else:
                    ins_perm.insert(0,False)
                    pars_sets.insert(0, left_sets[i-1])
            else:
                pars_sets.insert(0, left_sets[i-1])
                
            i = i-1   

        align = np.concatenate((new_col, align), axis=1)
        
    tree.add_features(alignment = align)
    tree.add_features(parsimony_sets = pars_sets)  
    if indel_aware:
        tree.add_features(insertion_flags = ins_flags)
        tree.add_features(insertion_permanent = ins_perm)
        
def TraceBackAffine(T_M, T_X, T_Y, start_T, tree, indel_aware):
    
    '''
    Finds the alignment for the (sub-)tree and adds it to the (sub-)tree root. 
    Adds the parsimony sets to the (sub-)tree root.
    Parameters
    ----------
    T : numpy.ndarray
        Trace back matrix for the alognment
    tree : PhyloTree or PhyloNode
        Current (sub-)tree
    indel_aware : boolean, default = True
        If set to true possible insertions and deletions are treated 
        differently 
        
    Returns
    -------
    None.
    '''
    
    #get the alignemts from the left and right child
    left_alignment = tree.children[0].alignment
    right_alignment = tree.children[1].alignment
    
    #get the parsimony sets from the left and right child
    left_sets = tree.children[0].parsimony_sets
    right_sets = tree.children[1].parsimony_sets
    
    if indel_aware:
        left_flags = tree.children[0].insertion_flags
        right_flags = tree.children[1].insertion_flags
        
        left_insertions = tree.children[0].insertion_permanent
        right_insertions = tree.children[1].insertion_permanent

    number_of_rows = len(left_alignment) + len(right_alignment)
    
    align = np.empty((number_of_rows, 0), dtype=str)
    
    if indel_aware:
        ins_flags = []
        ins_perm = []
    
    i = left_alignment.shape[1]
    j = right_alignment.shape[1]
    
    current_T = setcurrentT(start_T, T_M, T_X, T_Y)
    
    pars_sets = []

    while i > 0 and j > 0:
        
        #permanent insertions 
        if current_T[i][j] == 0 and indel_aware:
            
            if left_insertions[i-1] and right_insertions[j-1]:
                new_col0 = np.empty((number_of_rows,1), dtype=str)
                for n in range(len(left_alignment)):
                    new_col0[n] = left_alignment[n][i-1]
                    
                for m in range(len(right_alignment)):
                    new_col0[len(left_alignment)+m] = '-' 
          
                pars_sets.insert(0, set('-'))
                ins_flags.insert(0,True)
                ins_perm.insert(0,True)
                    
                i = i-1 
                
                new_col1 = np.empty((number_of_rows,1), dtype=str)
                for n in range(len(left_alignment)):
                    new_col1[n] = '-'
                for m in range(len(right_alignment)):
                    new_col1[len(left_alignment)+m] = right_alignment[m][j-1]
               
                pars_sets.insert(0, set('-'))
                ins_flags.insert(0,True)
                ins_perm.insert(0,True)
            
                j = j-1  
                
                
                new_col = np.concatenate((new_col1, new_col0), axis=1)
            
            elif left_insertions[i-1]:
                new_col = np.empty((number_of_rows,1), dtype=str)
                for n in range(len(left_alignment)):
                    new_col[n] = left_alignment[n][i-1]
                    
                for m in range(len(right_alignment)):
                    new_col[len(left_alignment)+m] = '-' 
          
                pars_sets.insert(0, set('-'))
                ins_flags.insert(0,True)
                ins_perm.insert(0,True)
     
                i = i-1  
            
            elif right_insertions[j-1]:
                new_col = np.empty((number_of_rows,1), dtype=str)
                for n in range(len(left_alignment)):
                    new_col[n] = '-'
                for m in range(len(right_alignment)):
                    new_col[len(left_alignment)+m] = right_alignment[m][j-1]
               
                pars_sets.insert(0, set('-'))
                ins_flags.insert(0,True)
                ins_perm.insert(0,True)
                
                j = j-1 
                
        
        #move diagonal - match two colums with residues 
        elif np.array_equal(current_T, T_M):
            new_col = np.empty((number_of_rows,1), dtype=str)
            for n in range(len(left_alignment)):
                new_col[n] = left_alignment[n][i-1]
            for m in range(len(right_alignment)):
                new_col[len(left_alignment)+m] = right_alignment[m][j-1]
                
            if not left_sets[i-1].intersection(right_sets[j-1]):
                pars_sets.insert(0, left_sets[i-1].union(right_sets[j-1]))
            else:
                pars_sets.insert(0, left_sets[i-1].intersection(right_sets[j-1]))
            
            if indel_aware:
                ins_flags.insert(0, False)
                ins_perm.insert(0,False)
            current_T = setcurrentT(current_T[i][j], T_M, T_X, T_Y)
            
            i = i-1
            j = j-1 
            
                
        
        #move vertical - put a gap column in the left alignment and match 
        #with the column of the right alignment
        elif np.array_equal(current_T, T_Y):
            new_col = np.empty((number_of_rows,1), dtype=str)
            for n in range(len(left_alignment)):
                new_col[n] = '-'
            for m in range(len(right_alignment)):
                new_col[len(left_alignment)+m] = right_alignment[m][j-1]
           
            
            if indel_aware:
                ins_flags.insert(0, True)
                
                if right_flags[j-1]:
                    ins_perm.insert(0,True)
                    pars_sets.insert(0, set('-'))
                else:
                    ins_perm.insert(0,False)
                    pars_sets.insert(0, right_sets[j-1])
            else:
                pars_sets.insert(0, right_sets[j-1])
            current_T = setcurrentT(current_T[i][j], T_M, T_X, T_Y)
            
            j = j-1    
        
        #move horizontal - put a gap column in the right alignment and match 
        #with the column of the left alignment
        elif np.array_equal(current_T, T_X):
            new_col = np.empty((number_of_rows,1), dtype=str)
            for n in range(len(left_alignment)):
                new_col[n] = left_alignment[n][i-1]
                
            for m in range(len(right_alignment)):
                new_col[len(left_alignment)+m] = '-' 
      
            
            if indel_aware:
                ins_flags.insert(0, True)
                
                if left_flags[i-1]:
                    ins_perm.insert(0,True)
                    pars_sets.insert(0, set('-'))
                else:
                    ins_perm.insert(0,False)
                    pars_sets.insert(0, left_sets[i-1])
            else:
                pars_sets.insert(0, left_sets[i-1])
            
            current_T = setcurrentT(current_T[i][j], T_M, T_X, T_Y)
            
            i = i-1   

        align = np.concatenate((new_col, align), axis=1)
    
    while i == 0 and j > 0:
        new_col = np.empty((number_of_rows,1), dtype=str)
        for n in range(len(left_alignment)):
            new_col[n] = '-'
        for m in range(len(right_alignment)):
            new_col[len(left_alignment)+m] = right_alignment[m][j-1]
       
        
        if indel_aware:
            ins_flags.insert(0, True)
            
            if right_flags[j-1]:
                ins_perm.insert(0,True)
                pars_sets.insert(0, set('-'))
            else:
                ins_perm.insert(0,False)
                pars_sets.insert(0, right_sets[j-1])
        else:
            pars_sets.insert(0, right_sets[j-1])
        current_T = setcurrentT(current_T[i][j], T_M, T_X, T_Y)
        
        j = j-1
        
        align = np.concatenate((new_col, align), axis=1)

    while i > 0 and j == 0:
        new_col = np.empty((number_of_rows,1), dtype=str)
        for n in range(len(left_alignment)):
            new_col[n] = left_alignment[n][i-1]
            
        for m in range(len(right_alignment)):
            new_col[len(left_alignment)+m] = '-' 
  
        
        if indel_aware:
            ins_flags.insert(0, True)
            
            if left_flags[i-1]:
                ins_perm.insert(0,True)
                pars_sets.insert(0, set('-'))
            else:
                ins_perm.insert(0,False)
                pars_sets.insert(0, left_sets[i-1])
        else:
            pars_sets.insert(0, left_sets[i-1])
        
        current_T = setcurrentT(current_T[i][j], T_M, T_X, T_Y)
        
        i = i-1   
        
        align = np.concatenate((new_col, align), axis=1)
        
        
    tree.add_features(alignment = align)
    tree.add_features(parsimony_sets = pars_sets)  
    if indel_aware:
        tree.add_features(insertion_flags = ins_flags)
        tree.add_features(insertion_permanent = ins_perm)
        

def ParsAlign(sequence_file,  tree_file, alphabet, output_file=os.path.abspath(os.getcwd())+'/msa',
                      Q=None, gi_f=2.5, ge_f=0.5, indel_aware = True,
                      branch_length = True):
    '''
    
    Returns the indel-aware parsimony score and multiple sequence alignment 
    for a given tree and associated sequences.
    
    Parameters
    ----------
    
    sequence_file : str
        Path to file containing unaligned input sequences in fasta format. 
        If aligned sequences are given, gaps are removed and sequences are realigned.
    tree_file : str
        Path to file containing a guide tree to guide the progressive alignment in newick format.
    alphabet : str
        Choose between 'Protein' or 'DNA'.
    output_file : str
        Path and filename for the outputfile without suffix; 
        if left empty the file will save to the current working directory under the filename msa.fasta
    Q : numpy.ndarray
        Choose the substitution model; - Protein: WAG, blosum - DNA: JC69, K80(alpha,beta); 
        if no substitution model is specified each substitution is associated with the same cost. 
        You can specify your own model by giving a symmetric transition rate matrix with an average substitution rate of 1, 
        as two dimensional np.array enclosed in a tuple (np.array,), substitution rates are given as float.
    gi_f : float, default = 2.5
        Gives factor for gap opening panelty.
    ge_f : float, default = 0.5
        Gives factor for gap extension penalty.
    indel_aware : boolean, default = True
        If set to False the algorithm does not distinguish between insertion and deletion events.
    branch_length : boolean, default = True
        If set to False the cost matrix is calculated based on the transition probability matrix for branch length 0.5 for each branch; 
        if True the cost matrix is calculated based on one of the four distances [0.1,0.3,0.5,0.7], the distance closest to the branch length is chosen.
    
    Returns
    -------
    parsimony_score : numpy.float64
        parsimony score for the alignment on the tree
    alignment : numpy.ndarray
        Multiple sequence alignment for the given tree 
    '''
    
    tree = PhyloNode(newick = tree_file, alignment = sequence_file, format=1, quoted_node_names=True) 
    
    C_all = calculateC(alphabet, Q, branch_length)
    
    #Use Gotoh for affine gaps
    pars_s = 0
    if gi_f != ge_f:
           
        for node in tree.traverse('postorder'):
            if node.is_leaf():
                InitalizeSetsAndAlignment(node, alphabet, indel_aware)    
            else:
                pars_score, T_M, T_X, T_Y, start_T = GenerateMatricesAffine(node, alphabet, C_all, gi_f, ge_f,
                                                 indel_aware, branch_length)
                
                TraceBackAffine(T_M, T_X, T_Y, start_T, node, indel_aware)
                pars_s += pars_score
                node.add_features(parsimony_score = pars_score)
            
    
    #Use Needleman and Wunsch for linear gap cost
    else:
        for node in tree.traverse('postorder'):
            if node.is_leaf():
                InitalizeSetsAndAlignment(node, alphabet, indel_aware)    
            else:
                pars_score, T = GenerateMatrices(node, C_all, gi_f, ge_f, 
                                                  indel_aware, branch_length)
                TraceBack(T, node, indel_aware)
                pars_s = pars_s + pars_score
                node.add_features(parsimony_score = pars_score)
        
    alignment = tree.alignment 
    
    with open(output_file + '.fasta', 'w') as f:
        leaf_names = [leaf.name for leaf in tree.iter_leaves()]
        for i in range(len(tree.alignment)):
            print('>' + leaf_names[i] + '\n' + ''.join(tree.alignment[i]) + '\n', file=f)
        f.close()
        
    return pars_s, alignment

def get_arguments_from_CLI():
    parser = argparse.ArgumentParser(prog = 'Indel-aware parsimony alignment', description='The program progressively alignes protein or nucleotide sequences along a given guide tree under indel-aware parsimony.')
    parser.add_argument('--seq_file', '-s', help='Path to file containing unaligned input sequences in fasta format. If aligned sequences are given, gaps are removed and sequences are realigned.', required = True, type=str)
    parser.add_argument('--tree_file', '-t',  help='Path to file containing a guide tree to guide the progressive alignment in newick format.', required = True, type=str)
    parser.add_argument('--output_file', '-o', default=os.path.abspath(os.getcwd())+'/msa', help = 'Path and filename for the outputfile without suffix; if left empty the file will save to the current working directory under the filename msa.fasta', required = False, type=str)
    parser.add_argument('--alphabet', '-a',  help='Specify sequence alphabet, choose between DNA or Protein', type=str, required = True)
    parser.add_argument('--RateMatrix','-q', default = None, type=str, nargs='+', help = 'Choose the substitution model; - Protein: WAG, blosum - DNA: JC69, K80{alpha,beta}; if no substitution model is specified the default for DNA sequences is K80 with transition to transversion ratio of 2 and for Protein sequences WAG; if each substitution should be associated with the same cost you can type None. You can specify your own model by giving a symmetric transition rate matrix with an average substitution rate of 1. Columns have to be separated by a comma and rows with a colon, and transition rates are given as float. Example: -q=-1,0.3333333333333333,0.3333333333333333,0.3333333333333333:0.3333333333333333,-1,0.3333333333333333,0.3333333333333333:0.3333333333333333,0.3333333333333333,-1,0.3333333333333333:0.3333333333333333,0.3333333333333333,0.3333333333333333,-1')
    parser.add_argument('--gap_opening_factor','-go', default = 2.5, help = 'The gap opening cost is given by the gap_opening_factor*average_substitution_cost; default = 2.5', type=float, required = False)
    parser.add_argument('--gap_extension_factor', '-ge', default = 0.5, help = 'The gap extension cost is given by the gap_extension_factor*average_substitution_cost; default = 0.5', type=float, required = False)
    parser.add_argument('--indel-aware', '-ia', default=True, help = 'If set to False the algorithm does not distinguish between insertion and deletion events; default = True', type=bool, required = False)
    parser.add_argument('--branch_length', '-bl', default=True, help = 'If set to False the cost matrix is calculated based on the transition probability matrix for branch length 0.5 for each branch; if True the cost matrix is calculated based on one of the four distances [0.1,0.3,0.5,0.7], the distance closest to the branch length is chosen.', type=bool, required=False)
    args = parser.parse_args()
    return args

def main():
    
    args = get_arguments_from_CLI()
    
    tree_file = args.tree_file
    sequence_file = args.seq_file
    output_file = args.output_file
    alphabet = args.alphabet
    Q = args.RateMatrix

    if Q == 'None':
        q = None
    elif Q == None:
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
    
    parsimony_score, al = ParsAlign(sequence_file, tree_file, alphabet, output_file,
                          q, gi_f, ge_f, indel_aware,
                          branch_length)

    print('\nThe alignment file was saved to: '+ output_file +'.fasta\nIndel-aware parsimony score: '+str(parsimony_score)+'\n')



if __name__ == '__main__':
    main()
    
