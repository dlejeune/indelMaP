#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  6 15:08:35 2022

@author: claraiglhaut
"""
import random
from enum import Flag, auto
import numpy as np

random.seed(128)


class flags(Flag):
    no_gap = auto()
    gap_opening = auto()
    gap_extension = auto()


def determineT(min_S, S_M, S_X, S_Y):
    '''
    Sets the traceback pointer, in case of a tie the pointer is chosen uniformly 
    at random between the possibilities.

    Parameters
    ----------
    min_S : float
        Minimum value for the scoring matrix.
    S_M : float
        Score reached when matching two columns of the left and right alignment.
    S_X : float
        Score reached when matching the left alignment column with a gap column.
    S_Y : float
        Score reached when matching the right alignment column with a gap column.

    Returns
    -------
    T : int
        Returns the pointer.

    '''
    if (min_S == S_M and 
        min_S == S_X and 
        min_S == S_Y):
        T = random.sample([1,2,3],1)[0]
    elif (min_S == S_M and 
          min_S == S_X):
        T = random.sample([1,3],1)[0]
    elif (min_S == S_M and 
          min_S == S_Y):
        T = random.sample([1,2],1)[0]
    elif (min_S == S_X and 
          min_S == S_Y):
        T = random.sample([2,3],1)[0] 
    elif min_S == S_M:
        T = 1 #move diagonal
    elif min_S == S_X:
        T = 3 #move vertical
    elif min_S == S_Y:
        T = 2 #move horizontal
        
    return T

def setcurrentT(currentT_value, T_M, T_X, T_Y):
    '''
    Changes the current trace back matrix depending on the recorded value.

    Parameters
    ----------
    currentT_value : int
        The current pointer.
    T_M : np.array
        Trace back matrix for matching.
    T_X : np.array
        Trace back matrix for matching the left alignment column with a gap column.
    T_Y : np.array
        Trace back matrix for matching the right alignment column with a gap column..

    Returns
    -------
    current_T : np.array
        The new trace back matrix.

    '''
    if currentT_value == 1:
        current_T = T_M
    elif currentT_value == 2:
        current_T = T_Y
    elif currentT_value == 3:
        current_T = T_X
    return current_T

def gap_traceback(direction_gap, left_flags, right_flags, left_insertions, right_insertions, tmp_i, tmp_j, current_T, T_M, T_X, T_Y):  
    
    if direction_gap == 'left':
        gap = left_flags[tmp_i-1]  
    elif direction_gap == 'right':
        gap = right_flags[tmp_j-1]     
    count = 0
    while gap != flags.gap_opening and tmp_i>0 and tmp_j>0:
        count += 1
        if left_insertions[tmp_i-1] and right_insertions[tmp_j-1]:
            tmp_i -= 1 
            tmp_j -= 1
                        
        elif left_insertions[tmp_i-1]:
            tmp_i -= 1

        elif right_insertions[tmp_j-1]:
            tmp_j -= 1
                        
        #move diagonal 
        elif current_T[tmp_i][tmp_j] == 1:
            break
            
        #move horizontal
        elif np.array_equal(current_T, T_Y):
            current_T = setcurrentT(current_T[tmp_i][tmp_j], T_M, T_X, T_Y)
            tmp_j -= 1 
            
        #move vertical
        elif np.array_equal(current_T, T_X):
            current_T = setcurrentT(current_T[tmp_i][tmp_j], T_M, T_X, T_Y)
            tmp_i -= 1 
    print(count)
    return tmp_i, tmp_j, current_T

def determineC(C_all, dist, gi_f, ge_f, branch_length):
    
    close = 0.1
    intermediate = 0.3
    distant = 0.5
    very_distant = 0.7
    
    if len(C_all) != 2:
        
        if dist == 0.0 or not branch_length:
        
            C = C_all[distant][0]
            gi = gi_f*C_all[distant][1]
            ge = ge_f*C_all[distant][1]
        
        elif 0.0 < dist and dist <= 0.2:
            C = C_all[close][0]
            gi = gi_f*C_all[close][1]
            ge = ge_f*C_all[close][1]
        
        elif 0.2 < dist and dist <= 0.4:
            C = C_all[intermediate][0]
            gi = gi_f*C_all[intermediate][1]
            ge = ge_f*C_all[intermediate][1]
            
        elif 0.4 < dist and dist <= 0.6:
            C = C_all[distant][0]
            gi = gi_f*C_all[distant][1]
            ge = ge_f*C_all[distant][1]
            
        elif 0.6 < dist:
            C = C_all[very_distant][0]
            gi = gi_f*C_all[very_distant][1]
            ge = ge_f*C_all[very_distant][1]
    else:
        C = C_all[0]
        gi = gi_f*C_all[1]
        ge = ge_f*C_all[1]
    
    return C, gi, ge
        
#%%
