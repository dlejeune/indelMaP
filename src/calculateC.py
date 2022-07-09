#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 19 10:44:39 2022

@author: claraiglhaut
"""

import numpy as np
import math
from scipy.linalg import expm

def calculateC(alphabet, model, branch_length):
    '''
    

    Parameters
    ----------
    alphabet : str
        DNA or Protein.
    Q : np.array
        The symmetric rate matrix Q scaled such that the average substitution 
        rate is 1.
    d : float
        Branch length.
    
    Returns
    -------
    cost_matrix : double dict
        Specifies the cost for each character pair.
    average_cost: float
        The average cost which is used to calculate gap penalties.

    '''
    
    if alphabet == 'Protein':
        characters = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 
              'F', 'P', 'S', 'T', 'W', 'Y', 'V']
        
    elif alphabet == 'DNA':
        characters = ['T', 'C', 'A', 'G']
    
    
    # No rate matrix is given --> use cost of one for substitutions 
    if model == None:
        cost_matrix = {}
        for i in range(len(characters)):
            row = {}
            for j in range(len(characters)):
                if characters[i] == characters[j]:
                    row[characters[j]] = 0
                else:
                    row[characters[j]] = 1
                    
            cost_matrix[characters[i]] = row
            
        average_cost = 1
        C = (cost_matrix, average_cost)
    
    # Calculate the branch lenght specific cost matrix
    elif not branch_length:
        Q = model[0]
        d=0.5
        P = expm(Q*d)
        cost_matrix = {}
        sum_avg = 0
        for i in range(len(P)):
            row = {}
            for j in range(len(P)):
                row[characters[j]] = round(-math.log(P[i][j]))
                sum_avg += round(-math.log(P[i][j]))
                
            cost_matrix[characters[i]] = row
        
        average_cost = sum_avg/(len(characters)*len(characters))
        C = (cost_matrix, average_cost)
        
    else:   
        Q = model[0]
        # Calculate the Approximation for the transition probability matrix
        # if the distances is not given we use the distance 0.62
        C = {}
        for d in [0.1,0.3,0.5,0.7]:
            P = expm(Q*d)
            cost_matrix = {}
            sum_avg = 0
            for i in range(len(P)):
                row = {}
                for j in range(len(P)):
                    row[characters[j]] = round(-math.log(P[i][j]))
                    sum_avg += round(-math.log(P[i][j]))
                    
                cost_matrix[characters[i]] = row
            
            average_cost = sum_avg/(len(characters)*len(characters))
            C[d] = (cost_matrix, average_cost)
        
    return C


def calculateC1(alphabet, model, branch_length, d):
    '''
    

    Parameters
    ----------
    alphabet : str
        DNA or Protein.
    Q : np.array
        The symmetric rate matrix Q scaled such that the average substitution 
        rate is 1.
    d : float
        Branch length.
    
    Returns
    -------
    cost_matrix : double dict
        Specifies the cost for each character pair.
    average_cost: float
        The average cost which is used to calculate gap penalties.

    '''
    
    if alphabet == 'Protein':
        characters = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 
              'F', 'P', 'S', 'T', 'W', 'Y', 'V']
        
    elif alphabet == 'DNA':
        characters = ['T', 'C', 'A', 'G']
    
    
    # No rate matrix is given --> use cost of one for substitutions 
    if model == None:
        cost_matrix = {}
        for i in range(len(characters)):
            row = {}
            for j in range(len(characters)):
                if characters[i] == characters[j]:
                    row[characters[j]] = 0
                else:
                    row[characters[j]] = 1
                    
            cost_matrix[characters[i]] = row
            
        average_cost = 1
        C = (cost_matrix, average_cost)
    
    # Calculate the branch lenght specific cost matrix
    elif not branch_length:
        Q = model[0]
        d=0.5
        P = expm(Q*d)
        cost_matrix = {}
        sum_avg = 0
        for i in range(len(P)):
            row = {}
            for j in range(len(P)):
                row[characters[j]] = round(-math.log(P[i][j]))
                sum_avg += round(-math.log(P[i][j]))
                
            cost_matrix[characters[i]] = row
        
        average_cost = sum_avg/(len(characters)*len(characters))
        C = (cost_matrix, average_cost)
        
    else:   
        Q = model[0]
        # Calculate the Approximation for the transition probability matrix
        # if the distances is not given we use the distance 0.62
        C = {}
        
        P = expm(Q*d)
        cost_matrix = {}
        sum_avg = 0
        for i in range(len(P)):
            row = {}
            for j in range(len(P)):
                row[characters[j]] = round(-math.log(P[i][j]))
                sum_avg += round(-math.log(P[i][j]))
                
            cost_matrix[characters[i]] = row
        
        average_cost = sum_avg/(len(characters)*len(characters))
        C[d] = (cost_matrix, average_cost)
        
    return C
