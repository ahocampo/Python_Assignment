#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  8 17:45:20 2023

@author: ahocampo
"""

import pytest
import os
import pandas as pd
from itertools import product

file_path = '/Users/ahocampo/Documents/GitHub/Python_Assignment'
os.chdir(file_path)
import GeneticSequenceAssignment as gsa



def test_find_all_sizes():
    assert gsa.find_all_sizes(gene = pd.Series(['ATTTGGATT'])) == [1, 2, 3, 4, 5, 6, 7, 8, 9]

def test_all_combos():
    assert len(gsa.all_combos(letters = ['A', 'G', 'C', 'T'], k = 2)) == 16
    
def test_create_all_size_combos():
    assert len(gsa.create_all_size_combos(k_list = [1,2,3], gene = pd.Series(['ATTTGGATT']))) == 84
   
def test_ifexists_and_count():
    assert len(gsa.ifexists_and_count(combos = ['AT', 'GC'], gene = pd.Series(['ATGTCTGTCTGTA']))) 
    
def test_possible():
        assert gsa.possible(9, glen = 9) == 1
        
        


test_find_all_sizes()

test_all_combos()

test_create_all_size_combos()

test_ifexists_and_count()

test_possible()
