#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  8 17:38:28 2023

@author: ahocampo
"""

"""
Gene python assignment.
"""

import pandas as pd
from itertools import product


def find_all_sizes(gene = pd.Series(['ATTTGGATT'])):
    """Purpose is to find all of the k values needed for a gene sequence.

    Parameters
    ----------
    gene (series): Series containing the gene sequence. 
        The default is pd.Series(['ATTTGGATT']).

    Returns
    -------
    k_list (list): 
        List of possible sizes for the gene sequence. 

    """
    # Create a list of the size of the gene sequences in series.
    max_size_list = []
    for i in range(0, len(gene)):
        max_size_list.append(len(gene[i]))
    # Create empty lists for loops.
    k_list = []
    max_list= []
    gene_list = []
    # Loop through every gene length in the series.
    for i in range(0, len(max_size_list)):
        for j in range(0, max_size_list[i]):
            k_list.append(j+1)
            max_list.append(max_size_list[i])
            gene_list.append(gene[i])
    return k_list


def all_combos(letters = ['A', 'G', 'C', 'T'], k = 2):
    """Find all possible cominations of letters of size k.

    Parameters
    ----------
    letters (list): 
        The default is ['A', 'G', 'C', 'T'].
    k (integer):
        Size of combinations. The default is 2.

    Returns
    -------
    combos (DataFrame): 
        Contains all combinations of letters of size k. 
    """
    
    # Find all combinations of the letters list using replacement. 
    combos_list = product(letters, repeat=k)
    # Concatenate the tuples in the combos_list and return a Series.
    combos = pd.DataFrame(combos_list).agg(''.join, axis=1)
    return combos


        
    
def create_all_size_combos(k_list = [1,2,3], gene = pd.Series(['ATTTGGATT', 'ATGTCTGTCTGTA']), letters = ['A', 'G', 'C', 'T']):
    """Find all combinations for all possible sizes.

    Parameters
    ----------
    k_list (list):
        List of all possible sizes of sequences. The default is [1,2,3].
    gene (series): 
        Series with genetic sequence to analyze. The default
        is pd.Series(['ATTTGGATT', 'ATGTCTGTCTGTA']).
    letters (list): 
        The default is ['A', 'G', 'C', 'T'].

    Returns
    -------
    combos_df  (DataFrame):
        Contains all combinations of letters for various size k.

    """
    combos_df_list = []
    for i in range(0, len(k_list)):
        combos = all_combos(letters = letters, k = k_list[i])
        combos_df_list.append(combos)
    combos_df = pd.concat(combos_df_list)
    return combos_df


def ifexists_and_count(combos = ['AT', 'GC'], gene = pd.Series(['ATGTCTGTCTGTA'])):
    """Function 1. Find number of times a pattern occurs and if it occured.
    

    Parameters
    ----------
    combos (list): 
        List of patterns you want to find. The default is ['AT', 'GC'].
    gene (series): Series containing the gene sequence. 
        The default is pd.Series(['ATTTGGATT']).

    Returns
    -------
    df (DataFrame):
        Contains the gene strings, the pattern looked for, the number of times
        the pattern occured, and a binary column to flag if it existed in the
        gene string.
    """
    # Make the series a list. 
    gene_series_tolist = gene.to_list()
    
    # Create empty lists to iterate through with the loops.
    combo_list = []
    count_list = []
    combo_found_list = []
    gene_list = []
    
    # Loop through each string in the gene list.
    for j in range(0, len(gene_series_tolist)):
        gene_value = gene_series_tolist[j]
        # Count how many times the pattern occurs in the string.
        for i in range(0, len(combos)):
            count = gene_value.count(combos[i])
            count_list.append(count)
            combo_list.append(combos[i])
            gene_list.append(gene_value)
            # Create binary value to flag if the pattern exists in the string.
            if count > 0:
                combo_found_list.append(1)
            else:
                combo_found_list.append(0)
    # Create a dataframe to return from the function with the information about
    # the pattern.
    df = pd.DataFrame()
    df['GeneString'] = gene_list
    df['SearchFor'] = combo_list
    df['Count'] = count_list
    df['IfExists'] = combo_found_list
    search_len_list = []
    possible_strings = []
    for i in range(0, len(df)):
        search_len_list.append(len(combo_list[i]))
        if len(gene_list[i]) == len(combo_list[i]):
            if combo_found_list[i] == 1:
                possible_strings.append(1)
            else:
                possible_strings.append(0)
        else:
            possible_strings.append(1)
    df['SearchLen'] = search_len_list
    df['Possible'] = possible_strings
    return df



def possible(k, glen = 9):
    """Function 2. Find the number of possible substrings. 

    Parameters
    ----------
    k (integer):
        Size of combinations. The default is 2.
    glen (integer):
        Length of gene sequence. The default is 9.

    Returns
    -------
    p (integer):
        Possible substrings in a sequence with a substring length of k.
    """
    # If 4 to the power of k is less than length of g.
    if 4**k < glen:
        # If k is 1 then assign p as 4.
        if k == 1:
            p = 4
        # If p is not 1 and 4^k is less than length of the gene sequence, then
        # p is length of gene sequence^k.
        else:
            x = glen^k
            p = x
    # If 4^k is greater than length of gene sequence, calculate p.
    else:
        x = glen - (k - 1)
        p = x
    return p
    

def main_function_genes():
    """Function 3 and main function. Calculates linguistic complexity and the
    number of possible substrings and the number of observed substrings in an
    inputted gene sequence.

    Returns
    -------
    dffinal (DataFrame):
        Contains the summarized data for the gene sequence and provides the
        linguistic complexity.
    observedtable (DataFrame):
        Contains the summarized data for each k searched in a gene sequence.
    """
    
    # User enters genetic sequence.
    gene_input = input("Please enter the sequence: ")
    gene = pd.Series([gene_input])
    
    # Calculate the counts for the gene sequence for all possible combinations.
    letters = ['A', 'G', 'C', 'T']
    k_list = find_all_sizes(gene=gene)
    combos_df = create_all_size_combos(k_list = k_list, gene = gene)
    df = ifexists_and_count(combos = combos_df.to_list(),
                                      gene = gene)
    # Group the data to see how many observed substrings were found.
    final = df.groupby(by=['GeneString', 'SearchLen']).sum(numeric_only=True).reset_index()
    genestringlist = final['GeneString'].to_list()
    searchlenlist = final['SearchLen'].to_list()
    # Create a column to show how many possible substrings exists.
    z_possible = []
    for i in range(0, len(final)):
        z = possible(searchlenlist[i], len(genestringlist[i]))
        z_possible.append(z)
    final['Possible'] = z_possible
    
    # Subset the data and rename the columns.
    observedtable = final[['GeneString', 'SearchLen', 'IfExists', 'Possible']]
    observedtable.columns = ['Gene', 'k', 'Observed', 'Possible']
    
    # Calculate the linguistic complexity of the gene sequence.
    final = final[['GeneString', 'IfExists', 'Possible']]
    dffinal = final.groupby(by=['GeneString']).sum(numeric_only=True).reset_index()
    dffinal['LinguisticComplexity'] = dffinal['IfExists']/dffinal['Possible']
    dffinal.columns = ['Genetic Sequence', 'Observed Substrings', 'Possible Substrings', 'Linguistic Complexity']
    
    # Print findings.
    vargene = dffinal.iloc[0, 0]
    varobserved = dffinal.iloc[0, 1]
    varpossible = dffinal.iloc[0, 2]
    varcomplexity = dffinal.iloc[0, 3]
    print("Genetic Sequence: ", vargene)
    viewobserved = observedtable[['k', 'Observed', 'Possible']]
    print(viewobserved)
    print("Observed Substrings: ", varobserved)
    print("Possible Substrings: ", varpossible)
    print("Linquistic Complexity: ", varcomplexity)
    return dffinal, observedtable

# dffinal, observedtable = main_function_genes()
