# -*- coding: utf-8 -*-
"""
Created on Tue Sep 30 11:38:21 2025

@author: Maria Jimenez
Postgenomics: HW 5

"""

import pandas as pd
import os
# Load the CSV


def Select_Unique_Genes(path, filename, dest_name):
    filename = "AML_gene.csv"
    file = os.path.join(path, filename)
    
    
    df = pd.read_csv(file, sep=",")
    
    # Drop empty symbols and keep unique
    unique_genes = df["SYMBOL"].dropna().unique()
    
    # Save to new CSV
    pd.DataFrame(unique_genes, columns=["SYMBOL"]).to_csv("unique_genes.csv", index=False)
    return unique_genes

path = r"C:\Users\Maria\Documents\Bioinformatics_Maria\Classes_Fall_2025\Postgenomics\HWS\Lab5"
filename = "AML_gene.csv"

dest_name = filename[:-4]+"_unique.csv"
print(dest_name)
unique_genes =  Select_Unique_Genes(path, filename, dest_name)

# to find the all the paths from a directed graph, we can use a algorithm named
# Depth-First Search (DFS)
def all_paths(graph, start, end, path=[]):
    path = path + [start]
    if start == end:
        return [path]
    if start not in graph:
        return []
    paths = []
    for node in graph[start]:
        if node not in path:  # avoid cycles
            newpaths = all_paths(graph, node, end, path)
            for p in newpaths:
                paths.append(p)
    return paths


# define the directed graph:
graph = {
    'A': ['B'],
    'B': ['C', 'D'],
    'C': ['E', 'F'],
    'D': ['G', 'H'],
    'E': ['K'],
    'F': ['I', 'J'],
    'G': ['O'],
    'H': ['N'],
    'I': ['L'],
    'J': ['M'],
    'K': ['P'],
    'L': ['P'],
    'M': ['P'],
    'N': ['R'],
    'O': ['R'],
    'P': ['Q'],
    'Q': ['R'],
    'R': []
}
start = 'A'
end = 'R'   
    
paths = all_paths(graph, start, end)
print("All paths from", start, "to", end, ":", paths)

# Find the longest path
max_len = max(len(p) for p in paths)
longest_paths = [p for p in paths if len(p) == max_len]

def path_to_edges(path):
    return ", ".join(f"({path[i]},{path[i+1]})" for i in range(len(path)-1))
print("\ntotal of longest paths: ", len(longest_paths) )
for idx, path in enumerate(longest_paths, 1):
    print(f"\nLongest path {idx}: {path_to_edges(path)}")


##############################################################
#              nodes degree                                  #
##############################################################

def count_proteins_above_average(tsv_file, avg_degree):
    """
    Count how many proteins have node_degree greater than avg_degree.

    Parameters:
    tsv_file (str): path to the TSV file with protein node degrees
    avg_degree (float): average node degree from Network Stats

    Returns:
    int: number of proteins with node_degree > avg_degree
    """
    # Load the TSV file
    df = pd.read_csv(tsv_file, sep='\t')

    # Ensure node_degree column is numeric
    df['node_degree'] = pd.to_numeric(df['node_degree'], errors='coerce')

    # Count proteins above average
    count = (df['node_degree'] > avg_degree).sum()

    return count

# Example usage:
network_stats = 8.23  # STRING Network Stats degree: average node degree:8.23

path = r"C:\Users\Maria\Documents\Bioinformatics_Maria\Classes_Fall_2025\Postgenomics\HWS\Lab5"
tsv_file = "string_node_degrees.tsv"  # replace with your file path
file = os.path.join(path, tsv_file)
 
print(f"Average from the Network Stats: ", network_stats)
num_proteins = count_proteins_above_average(file, network_stats)
print(f"Number of proteins that have a node degree above the found average from the Network Stats: {num_proteins}")
