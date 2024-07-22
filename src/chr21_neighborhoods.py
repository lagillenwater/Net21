## Gene Neighborhoods of Chr21 Genes
# load the humanbase blood network

import networkx as nx

blood = nx.read_weighted_edgelist('/Users/lucas/OneDrive - The University of Colorado Denver/Projects/Net21/networks/blood_top', delimiter='\t', create_using=nx.DiGraph())

def network_summary(G):
    summary = {
        'Number of nodes': nx.number_of_nodes(G),
        'Number of edges': nx.number_of_edges(G),
        'Average degree': sum(dict(G.degree()).values()) / float(nx.number_of_nodes(G)),
        'Density': nx.density(G),
        'Is directed': nx.is_directed(G),
        'Number of selfloops': nx.number_of_selfloops(G),
        'Is weighted': nx.is_weighted(G)
    }
    return summary

# Get the summary
print(network_summary(blood))

import collections
import matplotlib.pyplot as plt

# Assuming 'G' is your directed NetworkX graph
in_degree_sequence = sorted([d for n, d in blood.in_degree()], reverse=True)  # in-degree sequence
in_degreeCount = collections.Counter(in_degree_sequence)
in_deg, in_cnt = zip(*in_degreeCount.items())

out_degree_sequence = sorted([d for n, d in blood.out_degree()], reverse=True)  # out-degree sequence
out_degreeCount = collections.Counter(out_degree_sequence)
out_deg, out_cnt = zip(*out_degreeCount.items())

plt.clf()
fig, ax = plt.subplots()
plt.bar(out_deg, out_cnt, width=0.80, color='b')

plt.title("Out-Degree Histogram")
plt.ylabel("Count")
plt.xlabel("Out-Degree")
ax.set_xticks([d + 0.4 for d in out_deg])
ax.set_xticklabels(out_deg)

plt.show()

# Repeat the plotting for out-degree histogram


import gseapy
import re
## read in the chromosome annotation file form migsigdb
chr_gmt = gseapy.read_gmt('/Users/lucas/OneDrive - The University of Colorado Denver/Projects/Net21/data/migsigdb/c1.all.v2023.2.Hs.symbols.gmt')

## find chromosome 21 genes
## use regex for 'chr21' key 
pattern = re.compile('chr21')

# Search for keys that match the pattern
chr21_dict = {key:value for key, value in chr_gmt.items() if pattern.match(key)}
chr21_gene_list = list(chr21_dict.values())
chr21_gene_list = [item for sublist in chr21_gene_list for item in sublist]


## the network has entrez ids, so we need to conver the gene symbols from the migsigdbs into entrez ids
from utility_functions import *

chr21_gene_list_entrez = hgnc_to_entrez(chr21_gene_list,'/Users/lucas/OneDrive - The University of Colorado Denver/Projects/Net21/data/hgnc_complete_set.txt')

## Find the nearest neighbors for all the chr21 genes
def find_nearest_neighbors(G, nodes_list):
    neighbors_dict = {}
    for node in nodes_list:
        try:
            neighbors_dict[node] = list(G.neighbors(node))
        except nx.NetworkXError as e:
            print(f"Node {node} not found in the graph.")
    return neighbors_dict

# Find the nearest neighbors
chr21_neighbors = find_nearest_neighbors(blood, chr21_gene_list_entrez)
chr21_neighbors_count = {}
for key, value in chr21_neighbors.items():
    chr21_neighbors_count[key] = len(chr21_neighbors[key])

# Create a bar plot from the dictionary values
plt.clf()
plt.bar(chr21_neighbors_count.keys(), chr21_neighbors_count.values(), color='g')
plt.xlabel("Ch21 genes")
plt.ylabel("number of neighbors")
plt.show()


    

chr21_neighbors_list = list(chr21_neighbors.values())
chr21_neighbors_list = [item for sublist in chr21_neighbors_list for item in sublist]
chr21_neighbors_list = list(dict.fromkeys(chr21_neighbors_list))

print('Number of neighbors of the 599 chr21 genes in migsigdb: ' + str(len(chr21_neighbors_list)))

## Convert entrez ids back to hgnc ids
chr21_neighbors_hgnc = entrez_to_hgnc(chr21_neighbors_list,'/Users/lucas/OneDrive - The University of Colorado Denver/Projects/Net21/data/hgnc_complete_set.txt')

# ## read in the hallmark and curated pathway gene sets from migsigdb
hallmark_gmt = gseapy.read_gmt('/Users/lucas/OneDrive - The University of Colorado Denver/Projects/Net21/data/migsigdb/h.all.v2023.2.Hs.symbols.gmt')

curated_gmt = gseapy.read_gmt('/Users/lucas/OneDrive - The University of Colorado Denver/Projects/Net21/data/migsigdb/c2.all.v2023.2.Hs.symbols.gmt')

## find the overlaps between the list of chr21 genes and each pathway

def find_overlaps(target_list, source_dict):
    # Initialize a dictionary to store the overlaps
    overlaps = {}
    
    # Convert the target list to a set for efficient comparison
    target_set = set(target_list)
    
    # Iterate over the dictionary items
    for key, values in source_dict.items():
        # Find the intersection between the list and the current dictionary values
        overlap = target_set.intersection(values)
        # Store the overlap in the dictionary
        overlaps[key] = list(overlap)
    
    return overlaps

# Find the overlaps
## Hallmark
hallmark_overlaps = find_overlaps(chr21_neighbors_hgnc, hallmark_gmt)

## Find the percentage of the pathway representation
hallmark_overlaps_percentage = {}

          
# Calculate the length of each list value in the dictionary
for key, value in hallmark_overlaps.items():
    hallmark_overlaps_percentage[key] = len(value)/len(hallmark_gmt[key])

hallmark_overlaps_percentage = {k: v for k, v in sorted(hallmark_overlaps_percentage.items(), key=lambda item: item[1], reverse=True)}

# Iterate over the dictionary and print each key-value pair on a new line
for key, value in hallmark_overlaps_percentage.items():
    print(f'{key}: {value}')

# # Filter the dictionary to include only pairs with values greater than the threshold
# hallmark_overlaps_top = dict(list(hallmark_overlaps_count.items())[:10])

# # Create a bar plot
# plt.clf()
# plt.barh(range(len(hallmark_overlaps_top)), list(hallmark_overlaps_top.values()), align='center')
# plt.yticks(range(len(hallmark_overlaps_top)), list(hallmark_overlaps_top.keys()), fontsize = 8)
# plt.xlabel('Overlap')
# plt.title('Chr21 Overlap of Hallmark Gene Sets')
# plt.show()
