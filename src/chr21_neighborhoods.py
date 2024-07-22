## Gene Neighborhoods of Chr21 Genes
# load the humanbase blood network

import networkx as nx

blood = nx.read_weighted_edgelist('/Users/lucas/OneDrive - The University of Colorado Denver/Projects/Net21/networks/blood_top', delimiter='\t')

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

## convert a list of gene symbols to entrez ids
def entrez_to_hgnc(entrez_ids,hgnc_file):
    hgnc = pd.read_csv(hgnc_file, sep = "\t", low_memory = False)
    hgnc['entrez_id'] = hgnc['entrez_id'].astype(pd.Int64Dtype()).astype(str)
    # filter by list of entrez_ids
    result_df = hgnc[hgnc['entrez_id'].isin(entrez_ids)]
    # Sort the DataFrame based on the 'order' column
    warnings.simplefilter("ignore") ## See the caveats in the documentation: https://pandas.pydata.org/pandas-docs/stable/user_guide/indexing.html#returning-a-view-versus-a-copy  result_df['entrez_id'] = result_df['entrez_id'].map({id: i for i, id in enumerate(entrez_ids)})/Users/lucas/Library/CloudStorage/OneDrive-TheUniversityofColoradoDenver/Projects/Net21/src/utility_functions.py:55: SettingWithCopyWarning: A value is trying to be set on a copy of a slice from a DataFrame
    result_df['entrez_id'] = result_df['entrez_id'].map({id: i for i, id in enumerate(entrez_ids)})
    warnings.simplefilter("ignore") ## See the caveats in the documentation: https://pandas.pydata.org/pandas-docs/stable/user_guide/indexing.html#returning-a-view-versus-a-cop
    result_df.sort_values(by='entrez_id', inplace=True)
    hgnc_symbols = result_df['symbol'].tolist()
    return hgnc_symbols


# ## read in the chromosome annotation file form migsigdb
# chr_gmt = gseapy.read_gmt('/Users/lucas/OneDrive - The University of Colorado Denver/Projects/Net21/data/migsigdb/c1.all.v2023.2.Hs.symbols.gmt')

# ## find chromosome 21 genes
# ## use regex for 'chr21' key 
# pattern = re.compile('chr21')

# # Search for keys that match the pattern
# chr21_dict = {key:value for key, value in chr_gmt.items() if pattern.match(key)}
# chr21_gene_list = list(chr21_dict.values())
# chr21_gene_list = [item for sublist in chr21_gene_list for item in sublist]

# ## read in the hallmark and curated pathway gene sets from migsigdb
# hallmark_gmt = gseapy.read_gmt('/Users/lucas/OneDrive - The University of Colorado Denver/Projects/Net21/data/migsigdb/h.all.v2023.2.Hs.symbols.gmt')

# curated_gmt = gseapy.read_gmt('/Users/lucas/OneDrive - The University of Colorado Denver/Projects/Net21/data/migsigdb/c2.all.v2023.2.Hs.symbols.gmt')

# ## find the overlaps between the list of chr21 genes and each pathway

# def find_overlaps(target_list, source_dict):
#     # Initialize a dictionary to store the overlaps
#     overlaps = {}
    
#     # Convert the target list to a set for efficient comparison
#     target_set = set(target_list)
    
#     # Iterate over the dictionary items
#     for key, values in source_dict.items():
#         # Find the intersection between the list and the current dictionary values
#         overlap = target_set.intersection(values)
#         # Store the overlap in the dictionary
#         overlaps[key] = list(overlap)
    
#     return overlaps

# # Find the overlaps
