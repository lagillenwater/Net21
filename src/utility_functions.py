## create a graph from an edgelist

# libraries
import argparse
import networkx as nx
import pandas as pd
import warnings
## arguments.
## 2024-05-28 - include edgelist file and source.

def define_arguments():
    parser=argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument( '-e', dest = "edgelist_file",help="Edgelist file. The format may be dependent upon the network. Default is [entrez gene id 1][entrez gene id 2][posterior prob.] based on the HumanBase full network format")
    parser.add_argument( '-s', dest = "source",default = "HumanBase", help ="Source of the edgelist. May be 'HumanBase' (default).")
    parser.add_argument( '-m', dest = "hgnc_file", help ="HGNC mapping file for converting entrez IDs to HGNC symbols")
    
    return parser

# wrapper function to generate the objects from arguments
def generate_arguments():

    #Generate argument parser and define arguments
    parser = define_arguments()
    args = parser.parse_args()
    edgelist = args.edgelist_file
    source = args.source
    hgnc_file = args.hgnc_file

    return edgelist,source,hgnc_file


# create graph from edgelist file edgelist
def read_edgelist(edgelist, source):
    if source=="HumanBase":
        G = nx.read_edgelist(edgelist, delimiter='\t', data = (("posterior_prob", float),))
        
        return(G)
                             
# print first first 3 edges from graph
def three_edges(G):
    # Get the first three edges and print them
    first_three_edges = list(G.edges(data = True))[:3]
    print("###### First three edges in the graph: ######")
    for edge in first_three_edges:
        print(edge)

## convert a list of entrez ids to gene symbols
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
    

    
