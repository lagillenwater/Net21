## create a graph from an edgelist

# libraries
import argparse
import networkx as nx
from PyEntrezId import Conversion

## arguments.
## 2024-05-28 - include edgelist file and source.

def define_arguments():
    parser=argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument( '-e', dest = "edgelist_file",help="Edgelist file. The format may be dependent upon the network. Default is [entrez gene id 1][entrez gene id 2][posterior prob.] based on the HumanBase full network format")
    parser.add_argument( '-s', dest = "source",default = "HumanBase", help ="Source of the edgelist. May be 'HumanBase' (default).")
    
    return parser

# wrapper function to generate the objects from arguments
def generate_arguments():

    #Generate argument parser and define arguments
    parser = define_arguments()
    args = parser.parse_args()
    edgelist = args.edgelist_file
    source = args.source

    return edgelist,source


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




def entrez_to_hgnc(entrez_id):
    Id = Conversion('dummyemail@dummybunny.info')
    hgnc_symbol = Id.convert_entrez_to_hgnc(entrez_id)
    return(hgnc_symbol)


