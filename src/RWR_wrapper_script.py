# libraries and functions
from utility_functions import *
from RWR import *
import numpy as np

# main function
def main():

    edgelist,source,hgnc_file = generate_arguments()
    print('\n')
    print("###### Creating graph from edgelist: " + edgelist + " #######")
    G = read_edgelist(edgelist,source)
    print('\n')

    print('###### Running random walk with restart (aka, PageRank) ######')
    pr =  RWR(G)
    print('\n')
    print('###### Top RWR scores:  ######')
    # Create a new dictionary with replaced keys
    hgnc_results = {k: v for k, v in zip(entrez_to_hgnc(list(pr.keys()),hgnc_file), pr.values())}
    print(dict(list(hgnc_results.items())[:10]))
    print('\n')
    print('###### Running personalized randome walk with restart (aka personalized PageRank) ######')
    print('\n')

    
    personalized_pr = RWR(G, personalization = {'1':.9})
    print('###### Top Personalized RWR scores:  ######')
    # Create a new dictionary with replaced keys
    personalized_hgnc_results = {k: v for k, v in zip(entrez_to_hgnc(list(personalized_pr.keys()),hgnc_file), personalized_pr.values())}
    print(dict(list(personalized_hgnc_results.items())[:10]))
    
    

if __name__ == '__main__':
    main()

