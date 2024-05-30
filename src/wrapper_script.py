# libraries and functions
from utility_functions import *
from RWR import RWR

# main function
def main():

    edgelist,source,hgnc_file = generate_arguments()
    print("###### Creating graph from edgelist: " + edgelist + " #######")
    G = read_edgelist(edgelist,source)
    print('###### Running random walk with restart (aka, PageRank) ######')
    pr = RWR(G)
    print('###### top RWR scores:  ######')

    # Create a new dictionary with replaced keys
    hgnc_results = {k: v for k, v in zip(entrez_to_hgnc(list(pr.keys()),hgnc_file), pr.values())}
    print(dict(list(hgnc_results.items())[:10]))

if __name__ == '__main__':
    main()

