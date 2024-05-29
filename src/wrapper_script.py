# libraries and functions
from utility_functions import *
from RWR import RWR

# main function
def main():

    edgelist,source = generate_arguments()
    print("###### Creating graph from edgelist: " + edgelist + " #######")
    G = read_edgelist(edgelist,source)
    print('###### Running random walk with restart (aka, PageRank) ######')
    pr = RWR(G)
    print('###### top RWR scores:  ######')
    print(dict(list(pr.items())[:10]))

    
if __name__ == '__main__':
    main()

