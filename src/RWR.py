## program for running random walks with restart over a network

import networkx as nx

def RWR(G):
    pr = nx.pagerank(G, weight = 'posterior_prob')
    # sort
    sorted_items = sorted(pr.items(), key=lambda item: item[1], reverse = True)
    sorted_dict_comp = {k: v for k, v in sorted_items}

    return(sorted_dict_comp)

