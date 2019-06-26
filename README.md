# Tools for motif discovery
This repository contains Python implementations of graph algorithms for fast motif seach in a network. Only compatible with the networkx library (for now). Installation instructions for networkx can be found [here](https://networkx.github.io/documentation/stable/install.html). Examples and tests can be found in [src/scratch.ipynb](https://github.com/UnifyAndConquer/motif-discovery/blob/master/src/scratch.ipynb).

### Onion decomposition algorithm
[Link to paper](https://www.nature.com/articles/srep31708). Use as follows:

    import functions as fn

    # G is the network to be onion-decomposed. 
    G = nx.fast_gnp_random_graph(10, 0.4, seed=12)
    
    # returns a dictionnary object containing node number, coreness and onion layer.
    fn.onionDecompose(G)
    
### Grochow-Kellis subgraph querying algorithm [in development]
[Link to paper](https://link.springer.com/chapter/10.1007/978-3-540-71681-5_7). Use as follows:
  
    import networkx as nx
  
    # H is the query graph and G the network to be searched
    G = nx.fast_gnp_random_graph(10, 0.4, seed=12)
    H = nx.fast_gnp_random_graph(5, 0.7, seed=14)
    
    # returns an array of objects of type Map from H to subgraphs of G up to symmetry
    instances = fn.findSubgraphInstances(H, G)
    
    # to view the mapped nodes
    [k.getMap() for k in instances]
    
### Shamir-Tsur subtree isomorphism [in development]

