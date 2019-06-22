# Tools for motif discovery
This repository contains Python implementations of graph algorithms for fast motif seach in a network. Only compatible with the networkx library (for now). Installation instructions for networkx can be found [here](https://networkx.github.io/documentation/stable/install.html). Examples and tests can be found in **src/scratch.ipynb**.

### Grochow-Kellis subgraph querying algorithm [in development]
[Link to paper](https://link.springer.com/chapter/10.1007/978-3-540-71681-5_7). Use as follows:
  
    import networkx as nx
    import functions as fn
  
    # H is the query graph and G the network to be searched
    G = nx.fast_gnp_random_graph(10, 0.4, seed=12)
    H = nx.fast_gnp_random_graph(5, 0.7, seed=14)
    
    # returns an array containing instances of H in G up to automorphisms of H
    fn.findSubgraphInstances(H, G) 
    
### Onion decomposition algorithm [in development]
[Link to paper](https://www.nature.com/articles/srep31708). Use as follows:

    # G is the network to be onion-decomposed. 
    # returns a structured array of triples containing node number, coreness and onion layer.
    fn.onionDecompose(G)
    
### Shamir-Tsur subtree isomorphism [in development]

