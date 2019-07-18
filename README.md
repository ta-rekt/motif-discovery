# Tools for motif discovery
This repository contains Python implementations of fast subgraph isomorphism algorithms to be used for searching network motifs. Only compatible with the networkx library (for now). Installation instructions for networkx can be found [here](https://networkx.github.io/documentation/stable/install.html). Examples and tests can be found in [src/scratch.ipynb](https://github.com/UnifyAndConquer/motif-discovery/blob/master/src/scratch.ipynb).

### Onion decomposition algorithm
[Link to paper](https://www.nature.com/articles/srep31708). Use as follows:

    import functions as fn

    # G is the network to be onion-decomposed. 
    G = nx.fast_gnp_random_graph(10, 0.4, seed=12)
    
    # returns a dictionnary object containing node number, coreness and onion layer.
    fn.onionDecompose(G)
    
### Grochow-Kellis subgraph querying algorithm
[Link to paper](https://link.springer.com/chapter/10.1007/978-3-540-71681-5_7). Use as follows:
  
    import networkx as nx
  
    # H is the query graph and G the network to be searched
    G = nx.fast_gnp_random_graph(10, 0.4, seed=12)
    H = nx.fast_gnp_random_graph(5, 0.7, seed=14)
    
    # returns an array of objects of type Map from H to subgraphs of G up to symmetry.
    # setting the third argument to `True` enables symmetry-breaking.
    instances = fn.findSubgraphInstances(H, G, True)
    
    # to view the mapped nodes
    [k.getMap() for k in instances]
    

### Counting rooted subtrees
This function counts the number of subtrees isomorphic to a rooted query tree T in a rooted tree S. Counting full 5-ary trees in 10-ary trees of sizes 120 and 500 respectively takes about 2 mins. Full k-ary trees are the worst case in terms of running time, which is super-exponential in k because the bulk of the computation consists in computing matrix permanents. However, performance can be greatly improved by using a more efficient implementation of the matrix permanent algorithm. In spite of this, the function does very well on random trees: on a random query tree of size 120 and a search tree of size 600, the number of subtrees is found in less than 9 seconds.

    # T and S are respectively full 5- and 10-ary trees rooted at 0
    T = nx.full_rary_tree(5, 60)
    S = nx.full_rary_tree(10, 200)
    
    fn.countRootedSubtrees(T, 0, S, 0)
