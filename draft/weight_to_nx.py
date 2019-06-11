#!/usr/bin/env python3
#import networkx as nx
from numpy import linalg as LA
import numpy as np
#import matplotlib.pyplot as plt
import ast
from sys import argv
import networkx as nx



if __name__ == "__main__":
    file_frag_name = argv[1]  # numpy weight matrix of frag_weak_splitter_v2d.py
    W=np.load(file_frag_name) 
    
    read_graph=nx.convert_matrix.from_numpy_matrix(W, parallel_edges=False)
    nx.write_edgelist(read_graph, "frag_node_frag.edgelist")

