#!/usr/bin/env python3
#import networkx as nx
from numpy import linalg as LA
import numpy as np
#import matplotlib.pyplot as plt
import ast
from sys import argv

from numpy import linalg as LA
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import ast



if __name__ == "__main__":
    file_frag_name = argv[1]   #file_frag_name ='frag2_edit_dic_a1.txt' #'old/frag_dic.txt'
    
    read_graph=nx.Graph()
    
    with open(file_frag_name, 'r') as file_frags:
            for frag in file_frags:
                frag_dic = ast.literal_eval(frag)
                read_graph.add_nodes_from(frag_dic.keys())  
                snp_in_frag=[i for i in list(frag_dic.keys())]  # list of SNPs within the fragment
                #tuples_SNP = list((x, y) for x in SNP_in_frag for y in SNP_in_frag) #
                #edges_list = list(edge for edge in tuples_SNP if edge[0] != edge[1])  # Adding the edges to the graph
                edges_list=[]
                for i in range(len(snp_in_frag)):
                    if i>0:
                        edges_list.append((snp_in_frag[i-1],snp_in_frag[i]))
                
                read_graph.add_edges_from(edges_list)
    print(len(read_graph)) 
    print(nx.is_connected(read_graph))
    nx.write_edgelist(read_graph, "frag.edgelist")
    
