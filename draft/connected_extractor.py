#!/usr/bin/env python3

from sys import argv
from igraph import *
import ast






if __name__ == "__main__":
        
    file_frag_name =  argv[1]  # 'frag_dic.txt'
    read_graph = Graph(17000)
    with open(file_frag_name, 'r') as file_frags:
            for frag in file_frags:
                #print(frag)
                frag_dic = ast.literal_eval(frag)
                #frags_dic = {line_num:frag_dic}
                # SNP_frag=[str(i) for i in list(frag_dic.keys())]  # list of SNPs within the fragment
                 # # add a node (vertex) for each SNP within the fragment
                # read_graph.add_vertices(SNP_frag)
                # tuples_SNP = list((x, y) for x in SNP_frag for y in SNP_frag) #
                # edges_list = list(edge for edge in tuples_SNP if edge[0] != edge[1])  # Adding the edges to the graph
                # read_graph.add_edges(edges_list)
                SNP_in_frag=[i for i in list(frag_dic.keys())]  # list of SNPs within the fragment
                # add a node (vertex) for each SNP within the fragment
                #read_graph.add_vertices(SNP_frag)
                tuples_SNP = list((x, y) for x in SNP_in_frag for y in SNP_in_frag) #
                edges_list = list(edge for edge in tuples_SNP if edge[0] != edge[1])  # Adding the edges to the graph
                read_graph.add_edges(edges_list)
                
    print('finish')
    #print(edges_list)
    
    comp=read_graph.components(mode=STRONG)
    comp_atleast=[]
    for i in range(len(comp)):
        if len(comp[i])>1:
            comp_atleast.append(comp[i])
            
    print(len(comp_atleast))
    print(len(comp_atleast[0]))
    
    print('**************')
    print(summary(read_graph))
    print(len(comp_atleast[1]))
    
