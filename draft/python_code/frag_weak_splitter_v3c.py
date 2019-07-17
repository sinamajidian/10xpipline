#!/usr/bin/env python3
#import networkx as nx
from numpy import linalg as LA
import numpy as np
import matplotlib.pyplot as plt
import ast
from sys import argv
import networkx as nx

thresh=.08



def check_is_connected(W):
    read_graph=nx.convert_matrix.from_numpy_matrix(W, parallel_edges=False)
    check=nx.is_connected(read_graph)
    return check


def bipartition_graph(W,ii):
    D=np.diag(W.sum(axis=0))
    D_half =np.diag(np.power(np.diag(D),-.5))# element wise for diagonal entries #L=D-A  #D degree matrix,  A  adjacency matrix
    mat=np.dot(np.dot(D_half,D-W),D_half)
    lam, V = LA.eig(mat)
    second_eig=V[:,1]
    print('second length is ',len(second_eig))
    #print(second_eig[:20])
    #print(second_eig[-20:])
    #print(np.mean(second_eig))
    plt.plot(second_eig)
    file_name_nct='second_eig_vector'+str(ii)+'.png'
    plt.savefig(file_name_nct)
    plt.close()
    #plt.show()    
    
    part1_indecis=np.where(second_eig>0)
    part2_indecis=np.where(second_eig<=0)
    part1_list_ind=list(part1_indecis[0])
    part2_list_ind=list(part2_indecis[0])
    #print('part 1 length ',len(part1_list_ind))
    #print('part 2 length ',len(part2_list_ind))
    W_part1=W[part1_list_ind][:,part1_list_ind] #A[[0,1]] extracts first and secnond row of matrix A
    W_part2=W[part2_list_ind][:,part2_list_ind]
    dic_2parts={'W_p1':W_part1,'W_p2':W_part2,'ind_p1':part1_list_ind,'ind_p2':part2_list_ind }
    return dic_2parts

def Ncut_graph(ind_A,ind_B,W): # eq 2 page 889 shi malik paper
    # W 2d numpy array
    W_sliced=W[ind_A][:, ind_B]
    cut=W_sliced.sum()
    W_restA=W[ind_A] # extracting those rows
    W_restB=W[ind_B]
    assoc_A=W_restA.sum()
    assoc_B=W_restB.sum()
    Ncut_graph=cut/assoc_A+cut/assoc_B
    return Ncut_graph

def parse_frag_calcu_weight(file_pfrag_name):
    with open(file_frag_name, 'r') as file_frags:
            for frag in file_frags:
                frag_dic = ast.literal_eval(frag)
                list_frag_dic.append(frag_dic)
    frag_num=len(list_frag_dic)
    print('number of fragments is ',frag_num)
    W=np.zeros((frag_num,frag_num))
    for i in range(frag_num):
        for j in range(frag_num):  # half an hour
            if i<j:
                f1=list_frag_dic[i]
                f2=list_frag_dic[j]
                shared_snp_num=len(set(f1.keys())&set(f2.keys()))
                shared_snp_allel_num= len(set(f1.items()) & set(f2.items()))
                #if shared_snp_num>0:
                #    W[i,j]=(2*shared_snp_allel_num-shared_snp_num)/shared_snp_num
                #else: W[i,j]=0
                #if W[i,j]<0:W[i,j]=0
                W[i,j]=shared_snp_num
                W[j,i]=W[i,j]
            if i==j: W[i,i]=1
    #print(np.shape(W)) # dimension
    name_out=file_frag_name+'weight_shared'
    np.save(name_out, W)
    return W

def partitioning(W):
    list_branch_working=[W]
    list_branch_glob_indices_working=[list(range(np.shape(W)[0]))]

    list_branch_finished=[]    
    list_branch_glob_indices_finished=[]
    ii=0
    while len(list_branch_working): 
        ii=ii+1
        print('*******',str(ii),'*******')
        print('number of fragments in working is',[np.shape(W1)[0] for W1 in list_branch_working]) 
        print('number of fragments in finsihed is ',[np.shape(W1)[0] for W1 in list_branch_finished]) 

        W=list_branch_working[0]
        glob_indices=list_branch_glob_indices_working[0]

        
        if check_is_connected(W):
            dic_2parts= bipartition_graph(W,ii) # dic_2parts['W_p1'] and dic_2parts['ind_p1']
            ncut=Ncut_graph(dic_2parts['ind_p1'],dic_2parts['ind_p2'],W)
            print('n cut between first and second is third',np.shape(dic_2parts['W_p1'])[0],np.shape(dic_2parts['W_p2'])[0],np.round(ncut,6))

            if ((ncut>thresh) | (np.shape(dic_2parts['W_p1'])[0] ==0) | (np.shape(dic_2parts['W_p2'])[0]==0)):
                list_branch_finished.append(W)
                list_branch_glob_indices_finished.append(glob_indices)
                filename_out='part'+str(len(list_branch_finished))+'.txt'
                with open(filename_out, 'w+') as file_out: 
                    for line_num in glob_indices:
                        file_out.write(str(line_num+1)+'\n')
                print('it goes to finsihed')
            else:
                list_branch_working.append(dic_2parts['W_p1'])
                glob_indices_p1=[glob_indices[i] for i in dic_2parts['ind_p1']] 
                #print('glob ind 1 is',glob_indices_p1)
                list_branch_glob_indices_working.append(glob_indices_p1)
                list_branch_working.append(dic_2parts['W_p2'])
                glob_indices_p2=[glob_indices[i] for i in dic_2parts['ind_p2']] 
                list_branch_glob_indices_working.append(glob_indices_p2)
                #print('glob ind 2 is',glob_indices_p2)
                #print(list_branch_glob_indices_working)
        else:
            list_branch_finished.append(W)
            list_branch_glob_indices_finished.append(glob_indices)
            filename_out='part'+str(len(list_branch_finished))+'.txt'
            with open(filename_out, 'w+') as file_out: 
                for line_num in glob_indices:
                    file_out.write(str(line_num+1)+'\n')
            print('it is not connected and goes to finished')
        
        del list_branch_glob_indices_working[0]
        del list_branch_working[0]
              
    return list_branch_finished


if __name__ == "__main__":
    file_frag_name =argv[1]  # 'old_data/small_exmp2.txt' #'old_data/frag2_edit_dic_a0.txt'# argv[1]   #file_frag_name ='frag2_edit_dic_a1.txt' #'old/frag_dic.txt'
    list_frag_dic=[]
    if file_frag_name[-3:]=='txt':
        W=parse_frag_calcu_weight(file_frag_name)
    if file_frag_name[-3:]=='npy':
        W=np.load(file_frag_name)
    a= partitioning(W)
    print('>>><<<')
    print('number of partitions is ',len(a))
