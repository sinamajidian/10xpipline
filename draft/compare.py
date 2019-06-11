#!/usr/bin/env python3

####usage python3 compare.py head_truth_21 out_algorithm.txt 

"""
This code is for comparing the estimated haplotype blocks with a grand truth.


My first code in python.

sina majdiain
10 Dec 

"""



from sys import argv
import subprocess
import numpy as np

def parse_hap_truth(file):
	""" extracting haplotype from haplotype truth file  in probhap website
	input: file name
	output: a dictionary, label:index , value: allele 0 or 1
	"""
	hap_file = open(file, "r")
	dic = {}
	for line in hap_file:
		line_slice=line.strip().split("\t")
		label = int(line_slice[0])
		dic[label] = int(line_slice[3])
	return dic 


def parse_hap_estimated_blocks(file):
	""" extracting estimated haplotype blocks from  algorithm's output	
	input: file name
	output: dictionary, label: block number,  value: a dictionary of idx and hap
	"""
	hap_block_file = open(file, "r")
	dic = {}
	block_idx=0;
	for line in hap_block_file:
			line=line.strip()
			if line.startswith('B'): # new block started
				block_idx+=1
				dic_in={}
				dic[block_idx]=dic_in
			elif  (not line.startswith('*')) & (not '-' in line): 
				line_slice=line.split("\t")
				idx = int(line_slice[0])
				hp = int(line_slice[1])
				dic_in[idx]=hp
	dic[block_idx]=dic_in
	return dic

def compare(dic_truth,dic_block_estimated):
	""" comaring estimated block of haplotype with truth haplotype
	input: two dicinary
	output: rr
	"""
	correct_num_1=0
	correct_num_2=0
	leng_block_intersct_truth=0;
	estimd_snp_isFrom_grandtruth=[]
	
	
	old_estimd_snp_isFrom_grandtruth=False   # current_estimd_snp_isFrom_
	counter_switch=0
	start=True
	for idx_estimated,val_estimated in dic_block_estimated.items():
		if idx_estimated in dic_truth:
			leng_block_intersct_truth+=1
		
			if val_estimated==dic_truth[idx_estimated]:
				current_estimd_snp_isFrom_grandtruth=True 				
			else:
				current_estimd_snp_isFrom_grandtruth=False
				
			if 	 (not start) & (current_estimd_snp_isFrom_grandtruth!=old_estimd_snp_isFrom_grandtruth):
				counter_switch+=1
			
			old_estimd_snp_isFrom_grandtruth=current_estimd_snp_isFrom_grandtruth
			start=False 		
				
				
				
			estimd_snp_isFrom_grandtruth.append(current_estimd_snp_isFrom_grandtruth)
			
			if val_estimated==dic_truth[idx_estimated]:
				correct_num_1+=1
			else:
				correct_num_2+=1
				
	#print(counter_switch)			
	#print(estimd_snp_isFrom_grandtruth)
	crrct=max(correct_num_1,correct_num_2)

	#if leng_block_intersct_truth!=0: 
	#	rr=float(crrct)/leng_block_intersct_truth
	#	counter_switch=float(counter_switch)/leng_block_intersct_truth

#	else:
#		rr=0  ###  ???

	return leng_block_intersct_truth, counter_switch, crrct


if __name__ == "__main__":
	Name_hap_truth=argv[1]
	Name_hap_matlab=argv[2]
	dic_hap_truth=parse_hap_truth(Name_hap_truth)
	dic_hap_estimated=parse_hap_estimated_blocks(Name_hap_matlab)
	
	leng_block_all=[]
	crrct_all=[]
	counter_switch_all=[]
	
	for block_num, dic_block in dic_hap_estimated.items():#enumerate
		leng_block_intersct_truth, counter_switch, crrct=compare(dic_hap_truth,dic_block)
		leng_block_all.append(leng_block_intersct_truth)
		counter_switch_all.append(counter_switch)
		crrct_all.append(crrct)

	
	rr_all=[float(x)/y for x, y in zip(crrct_all, leng_block_all) if y!=0]
	
#	print(np.round(rr_all,2))
	#print('max is ', max(rr_all))
	print('number of block estimated',len(dic_hap_estimated) )
#	print('length is')
#	print(len(crrct_all))

	mean_rr=np.mean(rr_all)
	print('rr mean is ',np.round(mean_rr,4))
	rr_weighted=[x*y for x,y in zip(rr_all,leng_block_all) ]
	mean_rr_weighted=np.mean(rr_weighted)
	print('rr weighted mean is ',np.round(mean_rr_weighted,3))

	
	swer_all=[float(x)/y for x, y in zip(counter_switch_all, leng_block_all) if y!=0]



	print('sum of block lengh is',sum(leng_block_all) )
	print('switch error rate',np.mean(swer_all))
#	print('dic truth leng is',len(dic_hap_truth))
	print('coverage of genome',float(sum(leng_block_all))/len(dic_hap_truth))
	
	print('switch error ',float(sum(counter_switch_all)/sum(leng_block_all) ))

	
	
	
	