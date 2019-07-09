#!/bin/bash


#cc_extr='/mnt/LTR_userdata/majid001/software/sdhap/hap_poly4'
#frag_weak='/mnt/LTR_userdata/majid001/software/code_py/frag_weak/frag_v3c2.py'
#sdhap='/mnt/LTR_userdata/majid001/software/sdhap/hap_poly'

#bamprocess='/mnt/LTR_userdata/majid001/software/bamprocess_chaning/bamprocess4e.py'
fragpoly='/mnt/LTR_userdata/majid001/software/code_py/FragmentPoly.py' # python2
sdhap='/mnt/LTR_userdata/majid001/software/sdhap/hap_poly'
split_molec='/mnt/LTR_userdata/majid001/software/code_py/split_v3b.py'
#ConvertAllelesSDhaP='/mnt/LTR_userdata/majid001/software/code_py/ConvertAllelesSDhaP.py' # python2
#hapcompare='/mnt/LTR_userdata/majid001/software/code_py/compare/hapcompare_v1.py' # python2


#time( 
cd $1  &&
k=$2   &&
cd frb   &&
python3 $split_molec frag.txt pos_freebayes.txt 50   &&
cd ..   
mkdir sdhap;  cp frb/frag_sp.txt sdhap; cd sdhap 
python2 $fragpoly -f frag_sp.txt  -o frag_sd.txt -x SDhaP   &&
$sdhap frag_sd.txt  out_sd_raw.hap  $k   &&
#)


#wait
pwd

