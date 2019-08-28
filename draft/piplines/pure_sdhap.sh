#!/bin/bash

fragpoly='/mnt/LTR_userdata/majid001/software/code_py/FragmentPoly.py' # python2
sdhap='/mnt/LTR_userdata/majid001/software/sdhap/hap_poly'

ConvertAllelesSDhaP='/mnt/LTR_userdata/majid001/software/code_py/ConvertAllelesSDhaP.py' # python2



cd $1  &&
k=$2   &&
   
mkdir pure_sdhap;  cp frb/frag.txt pure_sdhap/; cd pure_sdhap 
python2 $fragpoly -f frag.txt  -o frag_sd.txt -x SDhaP   &&
$sdhap frag_sd.txt  out_sd_raw.hap  $k   &&
python2 $ConvertAllelesSDhaP -p out_sd_raw.hap -o haplotype.hap -v ../frb/var_het.vcf 


pwd

