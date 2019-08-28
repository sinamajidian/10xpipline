#!/bin/bash


ConvertAllelesSDhaP='/mnt/LTR_userdata/majid001/software/code_py/ConvertAllelesSDhaP.py' # python2

cc_num=$(($(find ./ -maxdepth 1 -type d| wc -l)-1))

echo -n '' > raw.hap
for (( i=0; i<$cc_num; i++)); do  
    #ls
    cd ${i}
    part_num_i=$(ls | grep ".hap" | grep -v "raw" | wc -l)
    for ((j=1; j<=$part_num_i; j++)); do
        if [ -f frag${i}_${j}.hap ]; then
            cat  frag${i}_${j}.hap >> ../raw.hap
        fi
    done
    cd ..
done

python2 $ConvertAllelesSDhaP -p raw.hap -o haplotype.hap -v ../frb/var_het.vcf 
pwd
