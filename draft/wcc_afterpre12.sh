#!/bin/bash


cc_extr='/mnt/LTR_userdata/majid001/software/sdhap/hap_poly4'
frag_weak='/mnt/LTR_userdata/majid001/software/code_py/frag_weak/frag_v3c2.py'
sdhap='/mnt/LTR_userdata/majid001/software/sdhap/hap_poly'

bamprocess='/mnt/LTR_userdata/majid001/software/bamprocess_chaning/bamprocess4e.py'
fragpoly='/mnt/LTR_userdata/majid001/software/code_py/FragmentPoly.py' # python2
sdhap='/mnt/LTR_userdata/majid001/software/sdhap/hap_poly'
split_molec='/mnt/LTR_userdata/majid001/software/code_py/split_v3b.py'
#ConvertAllelesSDhaP='/mnt/LTR_userdata/majid001/software/code_py/ConvertAllelesSDhaP.py' # python2
#hapcompare='/mnt/LTR_userdata/majid001/software/code_py/compare/hapcompare_v1.py' # python2



k=$2 # ploidy level



cd $1  &&
mkdir wcc;cd wcc &&

cc_num=$(($(find ./ -maxdepth 1 -type d| wc -l)-1))

for (( i=(($cc_num-1)); i>=0; i--)); do  
    cd ${i} &&
    part_num_i=$(ls | grep "_sd.txt" | wc -l)
    for ((j=1; j<=$part_num_i; j++)); do
        if [ -f frag${i}_${j}_sd.txt ]; then
        $sdhap frag${i}_${j}_sd.txt  raw${i}_${j}.hap $k  # &&
        fi

        # stats
        #python2 $ConvertAllelesSDhaP -p raw${i}_${j}.hap -o p${i}_${j}.hap -v ../../frb/var_het.vcf  &&
        #python2 $hapcompare  ../../haplogen/hap.variants_het.txt  p${i}_${j}.hap -t -v > res_${i}_${j}.txt &&
        #grep "Allelic correlation" res_${i}_${j}.txt | head -n 1 | sed 's/Allelic correlation between the corresponding blocks (including variants with discordant dosages):/    /g' >> ../results_rr.txt &&
        #grep "Vector Error Rates" res_${i}_${j}.txt | head -n 1 | sed 's/Vector Error Rates for common variants of each block pair: /    /g' >> ../results_ver.txt &&
        #cat p${i}_${j}.hap >> ../hap.txt
        
        
    done
    cd ..
done



wait
pwd


#cat results_rr.txt | sed 's/    /\t/g' | cut -f 2 |sed 's/ ,/\n/g' | tr '\n' ,  > results_rr_pure.txt
#cat results_ver.txt| sed 's/    /\t/g' | cut -f 2 |sed 's/ ,/\n/g' | tr '\n' ,  > results_ver_pure.txt
#grep -n "lock" hap.txt | cut -f1 -d: | tr '\n' ,  >line_number
#cat hap.txt| wc -l >>line_number

#cat results_rr_pure.txt
#cat line_number
#cat results_ver_pure.txt
#cat results_rr_pure.txt | awk '{ sum+=$1 } END { print "Average:",sum/NR,"\nNumber of block:",NR}'
#wc -l hap.hap
#grep -v "Block" hap.hap | sort -u  > uniq.hap
#wc -l uniq.hap
