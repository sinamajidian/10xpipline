#!/bin/bash


cc_extr='/mnt/LTR_userdata/majid001/software/sdhap/hap_poly4'
frag_weak='/mnt/LTR_userdata/majid001/software/code_py/frag_weak/frag_v3c2.py'
sdhap='/mnt/LTR_userdata/majid001/software/sdhap/hap_poly'

bamprocess='/mnt/LTR_userdata/majid001/software/bamprocess_chaning/bamprocess4b.py'
fragpoly='/mnt/LTR_userdata/majid001/software/code_py/FragmentPoly.py' # python2
sdhap='/mnt/LTR_userdata/majid001/software/sdhap/hap_poly'
split_molec='/mnt/LTR_userdata/majid001/software/code_py/old/split_v2c.py'
ConvertAllelesSDhaP='/mnt/LTR_userdata/majid001/software/code_py/ConvertAllelesSDhaP.py' # python2
hapcompare='/mnt/LTR_userdata/majid001/software/code_py/compare/hapcompare_v1.py' # python2



k=3 # ploidy level



mkdir weakly_v2
cd weakly_v2
cp ../split/frag_sp.txt .

python2 $fragpoly -f frag_sp.txt  -o frag_sd.txt -x SDhaP

$cc_extr frag_sd.txt out_test $k # extracting connected componnets (part of sdhap code )
csplit -z connected_dic.txt -f cc -n 1 /block/ {*} 
cc_num=$(ls | grep "cc" | wc -l)

for (( i=(($cc_num-1)); i>=0; i--)); do  
    mkdir ${i}; cd ${i}
    mv ../cc${i} .
    grep -v "block" cc${i} >  cc${i}_pure
    awk 'NR == FNR{a[$0]; next};FNR in a' cc${i}_pure ../frag_sp.txt > frag${i}.txt
    python2 $fragpoly -f frag${i}.txt -o frag${i}_dic.txt -x HapTree
    
    #extracting weak parts 
    python3 $frag_weak frag${i}_dic.txt 
    part_num_i=$(ls | grep "part" | wc -l)
    for ((j=1; j<=part_num_i; j++)); do
        (
        awk 'NR == FNR{a[$0]; next};FNR in a' part${j}.txt frag${i}.txt > frag${i}_${j}.txt &&
        python2 $fragpoly -f frag${i}_${j}.txt -o frag${i}_${j}_sd.txt -x SDhaP &&
        $sdhap frag${i}_${j}_sd.txt  raw${i}_${j}.hap $k  &&
        
        
        # stats
        python2 $ConvertAllelesSDhaP -p raw${i}_${j}.hap -o p${i}_${j}.hap -v ../../frb/var_het.vcf  &&
        python2 $hapcompare  ../../haplogen/hap.variants_het.txt  p${i}_${j}.hap -t -v > res_${i}_${j}.txt &&
        grep "Allelic correlation" res_${i}_${j}.txt | head -n 1 | sed 's/Allelic correlation between the corresponding blocks (including variants with discordant dosages):/    /g' >> ../results_rr.txt &&
        grep "Vector Error Rates" res_${i}_${j}.txt | head -n 1 | sed 's/Vector Error Rates for common variants of each block pair: /    /g' >> ../results_ver.txt &&
        cat p${i}_${j}.hap >> ../hap.txt
        
        )&
    done
    cd ..
done



wait
pwd


cat results_rr.txt | sed 's/    /\t/g' | cut -f 2 |sed 's/ ,/\n/g' | tr '\n' ,  > results_rr_pure.txt
cat results_ver.txt| sed 's/    /\t/g' | cut -f 2 |sed 's/ ,/\n/g' | tr '\n' ,  > results_ver_pure.txt
grep -n "lock" hap.txt | cut -f1 -d: | tr '\n' ,  >line_number
cat hap.txt| wc -l >>line_number

#cat results_rr_pure.txt
#cat line_number
#cat results_ver_pure.txt
#cat results_rr_pure.txt | awk '{ sum+=$1 } END { print "Average:",sum/NR,"\nNumber of block:",NR}'
#wc -l hap.hap
#grep -v "Block" hap.hap | sort -u  > uniq.hap
#wc -l uniq.hap
