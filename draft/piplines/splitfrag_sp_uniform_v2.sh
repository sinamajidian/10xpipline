#!/bin/bash

bamprocess='/mnt/LTR_userdata/majid001/software/bamprocess_chaning/bamprocess4b.py'
fragpoly='/mnt/LTR_userdata/majid001/software/code_py/FragmentPoly.py' # python2
sdhap='/mnt/LTR_userdata/majid001/software/sdhap/hap_poly'
split_molec='/mnt/LTR_userdata/majid001/software/code_py/old/split_v2c.py'
ConvertAllelesSDhaP='/mnt/LTR_userdata/majid001/software/code_py/ConvertAllelesSDhaP.py' # python2
hapcompare='/mnt/LTR_userdata/majid001/software/code_py/compare/hapcompare_v1.py' # python2

k=3 # ploidy level
mkdir spf_sp2
cd spf_sp2

awk -vc=1 'NR%2500==0{++c}{print $0 > c".txt"}'  ../split/frag_sp.txt   &&  
part_num=$(ls | wc -l)

for (( i=(($part_num-1)); i>=1; i--)); do  
    (
    mkdir $i; cd $i
    mv ../${i}.txt . &&
    python2 $fragpoly -f ${i}.txt -o frag_sd.txt -x SDhaP  &&
    $sdhap frag_sd.txt raw.hap $k &&
    
    
    ## stats
    python2 $ConvertAllelesSDhaP -p raw.hap -o p.hap -v  ../../frb/var_het.vcf  &&
    python2 $hapcompare  ../../haplogen/hap.variants_het.txt p.hap -t -v > res.txt &&
    
    grep "Allelic correlation" res.txt | head -n 1 | sed 's/Allelic correlation between the corresponding blocks (including variants with discordant dosages):/    /g' >> ../results_rr.txt  &&
    grep "Vector Error Rates" res.txt | head -n 1 | sed 's/Vector Error Rates for common variants of each block pair: /    /g' >> ../results_ver.txt  &&
    cat p.hap   >> ../hap.txt  
    
    )&
done

wait
pwd

cat results_rr.txt | tr '\n' , > results_rr_p.txt
cat results_ver.txt | tr '\n' , > results_ver_p.txt
sed -i 's/inf/-1/g' results_ver_p.txt
grep -n "lock" hap.txt | cut -f1 -d: | tr '\n' ,  > line_number
cat hap.txt| wc -l >> line_number

