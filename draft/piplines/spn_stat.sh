#!/bin/bash


ConvertAllelesSDhaP='/mnt/LTR_userdata/majid001/software/code_py/ConvertAllelesSDhaP.py' # python2
hapcompare='/mnt/LTR_userdata/majid001/software/code_py/compare/hapcompare_v1.py' # python2


#k=$2 # ploidy level



cc_num=$(($(find ./ -maxdepth 1 -type d| wc -l)-1))
for (( i=(($cc_num)); i>0; i--)); do  
    cat ${i}/raw.hap >> raw.txt
done

mkdir hap_filtering; cd hap_filtering
csplit --quiet -z ../raw.txt -f hap -n 1 /Block/ {*}
num=$(ls | grep "hap" | wc -l)
echo -n '' > raw_filtered.hap
for ((j=0; j<$(($num-2)); j++)); do
    len_j=$(cat hap$j | wc -l )
    if [ $len_j -gt 5 ]; then 
    cat hap$j >> raw_filtered.hap
    fi
done

python2 $ConvertAllelesSDhaP -p raw_filtered.hap -o p.hap -v ../../frb/var_het.vcf
cd ..
cp hap_filtering/p.hap .
python2 $hapcompare  ../haplogen/hap.variants_het.txt  p.hap -t -v > res.txt &&


grep "Allelic correlation" res.txt | head -n 1 | sed 's/Allelic correlation between the corresponding blocks (including variants with discordant dosages):/    /g' > results_rr_pure.txt
grep "Vector Error Rates" res.txt | head -n 1 | sed 's/Vector Error Rates for common variants of each block pair: /    /g' > results_ver_pure.txt
grep -n "lock" hap.txt | cut -f1 -d: | tr '\n' ,  >line_number

cat hap.txt| wc -l >>line_number

#echo -n ' ' > results_rr.txt
#echo -n ' ' > results_ver.txt
#echo -n ' ' > hap.txt




#cc_num=$(($(find ./ -maxdepth 1 -type d| wc -l)-1))
#for (( i=(($cc_num)); i>0; i--)); do  
#    cd ${i}
#        (
        #python2 $ConvertAllelesSDhaP -p raw.hap -o p.hap -v ../../frb/var_het.vcf  &&
        #python2 $hapcompare  ../../haplogen/hap.variants_het.txt  p.hap -t -v > res.txt &&

        #grep "Allelic correlation" res.txt | head -n 1 | sed 's/Allelic correlation between the corresponding blocks (including variants with discordant dosages):/    /g' >> ../results_rr.txt &&
        #grep "Vector Error Rates" res.txt | head -n 1 | sed 's/Vector Error Rates for common variants of each block pair: /    /g' >> ../results_ver.txt &&
        #cat p.hap >> ../hap.txt
#        )&
#    cd ..
#done

wait
pwd


#cat results_rr.txt | sed 's/    /\t/g' | cut -f 2 |sed 's/ ,/\n/g' | tr '\n' ,  > results_rr_pure.txt
#echo '' >> results_rr_pure.txt
#cat results_ver.txt| sed 's/    /\t/g' | cut -f 2 |sed 's/ ,/\n/g' | tr '\n' ,  > results_ver_pure.txt
#echo '' >> results_ver_pure.txt
#grep -n "lock" hap.txt | cut -f1 -d: | tr '\n' ,  >line_number
#cat hap.txt| wc -l >>line_number

#cat results_rr_pure.txt
#cat line_number
#cat results_ver_pure.txt
#cat results_rr_pure.txt | awk '{ sum+=$1 } END { print "Average:",sum/NR,"\nNumber of block:",NR}'
#wc -l hap.hap
#grep -v "Block" hap.hap | sort -u  > uniq.hap
#wc -l uniq.hap
