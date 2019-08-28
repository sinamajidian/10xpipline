#!/bin/bash

# split based on N region in reference genome

#usage:  .
ConvertAllelesSDhaP='/mnt/LTR_userdata/majid001/software/code_py/ConvertAllelesSDhaP.py' # python2

bamprocess='/mnt/LTR_userdata/majid001/software/bamprocess_chaning/bamprocess4e.py' #bamprocess_10x/bamprocess.py' #
fragpoly='/mnt/LTR_userdata/majid001/software/code_py/FragmentPoly.py' # python2
sdhap='/mnt/LTR_userdata/majid001/software/sdhap/hap_poly'
split_molec='/mnt/LTR_userdata/majid001/software/code_py/split_v3b.py'
nonnregion='/mnt/LTR_userdata/majid001/software/code_py/nfinder_v3.py'


k=$3
cd $1

mkdir spn; cd spn
mean_shif_on=$2 # on or off

python3 $nonnregion ../ref/ref.fasta 


num=0
while IFS=$'\t' read -r a b; do  
    num=$(( $num + 1 ))
    mkdir $num
    samtools view -b ../frb/possorted_bam.bam "SolTub_3.0_chr1:${a}-${b}"  -o ${num}/part${num}.bam 
done < nonNregion.bed


p_num=$(cat nonNregion.bed| wc -l)
#for (( i=(($p_num)); i>0; i--)); do  
for (( i=1; i<=$p_num; i++)); do  
    (
    cd $i &&
    python $bamprocess part${i}.bam  ../../frb/var_het.vcf &&
    mv HapCUT_frags_out_longranger.matrix frag.txt &&
    if [ ${mean_shif_on} = 'on' ]; then
        python3 $split_molec  frag.txt ../../frb/pos_freebayes.txt 50 &&
        python2 $fragpoly -f frag_sp.txt -o frag_sd.txt -x SDhaP
    else
        python2 $fragpoly -f frag.txt -o frag_sd.txt -x SDhaP 
    fi
    $sdhap frag_sd.txt  raw.hap $k
    )&
done


wait

echo '' >raw_all.hap
for (( i=1; i<=$p_num; i++)); do    
    cat ${i}/raw.hap >> raw_all.hap
done

mkdir hap_filtering; cd hap_filtering  
csplit --quiet -z ../raw_all.hap -f hap -n 1 /Block/ {*}  &&
num=$(ls | grep "hap" | wc -l)
echo -n '' > raw_filtered.hap
for ((j=0; j<$(($num)); j++)); do
    len_j=$(cat hap$j | wc -l )
    if [ $len_j -gt 2 ]; then 
    cat hap$j >> raw_filtered.hap
    fi
done

python2 $ConvertAllelesSDhaP -p raw_filtered.hap -o p.hap -v ../../frb/var_het.vcf  &&
cd ..
cp hap_filtering/p.hap haplotype.hap 

pwd

