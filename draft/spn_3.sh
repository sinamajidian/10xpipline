#!/bin/bash

# split based on N region in reference genome

#usage:  ./sbu.sh out_folder 

sdhap='/mnt/LTR_userdata/majid001/software/sdhap/hap_poly'
bamprocess='/mnt/LTR_userdata/majid001/software/bamprocess_chaning/bamprocess4e.py' #bamprocess_10x/bamprocess.py' #
fragpoly='/mnt/LTR_userdata/majid001/software/code_py/FragmentPoly.py' # python2
sdhap='/mnt/LTR_userdata/majid001/software/sdhap/hap_poly'
split_molec='/mnt/LTR_userdata/majid001/software/code_py/split_v3b.py'
#ConvertAllelesSDhaP='/mnt/LTR_userdata/majid001/software/code_py/ConvertAllelesSDhaP.py' # python2
#hapcompare='/mnt/LTR_userdata/majid001/software/code_py/compare/hapcompare_v1.py' # python2
#lrsim='/home/majid001/Install/LRSIM/simulateLinkedReads_notRandom-m.pl'
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

for (( i=(($p_num)); i>0; i--)); do  
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
    
    ## stats
    #python2 $ConvertAllelesSDhaP -p raw.hap -o p.hap -v  ../../frb/var_het.vcf  &&
    #python2 $hapcompare  ../../haplogen/hap.variants_het.txt p.hap -t -v > res.txt &&
    
    #grep "Allelic correlation" res.txt | head -n 1 | sed 's/Allelic correlation between the corresponding blocks (including variants with discordant dosages):/    /g' >> ../results_rr.txt  &&
    #grep "Vector Error Rates" res.txt | head -n 1 | sed 's/Vector Error Rates for common variants of each block pair: /    /g' >> ../results_ver.txt  &&
    #cat p.hap   >> ../hap.txt  
    
    )&
done




wait
pwd


#cat results_rr.txt | sed 's/    /\t/g' | cut -f 2 |sed 's/ ,/\n/g' | tr '\n' ,  > results_rr_pure.txt
#cat results_ver.txt| sed 's/    /\t/g' | cut -f 2 |sed 's/ ,/\n/g' | tr '\n' ,  > results_ver_pure.txt
#grep -n "lock" hap.txt | cut -f1 -d: | tr '\n' ,  >line_number
#cat hap.txt| wc -l >>line_number










