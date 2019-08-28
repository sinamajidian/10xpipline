#!/bin/bash


cc_extr='/mnt/LTR_userdata/majid001/software/sdhap/hap_poly4'
frag_weak='/mnt/LTR_userdata/majid001/software/code_py/frag_weak/frag_v3c2.py'
sdhap='/mnt/LTR_userdata/majid001/software/sdhap/hap_poly'
bamprocess='/mnt/LTR_userdata/majid001/software/bamprocess_chaning/bamprocess4e.py'
fragpoly='/mnt/LTR_userdata/majid001/software/code_py/FragmentPoly.py' # python2
sdhap='/mnt/LTR_userdata/majid001/software/sdhap/hap_poly'
split_molec='/mnt/LTR_userdata/majid001/software/code_py/split_v3b.py'

origin=$(pwd)/$1 # input folder
mkdir $2; cd $2 # output folder
pr=$3  # downsaple rate
k=3

mkdir haplogen
cp $origin/haplogen/hap.variants_het.txt haplogen/ &&
cp -r  $origin/ref/  . &&
mkdir frb; cd frb &&
cp $origin/frb/var_het.vcf . &&
cp $origin/frb/pos_freebayes.txt . &&
cp $origin/frb/possorted_bam.bam  . &&

java -jar ~/Install/picard/picard.jar  DownsampleSam  I=possorted_bam.bam O=file PROBABILITY=$pr  &&
 
mv possorted_bam.bam full  &&
mv file possorted_bam.bam   &&
samtools index possorted_bam.bam  &&
samtools depth  possorted_bam.bam  |  awk '{sum+=$3} END { print "Average of coverage is ",sum/NR}'  &&

python $bamprocess possorted_bam.bam var_het.vcf  &&
mv HapCUT_frags_out_longranger.matrix frag.txt


pwd

