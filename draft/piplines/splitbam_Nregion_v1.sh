#!/bin/bash

bamprocess='/mnt/scratch/majid001/software/a/bamprocess_10x_orig/bamprocess.py'
fragpoly='/mnt/scratch/majid001/software/code_py/FragmentPoly.py'
sdhap='/mnt/scratch/majid001/software/sdhap/hap_poly'
split_molec='/mnt/scratch/majid001/software/code_py/old/split_v2c.py'

ConvertAllelesSDhaP='/mnt/scratch/majid001/software/code_py/ConvertAllelesSDhaP.py'
hapcompare='/mnt/scratch/majid001/software/code_py/compare/hapcompare_v1.py'





#python3 /mnt/scratch/majid001/software/code_py/Nfinder_v2.py ref.fasta > N_ref.bed
#head -n  5 N_ref.bed >N_ref_head5.bed

num=0
while IFS="    " read -r a b; do  
    num=$(( $num + 1 ))
    mkdir $num
    samtools view -b ../frb/possorted_bam.bam "SolTub_3.0_chr1:${a}-${b}"  -o ${num}/part${num}.bam &
done < N_ref.bed


for idx in {71..123}; do
cd $idx
python /mnt/scratch/majid001/software/a/bamprocess_10x_orig/bamprocess.py part${idx}.bam   ../../frb/var_het.vcf &
cd ..
done



for idx in {55..123}; do 
cd $idx
mv HapCUT_frags_out_longranger.matrix frag.txt
python2 /mnt/scratch/majid001/software/code_py/FragmentPoly.py -f frag.txt  -o frag_sd.txt -x SDhaP
/mnt/scratch/majid001/software/sdhap/hap_poly frag_sd.txt  out_sd_raw.hap  3 &
cd ..
done



for idx in {55..123}; do
cd $idx
python2 /mnt/scratch/majid001/software/code_py/ConvertAllelesSDhaP.py -p out_sd_raw.hap -o out_sdhap.hap -v ../../frb/var_het.vcf
python2 /mnt/scratch/majid001/software/code_py/compare/hapcompare_v1.py  ../../haplogen/hap.variants_het.txt  out_sdhap.hap   -t -v   > res_sdhap.txt &
cd ..
done







 
for i in {1..123}; do
grep "Vector Error Rates" $i/res_sdhap.txt  | head -n 1 >> results_ver.txt
done


 
for i in {1..123}; do
grep "Allelic correlatio" $i/res_sdhap.txt  | head -n 1 >> results_rr.txt
done

sed 's/Vector Error Rates for common variants of each block pair:/    /g' results_ver.txt > results_ver_pure.txt

 
for i in {1..123}; do
cat $i/out_sdhap.hap  >> hap.txt
done








 

python2 /mnt/scratch/majid001/software/code_py/FragmentPoly.py -f frag.txt  -o frag_dic.txt -x HapTree

python2 /mnt/scratch/majid001/software/code_py/FragmentPoly.py -f frag3.txt  -o frag_dic.txt -x HapTree

python3 /mnt/scratch/majid001/software/code_py/ a.txt

#Find N in refere


python3 /mnt/scratch/majid001/software/code_py/Nfinder.py ref.fasta 
python /mnt/scratch/majid001/software/a/bamprocess_10x_orig/bamprocess.py part2.bam   ../../../var_het.vcf 

mv HapCUT_frags_out_longranger.matrix frag.txt

python2 /mnt/scratch/majid001/software/code_py/FragmentPoly.py -f frag.txt  -o frag_sd.txt -x SDhaP

/mnt/scratch/majid001/software/sdhap/hap_poly frag_sd.txt  out_sd_raw.hap  3

time(/mnt/scratch/majid001/software/sdhap/hap_poly frag_sd.txt  out_sd_raw2.hap  4)




Compare
python2 /mnt/scratch/majid001/software/code_py/ConvertAllelesSDhaP.py -p out_sd_raw.hap -o out_sdhap.hap -v ../../../../../frb/var_het.vcf

python2 /mnt/scratch/majid001/software/code_py/hapcompare_v1.py ../../../../../haplogen/hap.variants_het.txt out_sdhap.hap   -t -v   > res_sdhap.txt

python2 /mnt/scratch/majid001/software/code_py/ConvertAllelesSDhaP.py -p out_sd_raw.hap -o out_sdhap.hap -v ../../../var_het.vcf


python2 /mnt/scratch/majid001/software/code_py/hapcompare_v1.py  ../../../../haplogen/hap.variants_het.txt out_sdhap.hap -t -v   > res_sdhap.txt


loop
for idx in {1..5}; do python2 /mnt/scratch/majid001/software/code_py/FragmentPoly.py -f frag${idx}.txt  -o frag_sd${idx}.txt -x SDhaP; done
m
for idx in {1..5}; do
/mnt/scratch/majid001/software/sdhap/hap_poly frag_sd${idx}.txt  out_sd_raw${idx}.hap  3 &
done

for idx in {1..5}; do
grep "Block" out_sd_raw${idx}.hap 
done

for idx in {1..5}; do
python2 /mnt/scratch/majid001/software/code_py/ConvertAllelesSDhaP.py -p out_sd_raw${idx}.hap -o out_sdhap${idx}.hap -v var_het.vcf
done


for idx in {1..5}; do
python2 /mnt/scratch/majid001/software/code_py/hapcompare_v1.py ../../../haplogen/hap.variants_het.txt  out_sdhap${idx}.hap   -t -v   > res_sdhap${idx}.txt &
done



for idx in {1..123}; do
#tail -n 1 ${idx}/res_sdhap.txt >> length.txt
grep "Allelic correlation between the corresponding blocks (including variants with discordant dosages)" ${idx}/res_sdhap.txt >> rr.txt
done


for idx in {1..123}; do
grep "Block" ${idx}/out_sd_raw.hap >> block_number.txt
done



for idx in {21..50}; do
/mnt/scratch/majid001/software/sdhap/hap_poly ${idx}/frag_sd.txt  ${idx}/out_sd_raw_new.hap  3 & >>time.txt
done









k=3 # ploidy level
mkdir split_bam_uniform
cd split_bam_uniform
step=300000
for i in {0..2}; do
    start_=$(($(($i*$step))+1))
    end_=$(($(($i+1))*$step))
    mkdir $i; cd $i
    samtools view -b ../frb/possorted_bam.bam  SolTub_3.0_chr1:${start_}-${end_} > par${i}.bam
    python $bamprocess par${i}.bam  ${addr}var_het.vcf &
    cd ..
done

for i in {0..2}; do
    cd ${i} 
    mv HapCUT_frags_out_longranger.matrix frag.txt
    python3 split_molec  frag.txt ../../frb/pos_freebayes.txt 50
    python2 $fragpoly -f frag_sp.txt -o frag_sd.txt -x SDhaP
    $sdhap frag_sd.txt raw.hap $k &
    cd ..
done

for i in {0..2}; do
    cd ${i} 
    python2 $ConvertAllelesSDhaP -p raw.hap -o p.hap -v  ../../frb/var_het.vcf 
    python2 $hapcompare  ../../haplogen/hap.variants_het.txt p.hap -t -v > res.txt &
    cd ..
done

for i in {0..2}; do
    grep "Allelic correlation" ${i}/res.txt | head -n 1 | sed 's/Allelic correlation between the corresponding blocks (including variants with discordant dosages):/    /g' >> results_rr.txt
    grep "Vector Error Rates" ${i}/res.txt | head -n 1 | sed 's/Vector Error Rates for common variants of each block pair: /    /g' >> results_ver.txt
    cat ${i}/p.hap   >> hap.txt
done
cat results_rr.txt | tr '\n' , > results_rr_p.txt
cat results_ver.txt | tr '\n' , > results_ver_p.txt
sed -i 's/inf/-1/g' results_ver_p.txt
grep -n "lock" hap.txt | cut -f1 -d: | tr '\n' ,  >line_number
cat hap.txt| wc -l >>line_number


