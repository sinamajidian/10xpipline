#!/bin/bash

#10x pipline
#triploid

#hap_generator_path='/mnt/scratch/majid001/software/haplogenerator/haplogenerator.py'


genomeGenerating()
{    
    mkdir ref
    # soltub_200k.fasta  soltub_1M_2_up.fasta
    #cp /mnt/LTR_userdata/majid001/genomes/solanum_tuberosum/solanum.fa ref/ref.fasta
    #head -n 1000 ref/ref1.fasta > ref/ref.fasta
    #samtools faidx ref/ref.fasta
    #longranger mkref  ref/ref.fasta
    
    #printf '=%.0s' {1..100}; printf "\n"; printf "haplogenerator started\n"; printf '=%.0s' {1..100}; printf "\n"
    #mkdir haplogen
    #python2 $hap_generator_path -f ref/ref.fasta -o haplogen/sim -m "{'A':'C','C':'A','G':'T','T':'G'}"  -p 3 -v --model poisson -s [.01,0,0] 
    #--dosage "[0.5,0.23,0.27]"
    #cd haplogen
    #grep -P "\tN\t" sim_varianthaplos.txt | wc -l;
    #grep -P "\t2\t" sim_varianthaplos.txt  | wc -l
    #grep -vP "1\t1\t1" sim_varianthaplos.txt | grep -vP "0\t0\t0" > hap.variants_het.txt # for triploid !!
    #cut -f 3  hap.variants_het.txt | grep -v "contig" > pos_truth.txt
    #mkdir edited; cp *.fa edited; cd edited
    #for idx in {1..3}; do sed -i "s/^>/>hap${idx}_/g" sim_hap$idx.fa; done # there is just one part in ref fasta
    #for idx in {1..3}; do cat sim_hap${idx}.fa >> genome.fa; done
    #cd ../..
}
    
ReadGenerating()
{
    printf '=%.0s' {1..100}; printf "\n"; printf "LRSIM started\n"; printf '=%.0s' {1..100}; printf "\n"
    mkdir out_lrsim
    #simulateLinkedReads_notRandom-m.pl  # number of molecules is not random    -f change , also change the spliter
    
    
    ~/Install/LRSIM/simulateLinkedReads_notRandom-m.pl -g haplogen/edited/genome.fa -p out_lrsim/sim \
     -o -x 9.3 -f 20 -m 1 -t 66.5
    
    
    
    cd out_lrsim; mkdir edited
    cp sim_S1_L00* edited/
    cd edited
    gunzip *
    sed 's/\/1/ /g' sim_S1_L001_R1_001.fastq > sim_ed_S1_L001_R1_001.fastq;
    sed 's/\/2/ /g' sim_S1_L001_R2_001.fastq > sim_ed_S1_L001_R2_001.fastq;
    gzip -k sim_ed_S1_L001_R1_001.fastq;
    gzip -k sim_ed_S1_L001_R2_001.fastq
    mkdir filegz
    mv *.gz filegz
    cd ../../
    printf '=%.0s' {1..100}; printf "\n"; printf "LRSIM finished\n"; printf '=%.0s' {1..100}; printf "\n"    

}

Aligning()
{
    printf '=%.0s' {1..100}; printf "\n"; printf "LongRanger started\n"; printf '=%.0s' {1..100}; printf "\n"    
    cd ./$2
    longranger align --id=out_longranger --fastqs=out_lrsim/edited/filegz/ --reference=refdata-ref
    printf '=%.0s' {1..100}; printf "\n"; printf "LongRanger finished\n"; printf '=%.0s' {1..100}; printf "\n"    

}

VariantCalling()
{
    printf '=%.0s' {1..100}; printf "\n"; printf "Variant calling  started\n"; printf '=%.0s' {1..100}; printf "\n"
    mkdir frb
    cd frb
    cp ../out_longranger/outs/possorted_bam.bam .
    samtools index possorted_bam.bam
    cp ../haplogen/pos_truth.txt .
    freebayes -f ../ref/ref.fasta -p 3  possorted_bam.bam  > var.vcf
    /mnt/scratch/majid001/software/break_vcf.sh var.vcf var_break.vcf
    cat var_break.vcf  | grep -v "1/1/1" | grep -v "0/0/0"  > var_het.vcf
    grep -v "#" var_het.vcf | cut -f 2  > pos_freebayes.txt
    cd ..
}


HaplotypingSdhap()
{
    cd frb
    python /mnt/scratch/majid001/software/a/bamprocess_10x_orig/bamprocess.py possorted_bam.bam var_het.vcf
    mv HapCUT_frags_out_longranger.matrix frag.txt
    cd ..;
    mkdir haplotyping/sdhap;  cp frb/frag.txt haplotyping/sdhap;  cd haplotyping/sdhap;
    #python2 /mnt/scratch/majid001/software/code_py/FragmentPoly.py -f frag.txt  -o frag_sd.txt -x SDhaP
    #/mnt/scratch/majid001/software/sdhap/hap_poly frag_sd.txt  out_sd_raw.hap  3

    #python2 /mnt/scratch/majid001/software/code_py/ConvertAllelesSDhaP.py -p out_sd_raw.hap -o out_sdhap.hap -v ../../frb/var_het.vcf
    #cut -f 2 out_sdhap.hap > pos_sdhap.txt 
    #sed -i '1d' pos_sdhap.txt 
    #cd ../..
    #python2 /mnt/scratch/majid001/software/code_py/compare/hapcompare_v1.py haplogen/hap.variants_het.txt haplotyping/sdhap/out_sdhap.hap   -t -v   > haplotyping/res_sdhap.txt
}





######  main
#usage:  ./code.sh 31

#mkdir $1
#cd $1
cd  65c_m10

#genomeGenerating
ReadGenerating
Aligning
VariantCalling

mkdir haplotyping
HaplotypingSdhap  



#HaplotypingSdhapSplited



