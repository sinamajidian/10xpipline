#!/bin/bash

#10x pipline /32/60



#triploid


genomeGenerating()
{    
    mkdir ref
    # soltub_200k.fasta  soltub_1M_2_up.fasta
    cp data/soltub_200k.fasta ref/
    cd ref
    cd ..
    samtools faidx ref/ref.fasta
    longranger mkref  ref/ref.fasta
    
    printf '=%.0s' {1..100}; printf "\n"; printf "haplogenerator started\n"; printf '=%.0s' {1..100}; printf "\n"
    mkdir haplogen
    python /mnt/nexenta/majid001/nobackup/software/haplogenerator/haplogenerator_v1.py -f ref/ref.fasta -o haplogen/sim \
    -m "{'A':'C','C':'A','G':'T','T':'G'}" --dosage "[0.5,0.23,0.27]" -p 3 -v --model poisson -s [.01,0,0]
    cd haplogen
    grep -P "\tN\t" sim_varianthaplos.txt | wc -l;
    grep -P "\t2\t" sim_varianthaplos.txt  | wc -l
    grep -vP "1\t1\t1" sim_varianthaplos.txt | grep -vP "0\t0\t0" > hap.variants_het.txt # for triploid !!
    cut -f 3  hap.variants_het.txt | grep -v "contig" > pos_truth.txt
    mkdir edited; cp *.fa edited; cd edited
    for idx in {1..3}; do sed -i "s/NW_006238927.1/NW_006238927\.1_$idx /g" sim_hap$idx.fa; done 
    for idx in {1..3}; do cat sim_hap${idx}.fa >> genome.fa; done
    cd ../..
}
    
ReadGenerating()
{
    printf '=%.0s' {1..100}; printf "\n"; printf "LRSIM started\n"; printf '=%.0s' {1..100}; printf "\n"
    mkdir out_lrsim
    #simulateLinkedReads_notRandom-m.pl  # number of molecules is not random    -f change , also change the spliter
    
    
    ~/Install/LRSIM/simulateLinkedReads_notRandom-m.pl -g haplogen/edited/genome.fa -p out_lrsim/sim \
     -o -x .15 -f 50 -m 1 -t 4.5 -b ../barcode15k.txt
    
    
    
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
    cp ../haplogen/pos_truth.txt .
    freebayes -f ../ref/ref.fasta -p 3  possorted_bam.bam  > var.vcf
    /mnt/nexenta/majid001/nobackup/software/break_vcf.sh var.vcf var_break.vcf
    cat var_break.vcf  | grep -v "1/1/1" | grep -v "0/0/0"  > var_het.vcf
    grep -v "#" var_het.vcf | cut -f 2  > pos_freebayes.txt
    cd ..
}


HaplotypingHapcompass()
{
    mkdir haplotyping/hapcom
    java -Xmx8g -jar /mnt/nexenta/majid001/nobackup/software/hapcompass.jar --ploidy 3 --bam  frb/possorted_bam.bam --vcf frb/var_het.vcf -o haplotyping/hapcom/a
    python /mnt/nexenta/majid001/nobackup/software/hapcompare_v1.py haplogen/hap.variants_het.txt   haplotyping/hapcom/a_MWER_solution.txt   -t -v   > haplotyping/res_hapcomp.txt
}

######  main
#usage:  ./code.sh 31

mkdir $1
cd $1
genomeGenerating
ReadGenerating
Aligning
VariantCalling

mkdir haplotyping
HaplotypingHapcompass 
