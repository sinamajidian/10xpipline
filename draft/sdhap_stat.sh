



ConvertAllelesSDhaP='/mnt/LTR_userdata/majid001/software/code_py/ConvertAllelesSDhaP.py'
hapcompare='/mnt/LTR_userdata/majid001/software/code_py/compare/hapcompare_v1.py'


cd sdhap  &&
python2 $ConvertAllelesSDhaP -p out_sd_raw.hap -o hap.txt  -v ../frb/var_het.vcf  &&
python2 $hapcompare ../haplogen/hap.variants_het.txt hap.txt    -t -v   > res_sdhap.txt   &&


grep "Allelic correlation" res_sdhap.txt | head -n 1 | sed 's/Allelic correlation between the corresponding blocks (including variants with discordant dosages):/    /g' > results_rr_pure.txt
grep "Vector Error Rates" res_sdhap.txt | head -n 1 | sed 's/Vector Error Rates for common variants of each block pair: /    /g' > results_ver_pure.txt
     

grep -n "lock" hap.txt | cut -f1 -d: | tr '\n' ,  >line_number  &&
cat hap.txt| wc -l >>line_number
