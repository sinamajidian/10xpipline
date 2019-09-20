




from sys import argv

def read_haplo(haplos_file_address):
    
    """ reading the output of haplogenerator
    
    input name of file 
    output list
    """
    haplos_file = open(haplos_file_address,'r'); 
    var_list=[]
    for line in haplos_file:
        line_strip=line.strip() 
        if line_strip.startswith('B'):
            header_line=line_strip
        else:
            line_split=line_strip.split('\t')
            chrom= line_split[1]
            genomic_position=int(line_split[2])
            ref_allele=line_split[3]
            alt_allele=line_split[4]
            hap_values=line_split[5:]

            var_list.append([chrom,genomic_position,ref_allele,alt_allele,hap_values])

    return var_list


def writ_vcf(vcf_file_address,out_unphased):

    
    """ writing a vcf file 
    
    input name of output file 
    output 1
    """
    
    out_unphased=True # if you need unphased,
    vcf_file = open(vcf_file_address,'w'); 

    vcf_file.write('##fileformat=VCFv4.2\n##source=Haplogenerator\n')
    vcf_file.write('##INFO=<ID=NS,Number=1,Type=Integer,Description="test">\n')
    vcf_file.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    vcf_file.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n')


    for var in var_list:
        [chrom,genomic_position,ref_allele,alt_allele,hap_values]=var
        ID='.'
        quality=str(100)
        filtr='PASS'
        formt='GT'
        info='NS=1'
        if out_unphased:
            hap='/'.join(hap_values)
        else:
            hap='|'.join(hap_values)

        var_out=[chrom,str(genomic_position),ID,ref_allele,alt_allele,quality,filtr,info,formt,hap]
        vcf_file.write('\t'.join(var_out)+'\n')
    vcf_file.close()
    return 1




if __name__ == "__main__":

    haplos_file_address=argv[1]  # input name 'sim_varianthaplos.txt'#
    var_list=read_haplo(haplos_file_address)

    vcf_file_address=argv[2]  # output name 'out2.vcf'#
    out_unphased=True # if you need unphased,

    result=writ_vcf(vcf_file_address,out_unphased)


    
    

