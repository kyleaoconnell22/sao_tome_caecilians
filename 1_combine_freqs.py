#!/usr/bin/env python
import argparse
'''
calculated geno freqs using vcftools
vcftools --vcf input_u.vcf --maf 0.2 --recode --keep North.txt (or South.txt) --freqs
keep file is just list of samples to use for calcs
output is two files named for the pops

vcftools outputs should be freq_out_pop1.frq and freq_out_pop2.frq 
Names need to match the inputs given in the arguments

Also do a control replace for 0 such that it is 0.00 and for 1 such that it is 1.00
I was having issues converting between string and floats

Desired output should be as follows
Locus	Allele	P1	P2
Locus_100014	1	0.8	0.2
Locus_100014	2	0.5	0.5
Locus_100015	1	1	0
Locus_100015	2	0	1

'''

#define arguments
def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-p1", "--pop1", required=True, help="REQUIRED: the name of the pop1 from the parentals file")
    parser.add_argument("-p2", "--pop2", required=True, help="REQUIRED: the name of the pop2 from the parentals file")
    parser.add_argument("-l", "--num_loci", required=True, help="REQUIRED: number of loci in the input file")

    return parser.parse_args()
    
def calc_freqs(pop1,pop2,num_loci):
	#get file names
	file1 = 'freq_out_' + pop1 + '.frq'
	file2 = 'freq_out_' + pop2 + '.frq'
	file_out = 'freq_out_combined.freq'
	
	#open input and output files
	fh_p1 = open(file1,'r')
	fh_p2 = open(file2,'r')
	fh_out = open(file_out,'a')
	
	#write output header
	fh_out.write("Locus"+'\t'+"Allele"+'\t'+'P1'+'\t'+'P2'+'\n')
	
	#do stuff with the contents of the inputs
	#write 4 lists one for each pop and each allele and make sure the order is the same
	pop1_allele1 = []
	pop2_allele1 = []
	pop1_allele2 = []
	pop2_allele2 = []
	
	#create list of the two input files
	#iterate over the list
	#I made two loops, one for each pop input
	for line in fh_p1:
		line=line.strip()
		line = line.split('\t')
		locus = line[0]
		freq1 = line[4].split(':')[1]
		freq2 = line[-1].split(':')[1]
					
		#write to the lists

		pop1_allele1.append(locus + '\t' + freq1)
		pop1_allele2.append(locus + '\t' + freq2)
			
	for line in fh_p2:
		line=line.strip()
		line = line.split('\t')
		locus = line[0]
		freq1 = line[4].split(':')[1]
		freq2 = line[-1].split(':')[1]
			
		pop2_allele1.append(locus + '\t' + freq1)
		pop2_allele2.append(locus + '\t' + freq2)
	

	#write to the output file
	for i,j,k,m in zip(pop1_allele1,pop2_allele1,pop1_allele2,pop2_allele2):
		fh_out.write(i.split('\t')[0] + '\t' + '1' + '\t' + i.split('\t')[1] + '\t' + j.split('\t')[1] + '\n' + 
		i.split('\t')[0] + '\t' + '2' + '\t' + k.split('\t')[1] + '\t' + m.split('\t')[1] + '\n')
	
def main():
	args = get_args()
	calc_freqs(args.pop1,args.pop2,args.num_loci)
		
if __name__ == '__main__':
    main()
	
	