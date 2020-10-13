#filter_vcf by chro
import sys

#in 1 is vcf to filter
#in 2 is locus list from script 2
vcf_in = sys.argv[1]
fh_in = open(vcf_in,'r')

vcf_out = sys.argv[1].split('.')[0] + 'filt' + '.vcf'
fh_out = open(vcf_out,'a')

loci_keep = []
with open(sys.argv[2], 'r') as filehandle:
    loci_keep = [locus.rstrip() for locus in filehandle.readlines()]

for line in fh_in:
	line = line.strip()
	if line.startswith('#'):
		fh_out.write(line+'\n')
	else:
		if line.split('\t')[0] in loci_keep:
			fh_out.write(line+'\n')