fh = open('freq_out_combined.freq','r')
fh_out_freq = open('freq_out_filt.freq','a')
fh_out_freq.write('Locus'+'\t'+'Allele'+'\t'+'P1'+'\t'+'P2'+'\n')
fh_out=open('loci_to_keep.txt','a')

low = 0.06
high = 0.94

loci = []

for line in fh:
	line = line.strip()
	if line.startswith('Locus'):
		pass
	else:
		x=line.split('\t')
		if float(x[2]) > high and float(x[3]) < low or float(x[2]) < low and float(x[3]) > high:
			loci.append(x[0])
			fh_out_freq.write(line+'\n')
			
loci_set = set(loci)
for locus in loci_set:
	fh_out.write(locus+'\n')
	