import os
import sys

#dataset number (equal to num of loci used in my naming convention, from the filtered plink)
ds = sys.argv[1]
old = str(int(ds)+1)

loci = []

map=open('HI.in.'+ds+'.map','r')
for line in map:
	line=line.strip()
	loci.append(line.split('\t')[0])
	
freq = open('freq_out_filt_'+old+'.freq','r')
fh_out= open('freq_out_'+ds+'.freq','a')
fh_out.write("Locus"+'\t'+"Allele"+'\t'+'P1'+'\t'+'P2'+'\n')

for line in freq:
	if line.split('\t')[0] in loci:
		fh_out.write(line)
	