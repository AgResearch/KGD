#input:
#           infile.vcf
#           vcf file, tab delimited
#
#output:
#           infile.vcf.ra.tab
#           RA (ReferenceAlternative) file, tab-delimited with columns, CHROM, POS, SAMPLES
#           CHROM   POS     999220   999204
#           1       415     0,0     0,0
#           1       443     1,0     9,0
#           1       448     0,0     0,0
#
#processing:
#           indels are removed and reported in infile.vcf.indel
#           SNP that are more than biallelic are removed and reported in infile.vcf.polyallele
#           all redundant sites are removed and reported in infile.vcf.posred
#           ./. is translated into 0,0

import sys, os

if len(sys.argv) != 2:
	sys.exit('Usage: python %s vcf_file\n' % sys.argv[0])

if not os.path.exists(sys.argv[1]):
	sys.exit('ERROR: vcf_file \'%s\' was not found!\n' % sys.argv[1])

infile = sys.argv[1]
infh = open(infile)

outfile = infile + '.ra.tab'
ofh = open(outfile, 'w')

indelfile = infile + '.indel'			#indels
indelfh = open(indelfile, 'w')
indels = 0					#rows with indels

polyallelefile = infile + '.polyallele'		#sites that are more than biallelic
polyallelefh = open(polyallelefile, 'w')
polyallele = 0					#polyallelic rows

posredfile = infile + '.posred'			#sites with redundant positions
posredfh = open(posredfile, 'w')
posred = 0					#rows with redundant positions

#scan input file for redundant positions
########################################

print('Scanning SNP positions')
line = infh.readline()
pos_seen = {}

while line:
	line = line.strip()
	if line.startswith('#'):
		pass
	else:
		line = line.split('\t')
		chrom = line[0]
		if not chrom in pos_seen:
			pos_seen[chrom] = {}
		pos = line[1]
		if not pos in pos_seen[chrom]:
			pos_seen[chrom][pos] = 0
		pos_seen[chrom][pos] += 1
	line = infh.readline()
infh.close()

counter = 0
for i in pos_seen:
	for j in pos_seen[i]:
		if pos_seen[i][j] > 1:
			counter += 1

print('\tFound %s redundant positions\n' %counter)

#filter input file for indels, polyallelic SNPs, redundant positions
#####################################################################

infh = open(infile)

headerlist = ['CHROM', 'POS']
line = infh.readline()

snp_counter = 0
line_counter = 0

#initialize various options for depth related fields
ad_pos = ''
ro_pos = ''
ao_pos = ''
dp4_pos = ''

while line:
	line = line.strip()
	outlist = []
	annotlist = []
	empty_genotypes = ['./.', '.,.', '.', '.|.']
	if line.startswith('##'):
		pass
	elif line.startswith('#CHROM'):
		line = line.split('\t')
		headerlist += line[9:]
		ofh.write('%s\n' %('\t'.join(headerlist)))
		print('Found %s samples' %(len(headerlist) - 2))
		for i in headerlist:
			if ' ' in i:
				print('WARN: spaces in sample names are discouraged %s' %i)
			else:
				pass
	else:
		line_counter += 1
		line = line.split('\t')
		chrom = line[0]
		outlist.append(chrom)
		annotlist.append(chrom)
		pos = line[1]
		if pos_seen[chrom][pos] > 1:				#filter out redundant positions
			posredfh.write('%s\n' %('\t'.join(line)))
			posred += 1
		else:
			outlist.append(pos)
			annotlist.append(pos)
			ref = line[3]
			alt = line[4]
			if ref == '.' or alt == '.':			#filter out indels
				indelfh.write('%s\n' %('\t'.join(line)))
				indels += 1
			elif alt.count(',') > 0:			#filter sites that are more than biallelic
				polyallelefh.write('%s\n' %('\t'.join(line)))
				polyallele += 1
			elif len(ref) > 1 or len(alt) > 1:
				indelfh.write('%s\n' %('\t'.join(line)))
				indels += 1
			else:
				if line_counter == 1:			#only look at first SNP to determine depth field
					format = line[8].split(':')
					if "AD" in format:
						ad_pos = format.index('AD')
						print('Using AD field for depth information')
					elif "RO" and "AO" in format:
						ro_pos = format.index('RO')
						ao_pos = format.index('AO')
						print('Using RO and AO fields for depth information')
					elif "DP4" in format:
						dp4_pos = format.index('DP4')
						print('Using DP4 field for depth information')
					else:
						print("\nERROR: We can't use this vcf file. AD (Allelic Depth), or RO (Reference allele observation count) and AO (Alternate allele observation count) information, or DP4 is needed.\n")
						sys.exit()

				for i in line[9:]:
					if i in empty_genotypes:	#translate uncovered to 0,0
						outlist.append('0,0')
					else:
						i = i.split(':')
						if ad_pos or ad_pos==0:
							if i[ad_pos] in empty_genotypes:
								outlist.append('0,0')
							else:
								outlist.append(i[ad_pos])
						elif (ro_pos and ao_pos) or (ro_pos == 0 and ao_pos==0):
							ad = str(i[ro_pos]) + ',' + str(i[ao_pos])
							outlist.append(ad)
						elif (dp4_pos) or dp4_pos==0:
							counts = i[dp4_pos].split(',')
							allele1 = int(counts[0]) + int(counts[1])
							allele2 = int(counts[2]) + int(counts[3])
							ad = str(allele1) + ',' + str(allele2)
							outlist.append(ad)
						else:
							##Should never really get here, but if AD, AO and RO are all null, it will break the script
							print("\nERROR: Can't find either an Allele Depth (AD) or  RO (Reference allele observation count) and AO (Alternate allele observation count) or DP4 at this position.\n")
							sys.exit()
				ofh.write('%s\n' %('\t'.join(outlist)))
				snp_counter += 1
	line = infh.readline()

infh.close()
ofh.close()

print('\n')
print('%12s SNPs written to %s' %(snp_counter, outfile))
print('%12s SNPs (rows) ignored because of indels and written to %s' %(indels, indelfile))
print('%12s SNPs (rows) ignored because of more than two alleles and written to %s' %(polyallele, polyallelefile))
print('%12s SNPs (rows) ignored because of redundant genomic positions and written to %s' %(posred, posredfile))

