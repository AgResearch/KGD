#input: vcf file, tab delimited
#output: RA (ReferenceAlternative) file, tab-delimited with columns, CHROM, POS, SAMPLES
# where SAMPLES consists of colon-delimited sampleID, flowcellID, lane, seqlibID
# CHROM   POS     999220:C4TWKACXX:7:56   999204:C4TWKACXX:7:56
# 1       415     0,0     0,0
# 1       443     1,0     9,0
# 1       448     0,0     0,0
#
#filters: indels are removed, multiple alternative alleles are removed
#./. is translated into 0,0

#NB: Ammendments made by Rachael Ashby to vcf2ra.py script to allow for 'RO' and 'AO' values in a 
#VCF file to be used in place of 'AD'
#Date: 15/03/2016

import sys, os

#check python version
PY2 = sys.version_info[0] == 2
PY3 = sys.version_info[0] == 3

##Should this be greater than?
if len(sys.argv) < 2:
	sys.exit('Usage: python %s vcf_file\n' % sys.argv[0])

if not os.path.exists(sys.argv[1]):
	sys.exit('ERROR: vcf_file \'%s\' was not found!\n' % sys.argv[1])

infile = sys.argv[1]
infh = open(infile)

outfile = infile + '.ra.tab'
ofh = open(outfile, "w")

headerlist = ['CHROM', 'POS']
line = infh.readline()
snp_counter = 0

while line:
	line = line.strip()
	outlist = []
	annotlist = []
	ad_pos = "" ##Needs to be defined so can check if the string is NULL
	ro_pos = "" ##Needs to be defined so can check if the string is NULL
	ao_pos = "" ##Needs to be defined so can check if the string is NULL
	dp4_pos = "" ##Needs to be defined so can check if the string is NULL
	empty_genotypes = ['./.', '.,.', '.', '.|.']
	if line.startswith('##'):
		pass
	elif line.startswith("#CHROM"):
		line = line.split('\t')
		headerlist += line[9:]
		if PY2:
			print >> ofh, '\t'.join(headerlist)
			print "found", len(headerlist), "samples"
		if PY3:
			print('\t'.join(headerlist), file=ofh)
			print("found", len(headerlist), "samples")
		for i in headerlist:
			if " " in i:
				if PY2:
					print "WARN: spaces in sample names are discouraged", i
				if PY3:
					print("WARN: spaces in sample names are discouraged", i)
			else:
				pass
	else:
		line = line.split('\t')
		chrom = line[0]
		outlist.append(chrom)
		annotlist.append(chrom)
		pos = line[1]
		outlist.append(pos)
		annotlist.append(pos)
		ref = line[3]
		alt = line[4]
		if ref == '.' or alt == '.' or len(alt) > 1:	#filter out indels and multiple alternative alleles
			pass
		else:
			format = line[8].split(':')
			if "AD" in format:
				ad_pos = format.index('AD')     #where is AD?
			elif "RO" and "AO" in format:		#Added line here to look for RO and AO as well as AD. Could possibly combine with first if?
				ro_pos = format.index('RO')	#Where is RO?
				ao_pos = format.index('AO')	#Where is AO?
			elif "DP4" in format:
				dp4_pos = format.index('DP4')
			else:
				if PY2:
					print "\nERROR: We can't use this vcf file. AD (Allelic Depth) or RO (Reference allele observation count) and AO (Alternate allele observation count) information is needed.\n"
				if PY3:
					print("\nERROR: We can't use this vcf file. AD (Allelic Depth) or RO (Reference allele observation count) and AO (Alternate allele observation count) information is needed.\n")
				sys.exit()
			for i in line[9:]:
				if i in empty_genotypes:	#translate uncovered to 0,0. Added translations for '.,.' and '.' to allow for
                                                                #VCF files with other formats
					outlist.append('0,0')
				else:
					i = i.split(':')
					if ad_pos or ad_pos==0:				#IF ad_pos is not null, i.e., there is a value for it, append it to outlist
						if i[ad_pos] in empty_genotypes:        #gatk vcf4.2 will fill out genotype fields, even for uncovered data
							outlist.append('0,0')
						else:
							outlist.append(i[ad_pos])
					elif (ro_pos and ao_pos) or (ro_pos == 0 and ao_pos==0):	#ELSE IF ro_pos and ao_pos are not null or are equal to 0
						ad = str(i[ro_pos]) + "," + str(i[ao_pos])	#Create a string of 'RO,AO', in the same format as AD.	
						outlist.append(ad)
					elif (dp4_pos):
						counts = i[dp4_pos].split(",")
						allele1 = int(counts[0]) + int(counts[1])
						allele2 = int(counts[2]) + int(counts[3])
						ad = str(allele1) + "," + str(allele2)
						outlist.append(ad)
					else:
						##Should never really get here, but if AD, AO and RO are all null, it will break the script
						if PY2:
							print "\nERROR: Can't find either an Allele Depth (AD) or  RO (Reference allele observation count) and AO (Alternate allele observation count) at this position.\n"
						if PY3:
							print("\nERROR: Can't find either an Allele Depth (AD) or  RO (Reference allele observation count) and AO (Alternate allele observation count) at this position.\n")
						sys.exit()
			if PY2:

				print >> ofh, '\t'.join(outlist)
			if PY3:
				print('\t'.join(outlist), file=ofh)
			snp_counter += 1
	line = infh.readline()

infh.close()
ofh.close()

if PY2:
	print snp_counter, "SNPs written."
if PY3:
	print(snp_counter, "SNPs written.")
