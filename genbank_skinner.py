#!/usr/bin/python
###Genbank skinner
from glob import glob
from Bio import SeqIO
import re
from sys import argv
#import timeit

delim = "|"
# gb_list = glob("*.gbff")
gb_list = argv[1:]
#start_time = timeit.default_timer()
for gbfile in gb_list:
	gb = SeqIO.parse(gbfile, "genbank")					###parse needed for multirecord
	for rec in gb:
		if len(rec.seq) < 10000:						###Kill shitty contigs
			break
		for feat in rec.features:						###note: feature parsing works for refseq records from ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/
			if 'protein_id' in feat.qualifiers:			###Not funtional with EG JGI genomes, yet.
				organism   = re.sub('[$%^&*#|]', "", rec.annotations['source'].replace(" ", "_"))
				accession  = re.sub('[$%^&*#|]', "", rec.id)
				wp_id      = re.sub('[$%^&*#|]', "", feat.qualifiers['protein_id'][0])
				annotation = re.sub('[$%^&*#|]', "", feat.qualifiers['product'][0].replace(" ", "_"))
				ss         = str(feat.location.start) + "-" + str(feat.location.end)
				print(">" + delim.join([organism, accession, wp_id, annotation, ss]))
				print((feat.extract(rec.seq).translate()))                        #use [:-1]to get rid of the * for stop codon, don't know if blastdb will like that

# code you want to evaluate
#elapsed = timeit.default_timer() - start_time
#print(elapsed)

#: re.sub('[$%^&*#|]', "", astring)

#Mooched from multigeneblast, should include similar to sanitize the input so as to not get aberant parsing:
#Remove illegal chars
#    illegal_chars  = '''!"#$%&()*+,:;=>?@[]^`'{|} '''
#    genename = "".join([char for char in genename if char not in illegal_chars])
#    if len(genename) < 2:
#      genename = "orf" + str(nr) + "_" + str(genenr)
