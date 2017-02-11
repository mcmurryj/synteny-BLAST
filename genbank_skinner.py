#!/usr/bin/python
###Genbank skinner
from glob import glob
from Bio import SeqIO
from sys import argv
from os.path import abspath
#import timeit

delim = "|"
dir = abspath(sys.argv[1])
gb_list = glob(dir + "*.gbff")
#gb_list =["GCF_000242715.1_ASM24271v2_genomic.gbff"]
#start_time = timeit.default_timer()
for gbfile in gb_list:
	gb = SeqIO.parse(gbfile, "genbank")					###parse needed for multirecord
	for rec in gb:
		if len(rec.seq) < 10000:						###Kill shitty contigs
			break
		for feat in rec.features:						###note: feature parsing works for refseq records from ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/
			if 'protein_id' in feat.qualifiers:			###Not funtional with EG JGI genomes, yet.
				#organism = rec.annotations['source']	#other note-"protein_id" in feat.qualifiers is faster vs. feat.qualifiers.keys() cause that makes a new list.
				print(">"+ rec.annotations['source'].replace(" ", "_")   #replacing whitespace helps to parse for hmmscan
				 + delim + rec.id + delim + feat.qualifiers['protein_id'][0]
				 + delim + feat.qualifiers['product'][0].replace(" ", "_")
				 + delim + str(feat.location.start)
				 + "-" + str(feat.location.end))#(">"+gbfile+delim+ pid+delim+locstart+"-" +locend)
				print((feat.extract(rec.seq).translate()))                        #use [:-1]to get rid of the * for stop codon, don't know if blastdb will like that

# code you want to evaluate
#elapsed = timeit.default_timer() - start_time
#print(elapsed)
