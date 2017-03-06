#!/usr/bin/python
#python genbank_skinner.py -gb_dir /path/to/gb/file/directory/ > fastadump.faa
#Processes  Genbank files into multiple fasta protein sequences that contain rudimentary information about genomic context.
#I believe a similar approach is taken by multigeneblast; I've taken inspiration from their routine for purging special characters.
from glob import glob
from Bio import SeqIO
import re
import argparse

#Get the dir with the gb files
parser = argparse.ArgumentParser()
parser.add_argument("-gb_dir",
                    help = "Directory with genbenk files, .gb or .gbff suffix.")
args = parser.parse_args()

#Initialize misc variables and read the files
delim = "|"
spec_counts = {}
gb_list = glob(args.gb_dir + "*.gb*")

#Loops loops loops
for gbfile in gb_list:
    gb    = SeqIO.parse(gbfile, "genbank")					###parse needed for multirecord
    #Assign variable to tell script if it is the first run thru a genbank file.
    first = True
    for rec in gb:
        #stuff you only have to do for first record, stays the same afterwards.
        if first == True:
            first = False
            #Retrieve organism name; sanitize illegal_chars; fill spaces with _
            organism   = re.sub('[$%^&*#|]', "", rec.annotations['organism'].replace(" ", "_"))
            #Species is first 2x chunks of the organism name
            spec = "_".join(organism.split("_")[:2])
            #If we've never seen it before throw it in a dict.
            if spec not in spec_counts :
                spec_counts[spec] = 0
            #If we've seen it fewer than n times and species is not "sp.",
            if spec_counts[spec] > 0 and spec.split("_")[1] != "sp.":
                #bust out of the "record" loop and go to next gb
                break
        #Check that the contig is of a reasonable length.  Need to do for every record.
        if len(rec.seq) < 10000:
            #If it sucks, skip the bottom stuff and head to next record
            continue
        #Pull the contig name.  Need to do for every record.
        accession  = re.sub('[$%^&*#|]', "", rec.id)
        for feat in rec.features:						###note: feature parsing works for refseq records from ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/
            if 'protein_id' in feat.qualifiers:			###Not funtional with EG JGI genomes, yet.
                wp_id          = re.sub('[$%^&*#|]', "", feat.qualifiers['protein_id'][0])
                if 'product' in feat.qualifiers :
                    annotation = re.sub('[$%^&*#|]', "", feat.qualifiers['product'][0].replace(" ", "_"))
                else:
                    annotation = "No annotation."
                ss             = str(feat.location.start) + "-" + str(feat.location.end)
                print(">" + delim.join([organism, accession, wp_id, annotation, ss]))
                print((feat.extract(rec.seq).translate()))                        #use [:-1]to get rid of the * for stop codon, don't know if blastdb will like that
    #Increment the species counter.
    spec_counts[spec]    += 1
    #print(spec)
    #print(spec_counts[spec])

#print(spec_counts)
