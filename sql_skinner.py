#see: https://docs.python.org/2/library/sqlite3.html
import sqlite3
from os.path import abspath

#Get the dir with the gb files and the place to put the output.
parser = argparse.ArgumentParser()
parser.add_argument("-gb_dir",
                    help = "Directory with genbenk files, .gb or .gbff suffix."
                    "-output",
                    help = "Where to store the output sql database.  Should end in ".db".")
args = parser.parse_args()

#sqlite voodoo
cursor = initializeDb(abspath(args.output))

#get genbank files
gb_list = glob(args.gb_dir + "*.gb")
gb_list += glob(args.gb_dir + "*.gbff")

for gbfile in gb_list:
    gb    = SeqIO.parse(gbfile, "genbank")					###parse needed for multirecord
    for rec in gb:
        if len(rec.seq) < 10000: #If it sucks, skip the bottom stuff and head to next record
            continue
        else: #otherwise get data on the contig
            recData = getRecData(gb) #tuple of data, same order as for input into sql db
            cursor.execute('''INSERT INTO
                           nucleotidedb,
                           VALUES (?)''', gbData)
        ###note: feature parsing works for refseq records from ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/
        for feat in rec.features:
            if 'protein_id' in feat.qualifiers:			###Not funtional with EG JGI genomes, yet.
                featureData = getFeatureData(feat, recData[0]) #tuple
                cursor.execute('''INSERT INTO
                               proteindb,
                               VALUES (?)''', featureData)

def initializeDb(outputfile):
    conn = sqlite3.connect(outputfile)
    cursor = conn.cursor()

    #Make table for nucleotide and protein data.
    cursor.execute('''CREATE TABLE
                nucleotidedb(
                contig_id text PRIMARY KEY,
                sequence text,
                species text
                ''');

    #bear in mind that protein ID alone is not a unique identifier.
    #In fact, no guarantee that protein_id and contig_id together are unique.
    #That would be freaky, though.
    cursor.execute('''CREATE TABLE
                   proteindb,
                   protein_id text,
                   start integer,
                   stop integer,
                   contig_id text,
                   sequence text,
                   annotation text,
                   PRIMARY KEY (contig_id, protein_id)
                   FOREIGN KEY (contig_id) REFERENCES nucleotidedb(contig_id)
                   ''');
    return cursor

def getRecData(rec):
    organism   = re.sub('[\/$%^&*#|\[\]]', "", rec.annotations['organism'].replace(" ", "_"))
    #Species is first 2x chunks of the organism name
    species = "_".join(organism.split("_")[:2])
    #Pull the contig name.  Need to do for every record.
    contig_id  = re.sub('[\/$%^&*#|]', "", rec.id)
    #pull the sequence TODO
    sequence = ""
    return (contig_id, sequence, species)

def getFeatureData(feat, contig):
    protein_id = re.sub('[$%^&*#|]', "", feat.qualifiers['protein_id'][0])
    #Is feature on sense or antisense strand?
    strand = feat.strand
    if strand == 1:   #If it is on sense, record as start-end
        ss = (int(feat.location.start), int(feat.location.end))
    elif strand == -1: #If antisense, record as end-start
        ss = (int(feat.location.end), int(feat.location.start))
    sequence = "" #TODO
    if 'product' in feat.qualifiers :
        annotation = re.sub('[\/$%^&*#|]', "", feat.qualifiers['product'][0].replace(" ", "_"))
    else:
        annotation = "No annotation."
    return(protein_id, ss[0], ss[1], contig_id, sequence, annotation, species)
