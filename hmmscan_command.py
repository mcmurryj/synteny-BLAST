#!/usr/bin/python
    #problem-to-deal-with: hmmscan parses spaces in faa sequence titles.
    #will have to call iteratively on each sequence or reformat titles.
def run_hmmscan(fasta, output_dir,
                pfam_db = os.path.abspath("/media/mchanglab/shared-disk/shared-drive/Jon/pfam/Pfam-A.hmm"))
    """Run the hmmscan program on the specified fasta file, and put the results in the specified output location"""
    #Set directory names and make dirs:
    hmmscan_dir = os.path.abspath(output_dir + "/hmmscan_data")
    os.makedirs(hmmscan_dir)
    hmmscan_file = os.path.abspath(hmmscan_dir + "/hmmscan_out.tab")
    print("Made dir " + hmmscan_dir + " to hold file " + hmmscan_file)

    #define the command to be run:
    hmmscan_cmd = (
    #program name; compact output; allow very long lines of text; output tabular file
    "hmmscan --noali --notextw --tblout "
    #name of file to be outputted; pad with a space
    + hmmscan_file + " "
    #location of PFAM db
    + pfam_db
    #Location of target sequences:
    + secondary_BLAST_seq_file
    )

    print("Now running command:  " + hmmscan_cmd)
    subprocess.call(make2ndBLASTdbcmd, shell = True)
    print("Done running hmmscan!!!")
