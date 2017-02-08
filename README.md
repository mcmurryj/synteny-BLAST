��# s y n t e n y - B L A S T  Synteny blast is a tool to identify conserved genetic context(s) around a homologues of a protein of interest. It functions as follows:
1.  Print multiple fasta file from genbank files using genbank_skinner.py.  This preserves context info associated with the CDS.
2.  Make the resulting fasta into a blast db.  Not automated currently.
3.  Search query seq. against DB.  Important! Warning! If the query seq is not in the DB, results may be messed up.
4.  Parse results.  For all good hits, store proteins on either side of the hit.


 #TODO:
#Add call to hmmscan at the end of the script to enhance annotation.
#Switch making of the networks to before the printing the tables
#
#Make new DBs:
  #full bacteria DB
  #test db for running quick checks
#Usability and beautification:
  #split script into chunks
  #Add writing of intermediate outputs

#Performance enhancement:
  #DIY fasta parser?

#Aesthetics:
  #alternate strategy vs. the ungodly ID storage approach?
  #Does multigeneblast do similar? Can I use their custom DBs?
