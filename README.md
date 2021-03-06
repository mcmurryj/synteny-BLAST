## synteny-BLAST 

Synteny blast is a tool to identify conserved genetic context(s) around a homologues of a protein of interest.



### Usage:



*genbank_skinner.py /a/big/list/of/genbank/files/*\*.*gbff > MY_DATABASE.faa*



*makeblastdb -in  MY_DATABASE.faa -input_type fasta -dbtype prot*



*python .\FIND_ME_NEW_PROTEINS.py -h*

*python .\FIND_ME_NEW_PROTEINS.py -input /path/to/sequence.faa -database /path/to/DB.faa -output /path/to/folder  --rps_db /path/to/rpsBLSTdb*



### Notes:

It takes ~5-10 minutes to run a normal query with all streptomyces as the DB.  Much of that time is in running rpsBLAST so if you don't specify that DB it will run faster.



For the curious or frustrated, it functions as follows:

1.  Print multiple fasta file from genbank files using genbank_skinner.py.  This preserves context info associated with the CDS.



2.  Make the resulting fasta into a blast db.  Not automated currently.



3.  Search query seq. against DB.  Important! Warning! If the query seq is not in the DB, results may be messed up.



4.  Parse results.  For all good hits, store proteins on either side of the hit.  Write these to a new fasta file, where the names contain info on the protein itself as well as the progenitor



5.  Make a secondary blast db and run all v. all blast



6.  Optionally, run rpsBLAST to tag proteins with CDD, pfam, or other domain annotations.



7.  Parse those results, keeping the context info associated with each sequence.



8.  Print tabular output.



9.  Make cluster-level and protein-level networks; write gml files.




### Terminology:


*seed sequence*:

        The original protein used as the query for the first blast search.  Generally, a protein you are interested in.  Potentially YFP.



*cluster_ID*:

            ID of protein homologous to the seed sequence.



*cluster_member_ID*:

     ID of protein located within n coding sequences of cluster_ID protein.



*homologue_ID*:

          ID of protein with homology to cluster_member_ID.



### Todo: functionality:

- Output:  print contingency table describing frequency of PFAM domain co-occurrence.



- Output:  Add function to generate graphical output.  Ideally, can generate aligned graphical output.  One possibility is to run multigeneblast, using conserved proteins as queries.  Another options is to plot the top 10 or so hits for each cluster, and color CDS by PFAM domain. Plotting with Bio.Graphics seems OK.



- Clustering:  Switch to tabular output or XGMML output to allow more flexible node/edge attributes.



- Clustering:  Some sort of automated/unsupervised clustering would be save some time in Cytoscape.  Kmodes clustering based on PFAM domains, or alternately based on the protein level graph output.  Color the cluster network.



- Cytoscape:  Hyperlink from node to tabular output or graphics?

- Idea:  If I can find a mapping from WP_ ids to pfam domains, write a csv file with all results while doing parse_BLAST2.  Use pandas to summarize from that master file.

### Todo: DBs:



- full bacteria DB

- Diverse test db for running quick tests



### TODO: performance enhancement:

- Parsing the 2ndary BLAST results and especially writing the network files may also be a performance limit in cases with a large number of very similar clusters.  Hard solution: rewrite the write-networks routine, posibly by bypassing networkx alltogether and printing tabular output that is readable by Cytoscape. Hack job options:  just have a max output <200 set. Dereplication after first blast could be good, maybe using CD-HIT.



- Switch from BLAST to diamond.  Preliminary testing indicates that DIAMOND gets bunged up by the custom protein DBs used here and runs even slower than BLASTp.  It works fine on equivalent databases with normal protein sequence names, however, which is confusing.



- Alter structure of loop-doer.  If the list of BLAST hits can be arranged in the same order as the fasta sequences in the database file, this can eliminate the use of a nested for loop.   Alternately, using cython or storing the db on an SSD may help as well.



### Aesthetics:



- Make installer?  At least a file with a list of dependencies, so that they can be installed by pip.



- Alternate strategy vs. the ungodly ID storage approach?  Essentially I have a database that has been coerced into fasta format.  



- Does multigeneblast do similar? Can I use their custom DBs?




### CAVEATS/BUGFIX:

- Certain features of genbank_skinner only work in python 3.  Modify print statements and make them write statements.



- Alter genbank_skinner to use write instead of print.
