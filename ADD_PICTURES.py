'''Script to add graphics to already computed SPLORT results.'''
#########Import standard stuff##############
from os.path import abspath
from os import makedirs
from os.path import isdir
import argparse
##########Import SPLORT routines############
from parse_BLAST2 import parse
from rpsBLAST_tools import add_annot
from write_networks import write_networks
from plasmid_plotting import printClusterPics

########function for getting parameters#####
def get_parm():
    '''Get parameters for plotting figures.'''
    parser = argparse.ArgumentParser()
    parser.add_argument("-blastdir",
                        help="path to the directory BLAST2_data.")
    parser.add_argument("--output",
                        default="./",
                        help="Folder to store the pdf files.")
    args = parser.parse_args()
    return args

#Get args, make output dir if it doesn't exist already
args = get_parm()
blastoutputfile = abspath(args.blastdir + "2ndBLAST_output.XML")
rpsBLAST_output = abspath(args.blastdir + "rpsBLAST_output.XML")
output_dir = abspath(args.output)
if isdir(output_dir):
    pass
else:
    makedirs(output_dir)

#Parse the BLAST output and stor in data_box
data_box = parse(blastoutputfile)

#Grab the RPS-blast annotations
annot_dict, annot_def_dict = add_annot(data_box, rpsBLAST_output)

#Get the network results; could rewrite to just pull the network results from the gml file but this is probably fine.
cluNet = write_networks(data_box)

#put the output pdfs somewhere.
printClusterPics(cluNet, data_box, annot_def_dict, output_dir, ncutoff=8)
