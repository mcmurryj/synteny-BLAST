#!/usr/bin/python
import os
from shutil import rmtree

def make_folders(path):
    ###MAKE FOLDER FOR OUTPUT###
    if os.path.exists(path):
    	rmtree(path)
    os.makedirs(path)
