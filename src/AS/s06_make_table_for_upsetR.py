#!/usr/bin/python

##############################################
#	Purpose: Make a 1 and 0 matrix for each gene that is 
#		readable by upsetR
#	Author: Manishi Pandey
#	Date: October 6, 2018
################################################

import re, sys, os
from os.path import basename
from collections import defaultdict
import argparse
import numpy as np
import time

def parse_tsv_file(fullpath):
	geneList = list()
	with open(fullpath, "r") as fIN:
		for line in fIN:
			if not(line.startswith("#")):
				line = line.rstrip("\n")
				data = line.split("\t")
				geneList.append(data[0])
	return geneList

if (__name__ == "__main__"):
	date = time.strftime("%m%d")
	parser = argparse.ArgumentParser(prog="Make upsetR matrix from tsvFiles of MAJIQ output")
	parser.add_argument("-tsvpath", dest="tsvpath", help="Enter path to the tsv files")
	
	args = parser.parse_args()
	tsv_path = args.tsvpath

	diffGenes = defaultdict(list)
	fileList = os.listdir(tsv_path)
	print(fileList)
	uniqGene = list()
	for filename in fileList:
		fullpath = "{}/{}".format(tsv_path,filename)
		pts = re.search(r'TP_(.*)_(TP_.*).deltapsi.tsv', filename)
		tp2 = pts.group(2)
		geneL = parse_tsv_file(fullpath)
		diffGenes[tp2] = geneL
		for gene in geneL:
			if (gene not in uniqGene):
				uniqGene.append(gene)
	tpIDs = sorted(diffGenes.keys())
	mat = np.zeros((len(uniqGene),len(tpIDs)))
	
	for i,gene in enumerate(uniqGene):
		for j,val in enumerate(tpIDs):
			if (gene in diffGenes[val]):
				mat[i][j] = 1
	
	fout = open("../output/s06_output/{}-outfile.txt".format(date), "w")
	tpStr = "\t".join(tpIDs)
	fout.write("GeneName\t{}\n".format(tpStr))
	for i,arr in enumerate(mat):
		arrStr = "\t".join(map(str,arr))
		fout.write("{}\t{}\n".format(uniqGene[i],arrStr))	
		
