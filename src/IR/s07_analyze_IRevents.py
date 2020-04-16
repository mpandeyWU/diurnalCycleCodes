#!/usr/bin/python

import re, sys, os
from collections import defaultdict
import argparse

def get_coord(sscoord):
	ssInfo = re.search(r'(.*):(.*)-(.*)', sscoord)
	ssArr = [int(ssInfo.group(i)) for i in range(2,4) ]
	ssArr.append(ssInfo.group(1))
	return(ssArr)

def parse_IR_file(filename):
	fout = open("./Only_IR_hits.bed", "w")
	with open(filename, "r") as fIN:
		for line in fIN:
			if (line.startswith("Cre")):
				line = line.rstrip("\n")
				data = line.split(",")
				coord = get_coord(data[1])
				if ((coord[1] - coord[0]) == 1):
					fout.write("{}\t{}\t{}\t{}\n".format(coord[2],coord[0],coord[1],data[0]))
	fout.close()
				

if (__name__ == "__main__"):
	parser = argparse.ArgumentParser(prog="Analyze IR events")
        parser.add_argument('-infile', dest="infile", help="Enter TSV path to all tsv files")
        #parser.add_argument('-rcfile', dest="rcfile", help="File with readcounts at each timepoint")
        #parser.add_argument('-outdir', dest="outdir", help="File with readcounts at each timepoint")
	args = parser.parse_args()

	inFile = args.infile

	parse_IR_file(inFile)	
