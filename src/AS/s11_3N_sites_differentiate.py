#!/usr/bin/python

import re, sys, os
from collections import defaultdict
import argparse

def parse_coord_file(filename):
	ssDict = defaultdict(float)
	with open(filename, "r") as fIN:
		for line in fIN:
			if not(line.startswith("#")):
				line = line.rstrip("\n")
				data = line.split("\t")
				scoord = data[1].split(";")
				for i,val in enumerate(sscoord):
					ssID = "{}@{}".format(val,data[0]) #data[2] is LSV-ID
					ssDict[ssID] = val
	return ssDict

def parse_PSI_matrix(filename, lsvlist):
	psiMat = defaultdict(list)
	with open(filename, "r") as fIN:
		for line in fIN:
			if not(line.startswith("#")):
				data = line.split(",")
				if (data[0] in lsvlist):
					psiMat[data[0]] = data[1:]
	return psiMat
	
#def get_3N_sites(filename):

if (__name__ == "__main__"):
	parser = argparse.ArgumentParser(prog = "Differentiate 3N sites from Non-3N")
	parser.add_argument('-lsvfile', dest='lsvfile', help="Selected LSV IDs")
	parser.add_argument('-psifile', dest='psifile', help="PSI matrix file")
	parser.add_argument('-coordfile', dest='coordfile', help="File with LSV IDs and coordinates")
	
	args = parser.parse_args()
	
	lsvfile = args.lsvfile
	psifile = args.psifile
	coordfile = args.coordfile

	parse_coord_file(coordfile)
	parse_PSI_matrix(psifile)
	

		
