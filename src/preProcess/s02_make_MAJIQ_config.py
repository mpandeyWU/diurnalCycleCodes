#!/usr/bin/python

#################################################
#	Program: To make a majiq config file to run the build command
#	Author: Manishi Pandey
#	Date: October 1, 2018
#################################################


import re, sys, os
from collections import defaultdict
import argparse
import time

def read_timepoints_file(filename):
	srrids = defaultdict(list)
	with open(filename, "r") as fIN:
		for line in fIN:
			line = line.rstrip("\n")
			data = line.split("\t")
			tp = data[1].split("_")
			if (tp[0] in srrids):
				srrids[tp[0]].append(data[0])
			else:
				srrids[tp[0]] = [ data[0] ]
	return srrids

if (__name__ == "__main__"):
	date = time.strftime("%m%d")
	parser = argparse.ArgumentParser(prog="Make MAJIQ config file")
	parser.add_argument('-file', dest="inputFile", help="Enter IDs and timepoint filename")
	parser.add_argument('-samdir', dest="samdir", help="path to sorted bam files in one directory")
	
	args = parser.parse_args()

	inputFile = args.inputFile
	samdir = args.samdir
	srr = read_timepoints_file(inputFile)
	fout = open("./{}-majiq_config.txt".format(date), "w")
	fout1 = open("./r07")
	fout.write("[info]\nreadlen=101\n")
	fout.write("samdir={}\n".format(samdir))
	fout.write("genome=v5\ngenome_path=/scratch/gslab/mpandey/lib/PhytozomeV12/Creinhardtii/assembly/\n\n")
	fout.write("[experiments]\n")
	for tp in sorted(srr):
		fout.write("Timepoint_{}={}.sorted,{}.sorted\n".format(tp,srr[tp][0],srr[tp][1]))
		
		
		
