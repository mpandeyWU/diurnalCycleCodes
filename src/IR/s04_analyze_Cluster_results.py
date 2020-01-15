#!/usr/bin/python

import re, sys, os
from collections import defaultdict
from collections import OrderedDict
import argparse

def parse_cluster_file(filename):
	clustInfo = defaultdict(dict)
	tsvOrder = ["TP_1", "TP_2", "TP_3", "TP_4", "TP_5", "TP_6", "TP_7", "TP_8", "TP_9", "TP_10", "TP_11", "TP_11.5", "TP_12", "TP_12.5", "TP_13.5", "TP_14", "TP_14.5", "TP_15", "TP_16", "TP_17", "TP_18", "TP_19","TP_20", "TP_21", "TP_22", "TP_23", "TP_24"]
	with open(filename, "r") as fIN:
		for line in fIN:
			if not(line.startswith("intronID")):
				line = line.rstrip("\n")
				data = line.split(",")
				pir = map(float,data[2:])
				tpMax = pir.index(max(pir))
				print(tsvOrder[tpMax])
				clustInfo[data[1]][data[0]] = data[2:]
			
	return clustInfo

if (__name__ == "__main__"):
	parser = argparse.ArgumentParser(prog="Get IR events for all important hits")
        parser.add_argument('-clustfile', dest="clustfile", help="Enter path to bed files that are mapped within a given intron")
	#parser.add_argument('-covpath', dest="covpath", help="Enter the IR coverage path")
        #parser.add_argument('-timepoints', dest="timepoints", help="Timepoints list file")
        #parser.add_argument('-outdir', dest="outdir", help="File with readcounts at each timepoint")
        args = parser.parse_args()
	
	clustfile = args.clustfile
	
	clustInfo = parse_cluster_file(clustfile)
	
	
