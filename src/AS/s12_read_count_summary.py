#!/usr/bin/python

#Files Needed: 1015_Majiq_All_RC, ../data/timepoints.txt
#
#python s12_read_count_summary.py -tpfile ../data/timepoints.txt -rcfile ../output/1015_Majiq_All_RC.txt

from __future__ import division
import re, sys, os
from collections import defaultdict
from collections import OrderedDict
import argparse
import time

def timepoint_mapping(filename):
	timepts = defaultdict(list)
	tpArray = list()
	with open(filename, "r") as fIN:
		for line in fIN:
			line = line.rstrip("\n")
			data = line.split("\t")
			timepts[data[0]] = "{}_1".format(data[2])
			timepts[data[1]] = "{}_2".format(data[2])
			tpArray.append(data[2])
	print(tpArray)
	return (timepts, tpArray)

def parse_reads_file(filename, timepts):
	readDict = defaultdict(dict)
	with open(filename, "r") as fIN:
		for line in fIN:
			line = line.rstrip("\n")
			data = line.split("\t")
			tpVal = timepts[data[4]]
			gen = re.search(r'(.*).v5.5', data[3])
			geneName = gen.group(1)
			chrID = "{}@{}:{}-{}".format(geneName,data[0],data[1],data[2])
			readDict[chrID][tpVal] = data[5]
	return readDict

if (__name__ == "__main__"):
	parser = argparse.ArgumentParser(prog="All New things in the world")
	parser.add_argument('-tpfile', dest='tpfile', help="Timepoints file that map timepoint to seqID")
	parser.add_argument('-rcfile', dest='rcfile', help="Read Counts file")
	
	args = parser.parse_args()

	tpfile = args.tpfile
	rcfile = args.rcfile
	date = time.strftime("%m%d")

	tps, tpArr = timepoint_mapping(tpfile)
	rcDict = parse_reads_file(rcfile, tps)
	
	print(len(rcDict.keys()))
	readsArr = defaultdict(list)
	for ssID in rcDict:
		tpSS = rcDict[ssID].keys()
		tpAvgArr = list()
		for tp in tpArr:
			tp1 = "{}_1".format(tp); tp2 = "{}_2".format(tp)
			if (tp1 in tpSS) & (tp2 in tpSS):
				tpAvg = (int(rcDict[ssID][tp1]) + int(rcDict[ssID][tp2]))/2;
			elif (tp1 in tpSS):
				tpAvg = int(rcDict[ssID][tp1])/2
			elif (tp2 in tpSS):
				tpAvg = int(rcDict[ssID][tp2])/2
			else:
				tpAvg = 0
			tpAvgArr.append(tpAvg)
		readsArr[ssID] = tpAvgArr
	print(len(readsArr.keys()))
	fout = open("../output/{}_ReadCount_In_TP_format.csv".format(date), "w")
	for ssID in readsArr:
		gene = ssID.split("@")
		rcArr = readsArr[ssID]
		#if (max(rcArr) > 1):
		rcStr = "\t".join(map(str,readsArr[ssID]))
		fout.write("{}\t{}\t{}\n".format(gene[0],gene[1],rcStr))
		
