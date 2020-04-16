#!/usr/bin/python

import re, sys, os
from collections import defaultdict
import argparse

def parse_readcount_file(filename):
	rcDict = defaultdict(list)
	with open(filename, "r") as fIN:
		for line in fIN:
			line = line.rstrip("\n")
			data = line.split("\t")
			rcDict[data[1]] = map(int,data[2:])
	return rcDict

def parse_majiq_tsv(filename):
	ssDict = defaultdict(dict)
	ssType = defaultdict(dict)
	ssNames = defaultdict(list)
	with open(filename, "r") as fIN:
		for line in fIN:
			if not(line.startswith("#")):
				line = line.rstrip("\n")
				data = line.split("\t")
				psi = list(map(float,data[3].split(";")))
				sscoord = data[14].split(";")
				chrom = data[12]; sign = data[13]
				#Index of entries that have PSI greater than 0.05
				index = [ i for i,k in enumerate(psi) if (k > 0.05)]
				if (len(index) > 1):
					for i in index:
						SSt = "{}:{}".format(chrom,sscoord[i])
						psi[i] = round(float(psi[i]), 4)
						#print(data[2],SSt,psi[i])
						ssDict[data[2]][SSt] = psi[i]
					ssType[data[2]]["Alt5"] = data[6] 
					ssType[data[2]]["Alt3"] = data[7] 
					ssType[data[2]]["ES"] = data[8] 
	
	return (ssDict, ssType)



if (__name__ == "__main__"):
	parser = argparse.ArgumentParser(prog="Get all alternative SS events from Diurnal Cycle")
	parser.add_argument('-tsvpath', dest="tsvpath", help="Enter TSV path to all tsv files")
	parser.add_argument('-rcfile', dest="rcfile", help="File with readcounts at each timepoint")
	args = parser.parse_args()

	tsvpath = args.tsvpath
	rcfile = args.rcfile

	#Read the readcounts file and store information in dictionary
	rcDict = parse_readcount_file(rcfile)
	
	#Parse tsv files
	tsvFiles = os.listdir(tsvpath)
	tsvOrder = ["TP_1", "TP_2", "TP_3", "TP_4", "TP_5", "TP_6", "TP_7", "TP_8", "TP_9", "TP_10", "TP_11", "TP_11.5", "TP_12", "TP_12.5", "TP_13.5", "TP_14", "TP_14.5", "TP_15", "TP_16", "TP_17", "TP_18", "TP_19","TP_20", "TP_21", "TP_22", "TP_23", "TP_24"]
	allLSV = defaultdict(dict)
	summary = defaultdict(list)
	listOfZero = [0] * len(tsvOrder)
	for k,sampleID in enumerate(tsvOrder):
		filename = "{}.psi.tsv".format(sampleID)
		print(filename)
		#fname = re.search(r'(.*)psi.tsv', filename)
		#sampleID = fname.group(1)
		#print(fname.group(1),filename)
		filepath = "{}/{}".format(tsvpath,filename)
		ssDict,ssType = parse_majiq_tsv(filepath)
		
		qualLSV = list()
		for lsvid in ssDict:
			flag = 0
			ssCount = [ssID for ssID in ssDict[lsvid] if (ssID in rcDict)]
			if (len(ssCount) == len(ssDict[lsvid].keys())):
				if (lsvid not in allLSV):
					allLSV[lsvid] = dict()
					summary[lsvid] = [0] * k
				for ssID in ssDict[lsvid]:
					if (ssID not in allLSV[lsvid]):
						allLSV[lsvid][ssID] = ['NaN']*(k)
					if (rcDict[ssID][k] >= 2):
						allLSV[lsvid][ssID].append(ssDict[lsvid][ssID])
						flag += 1		
				if (flag > 1):
					summary[lsvid].append(1)
					qualLSV.append(lsvid)
				else:
					summary[lsvid].append(0)
		otherIDs = set(summary.keys()) - set(ssDict.keys())
		print(len(otherIDs))
		for oID in otherIDs:
			summary[oID].append(0)

	for lsvid in summary:
		if (len(summary[lsvid]) < 27):
			print(lsvid, summary[lsvid], len(summary[lsvid])) 
	print(len(summary.keys()))
	#for lsvid in allLSV:
	#	for ssid in allLSV[lsvid]:
	#		print(lsvid,ssid,allLSV[lsvid][ssid])
