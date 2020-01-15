#!/usr/bin/python

from __future__ import division
import re, sys, os
from collections import defaultdict
from collections import OrderedDict
import argparse
import time

def parse_AltSS_file(filename):
	ssDict = defaultdict(dict)
	with open(filename, "r") as fIN:
		for line in fIN:
			if (line.startswith("Cre")):
				line = line.rstrip("\n")
				data = line.split(",")
				lsv = re.search(r'(Cre.*).v5.5:(.*)', data[0])
				geneid = lsv.group(1)
				lsvType = data.pop()
				ssDict[geneid][data[0]] = data[2:]
	return ssDict

def parse_gff_file(gff_file):
	mRNAdict = defaultdict(dict)
	with open(gff_file, "r") as fIN:
		for line in fIN:
			if not(line.startswith("#")):
				line = line.rstrip("\n")
				data = line.split("\t")
				if (data[2] == "mRNA"):
					gene = re.search(r'ID=(Cre.*).v5.5;Name=(.*).t(.*);pacid=(.*)', data[8])
					#gene = re.search(r'ID=((.*).t(.*).v5.5.exon(.*));Parent=(.*)',data[8])
					geneName = gene.group(2)
					mRNA = gene.group(1)
					if (geneName not in mRNAdict):
						mRNAdict[geneName] = [ mRNA ]
					else:
						mRNAdict[geneName].append(mRNA)
	count = 0
	for geneid in mRNAdict:
		if (len(mRNAdict[geneid]) > 1):
			count += 1
	print(count)
	return mRNAdict

def parse_summary_file(filename):
	summary = defaultdict(list)
	with open(filename, "r") as fIN:
		for line in fIN:
			if not(line.startswith("#")):
				line = line.rstrip("\n")
				data = line.split(",")
				summary[data[0]]  = data[1:]
	return summary

def parse_IR_summary_file(filename):
	intronSum = defaultdict(dict)
	with open(filename, "r") as fIN:
		for line in fIN:
			if not(line.startswith("#")):
				line = line.rstrip("\n")
				data = line.split(",")
				gene = re.search(r'(.*).v5.5',data[1])
				intronSum[gene.group(1)][data[0]] = 1
	return intronSum
	
if (__name__ == "__main__"):
	parser = argparse.ArgumentParser(prog = "Compare AltSS genes with Phytozome transcripts")
	parser.add_argument('-ssfile', dest="ssfile", help="File with filtered SS LSVIDs and PSI values")
	parser.add_argument('-IRfile', dest="IRfile", help="File with filtered SS LSVIDs and PSI values")
	parser.add_argument('-infile', dest="infile", help="Summary file obtained from s02_ script file")
	#parser.add_argument('-gff', dest="gff", help="Chlamydomonas GFF file")
	args = parser.parse_args()

	date = time.strftime("%m%d")
	ssFile = args.ssfile
	infile = args.infile
	#gffFile = args.gff
	intronFile = "../output/IR_s03_output/Jan6/190106_PIR_filtered.csv"
	gffFile = "/scratch/gslab/mpandey/lib/PhytozomeV12/Creinhardtii/annotation/Creinhardtii_281_v5.5.gene_exons.gff3" 
	ssDict = parse_AltSS_file(ssFile)
	rnaDict = parse_gff_file(gffFile)
	summary = parse_summary_file(infile)
	intronDict = parse_IR_summary_file(intronFile)
	
	allSS = list(ssDict.keys())
	#allSS = list(ssDict.keys()) + list(set(intronDict.keys()) - set(ssDict.keys()))	
	print(len(allSS))
	multiHits = [s for s in rnaDict if len(rnaDict[s]) > 1 ]	
	#common = list(set(ssDict.keys()).intersection(multiHits))
	common = list(set(allSS).intersection(multiHits))
	uniqPhyto = list(set(multiHits) - set(common))
	uniqSS = list(set(allSS) - set(common))
	
	count = 0
	for geneID in allSS:
		if (geneID in intronDict):
			c = len(intronDict[geneID].keys())
			count = count + c
		if (geneID in ssDict):
			c = len(ssDict[geneID].keys())
			count = count + c
	print("Total Events:", count)
	print(len(common),len(uniqPhyto),len(uniqSS))
	
	#fout1 = open("../output/altSS_s05_output/{}_LSVID_Summary.csv".format(date), "w")
	#fout2 = open("../output/altSS_s05_output/{}_Gene_Summary.csv".format(date), "w")
	#for geneName in ssDict:
	#	inPhyto = None
	#	if (geneName in common):
	#		inPhyto = "Yes"
	#	else:
	#		inPhyto = "No"
	#	for ssid in ssDict[geneName]:
	#		summaryInfo = summary[ssid]
	#		sumStr = ",".join(summaryInfo)
	#		psiVals = ssDict[geneName][ssid]
	#		repPSI = [s for s in psiVals if (float(s) > 0.0001)]
	#		PSIper = len(repPSI)/27
	#		#print(summaryInfo[3])
	#		if (PSIper < 0.25 ) & (float(summaryInfo[3]) > 0.25) & (inPhyto == "No"):
	#			print("{},{},{},{}".format(ssid,sumStr,PSIper,inPhyto))
	#		fout1.write("{},{},{},{}\n".format(ssid,sumStr,PSIper,inPhyto))
	#		
		

