#!/usr/bin/python

import re, sys, os
from collections import defaultdict
import time
import argparse

#python s06_IR_analysis.py -tsvpath ../output/Majiq_PSI_All/tsvFiles/
#python s06_IR_analysis.py -tsvpath ../output/Majiq_PSI_All/tsvFiles/ -rcfile ../output/181206_ReadCount_Min2_formatted.csv -outdir ../output/altSS_s06_output/

def parse_readcount_file(filename):
	rcDict = defaultdict(list)
	count = 0
	with open(filename, "r") as fIN:
		for line in fIN:
			line = line.rstrip("\n")
			data = line.split("\t")
			rcDict[data[1]] = list(map(int,data[2:]))
			#print(data[2:])
	return rcDict

def get_coord(sscoord):
	ssInfo = re.search(r'(.*):(.*)-(.*)', sscoord)
	ssArr = [int(ssInfo.group(i)) for i in range(2,4) ]
	return(ssArr)

def parse_majiq_tsv(filename):
	ssDict = defaultdict(dict)
	ssType = defaultdict(list)
	ssNames = defaultdict(list)
	with open(filename, "r") as fIN:
		for line in fIN:
			if not(line.startswith("#")):
				line = line.rstrip("\n")
				data = line.split("\t")
				if (data[5].endswith("i")):
					psi = list(map(float,data[3].split(";")))
					sscoord = data[14].split(";")
					chrom = data[12]; sign = data[13]
					for i,k in enumerate(psi):
						SSt = "{}:{}".format(chrom,sscoord[i])
						k = round(float(k), 4)
						#print(data[2],SSt,psi[i])
						ssDict[data[2]][SSt] = k
						ssType[data[2]] = [ data[6],data[7],data[8] ] 
	
	return (ssDict, ssType)

def apply_filters(allLSV, rcDict):
	filteredIndex = defaultdict(dict)
	count = 0
	count1 = 0
	for lsvid in allLSV:
		flag = 0; intronFlag = 0
		ssIndex = defaultdict(list)
		for ssid in allLSV[lsvid]:
			psiVal = allLSV[lsvid][ssid]
			nk = [i for i,k in enumerate(psiVal) if (float(k) > 0.05)] #PSI cut-off applied here
			if (len(nk) > 0):
				if (ssid in rcDict):
					index = [k for k in nk if(rcDict[ssid][k] >= 2)] #Read cut-off applied here
					if (len(index) > 0):
						ssIndex[ssid] = index
						flag += 1
				else:
				#	#Since Intron Retention events are not present in rcDict file, 
				#	#I verified them based on coordinates and included them in the analysis
					ssInfo = get_coord(ssid)
					if ((ssInfo[1] - ssInfo[0]) == 1):
						ssIndex[ssid] = nk
						intronFlag += 1
						flag += 1
		if (intronFlag) & (flag >= 2):
			for ssid in ssIndex:
				filteredIndex[lsvid][ssid] = ssIndex[ssid] #index is the list of time-point index that qualify all filters: Used for binaryLSV
			count += 1
	print(count)
	#print("TP_13.5: {}".format(count1))
	return filteredIndex
	


if (__name__ == "__main__"):
	parser = argparse.ArgumentParser(prog="Get all alternative SS events from Diurnal Cycle")
	parser.add_argument('-tsvpath', dest="tsvpath", help="Enter TSV path to all tsv files")
	parser.add_argument('-rcfile', dest="rcfile", help="File with readcounts at each timepoint")
	parser.add_argument('-outdir', dest="outdir", help="File with readcounts at each timepoint")
	args = parser.parse_args()

	tsvpath = args.tsvpath
	rcfile = args.rcfile
	outdir = args.outdir
	
	date = time.strftime("%m%d")

	#Read the readcounts file and store information in dictionary
	rcDict = parse_readcount_file(rcfile)
	
	#Parse tsv files
	tsvFiles = os.listdir(tsvpath)
	tsvOrder = ["TP_1", "TP_2", "TP_3", "TP_4", "TP_5", "TP_6", "TP_7", "TP_8", "TP_9", "TP_10", "TP_11", "TP_11.5", "TP_12", "TP_12.5", "TP_13.5", "TP_14", "TP_14.5", "TP_15", "TP_16", "TP_17", "TP_18", "TP_19","TP_20", "TP_21", "TP_22", "TP_23", "TP_24"]
	allLSV = defaultdict(dict)
	lsvType = defaultdict(list)
	listOfZero = [0] * len(tsvOrder)
	for k,sampleID in enumerate(tsvOrder):
		filename = "{}.psi.tsv".format(sampleID)
		print(filename)
		filepath = "{}/{}".format(tsvpath,filename)
		#Get ssDict and ssType from parse_majiq_tsv function
		ssDict,ssType = parse_majiq_tsv(filepath)
		print(len(allLSV.keys()))	
		for lsvid in ssDict:
			if (lsvid in allLSV):
				for ssID in ssDict[lsvid]:
					if (ssID in allLSV[lsvid]):
						allLSV[lsvid][ssID].append(ssDict[lsvid][ssID])
					else:
						l = ['0.0001'] * (k+1)
						l[k] = ssDict[lsvid][ssID]
						allLSV[lsvid][ssID] = l
			else:
				allLSV[lsvid] = dict()
				lsvType[lsvid] = ssType[lsvid]
				for ssID in ssDict[lsvid]:
					l = ['0.0001'] * (k+1)
					l[k] = ssDict[lsvid][ssID]
					allLSV[lsvid][ssID] = l
		#All entries that are not filled with PSI value are filled with 'N'
		otherIDs = set(allLSV.keys()) - set(ssDict.keys())
		for oID in otherIDs:
			for ssID in allLSV[oID]:
				allLSV[oID][ssID].append('0.0001')
	count = 0
	for lsvid in allLSV:
		for ssid in allLSV[lsvid]:
			ssLen = len(allLSV[lsvid][ssid])
			print(lsvid,ssid,ssLen,allLSV[lsvid][ssid])
	filteredIndex = apply_filters(allLSV, rcDict)

	tpList = ",".join(tsvOrder)
	if not(os.path.isdir(outdir)):
		os.makedirs(outdir)
	
	fout1 = open("{}/{}_Qualified_LSV_Summary.csv".format(outdir,date), "w")
	fout1.write("LSV_ID,SS_ID,{}\n".format(tpList))
	for lsvid in allLSV:
		if(lsvid in filteredIndex):
			for ssid in allLSV[lsvid]:	
				ssArr = allLSV[lsvid][ssid]
				ssLine = ",".join(map(str,ssArr))
				#Write the PSI values for each LSV id in file_1
				fout1.write("{},{},{}\n".format(lsvid,ssid,ssLine))	
		
