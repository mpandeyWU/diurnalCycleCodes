#!/usr/bin/python

#########################################################
#	Purpose: Parse tsv output files and filter the relevant LSVs
#	Author: Manishi Pandey
#	Date: October 16, 2018
#	Usage: python s05_parse_PSI_vals.py -majiqpath ../output/Majiq_PSI_All/tsvFiles/ -outfile ../output/s05_output/all_LSV_list.txt
###########################################################

import re, sys, os
import argparse
from collections import defaultdict
import numpy as np
import scipy.stats as sp

def parse_majiq_tsv(filename):
	ssDict = defaultdict(float)
	with open(filename, "r") as fIN:
		for line in fIN:
			if not(line.startswith("#")):
				line = line.rstrip("\n")
				data = line.split("\t")
				psi = data[3].split(";")
				sscoord = data[14].split(";")
				chrom = data[12]; sign = data[13]
				for i,val in enumerate(sscoord):
					ssID = "{}:{}:{}@{}".format(chrom,val,sign,data[2]) #data[2] is LSV-ID
					psi[i] = round(float(psi[i]), 4)
					ssDict[ssID] = psi[i]
	return ssDict

if (__name__ == "__main__"):
	parser = argparse.ArgumentParser(prog="")
	parser.add_argument('-majiqpath', dest="majiqpath", help="Path to the output files of Majiq")
	parser.add_argument('-outfile', dest="outcsv", help="Path to the output files of Majiq")
	
	args = parser.parse_args()
	tsv_path = args.majiqpath
	csv_out = args.outcsv
	
	allSS = defaultdict(list)
	tsvFiles = os.listdir(tsv_path)
	tsvOrder = ["TP_1", "TP_2", "TP_3", "TP_4", "TP_5", "TP_6", "TP_7", "TP_8", "TP_9", "TP_10", "TP_11", "TP_11.5", "TP_12", "TP_12.5", "TP_13.5", "TP_14", "TP_14.5", "TP_15", "TP_16", "TP_17", "TP_18", "TP_19","TP_20", "TP_21", "TP_22", "TP_23", "TP_24"]
	print(sorted(tsvFiles))
	for k,sampleID in enumerate(tsvOrder):
		filename = "{}.psi.tsv".format(sampleID)
		print(filename)
		#fname = re.search(r'(.*)psi.tsv', filename)
		#sampleID = fname.group(1)
		#print(fname.group(1),filename)
		filepath = "{}/{}".format(tsv_path,filename)
		ssDict = parse_majiq_tsv(filepath)
		for ssid in ssDict:
			if (ssid in allSS):
				allSS[ssid].append(ssDict[ssid]) #allSS dictionary has summary of all ssID with lsvid with a list of psi-values for each time-point.
			else:
				l = [0.001]*(k+1)
				l[k] = ssDict[ssid]
				allSS[ssid] = l
		otherIDs = set(allSS.keys()) - set(ssDict.keys())
		for oID in otherIDs:
			allSS[oID].append(0.001)
	
	varArr = list()
	cvArr = list()
	filteredSS = defaultdict(dict)
	for sampleID in allSS:
		arrPSI = allSS[sampleID]
		#print(sampleID,arrPSI)
		if (all(x < 0.05 for x in arrPSI)) | (all(x > 0.95 for x in arrPSI)):
			next;
		else:
			varPSI = np.var(arrPSI)
			cvPSI = sp.variation(arrPSI)
			varArr.append(varPSI)
			cvArr.append(cvPSI)
			#print(arrPSI,varPSI)	
			if (varPSI >= 0.2):	#Need to choose a proper cut-off
				ids = sampleID.split("@")
				filteredSS[ids[1]][ids[0]] = arrPSI
	
	newFilter = defaultdict(dict)
	#np.histogram(varArr)
	#med = np.median(varArr)
	#mean = np.mean(varArr)
	#print(med, mean, min(varArr), max(varArr))
	#med1 = np.median(cvArr)
	#mean1 = np.mean(cvArr)
	#print(med1, mean1, min(cvArr), max(cvArr))
	allSSList = allSS.keys()
	for lsvid in filteredSS:
		#otherSS = [ x for x in allSSList if (re.search(r'(.*)@' + re.escape(lsvid), x))]
		#print(otherSS)
		#for oID in otherSS:
		#	sID = oID.split("@")
		#	newFilter[lsvid][sID[0]] = allSS[oID]
		for ssID in filteredSS[lsvid]:
			newFilter[lsvid][ssID] = filteredSS[lsvid][ssID]
	
	fout = open(csv_out, "w")
	tsvStr = ",".join(tsvOrder)
	fout.write("LSV_ID,Chr_Coord,{}\n".format(tsvStr))	
	for lsvid in newFilter:
		for ssID in newFilter[lsvid]:
			arrPSI = newFilter[lsvid][ssID]
			psistr = ",".join(map(str,arrPSI))
			fout.write("{},{},{}\n".format(lsvid,ssID,psistr))
	fout.close()
