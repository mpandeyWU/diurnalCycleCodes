#!/usr/bin/python

#########################################################
#	Purpose: Parse tsv output files and filter the relevant LSVs
#	Author: Manishi Pandey
#	Date: October 16, 2018
#	Usage: python s05_parse_PSI_vals.py -majiqpath ../output/Majiq_PSI_All/tsvFiles/ -outfile ../output/s05_output/all_LSV_list.txt
#	Update: To include variance and coefficient of variation in the analysis and then sorting the SS by these values
###########################################################

import re, sys, os
import argparse
from collections import defaultdict
import numpy as np
import scipy.stats as sp

def parse_majiq_tsv(filename):
	ssDict = defaultdict(float)
	ssNames = defaultdict(list)
	with open(filename, "r") as fIN:
		for line in fIN:
			if not(line.startswith("#")):
				line = line.rstrip("\n")
				data = line.split("\t")
				psi = data[3].split(";")
				sscoord = data[14].split(";")
				chrom = data[12]; sign = data[13]
				for i,val in enumerate(sscoord):
					SSt = "{}:{}:{}".format(chrom,val,sign)
					ssID = "{}@{}".format(SSt,data[2]) #data[2] is LSV-ID
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

	ssName = defaultdict(list)

	for k,sampleID in enumerate(tsvOrder):
		filename = "{}.psi.tsv".format(sampleID)
		print(filename)
		#fname = re.search(r'(.*)psi.tsv', filename)
		#sampleID = fname.group(1)
		#print(fname.group(1),filename)
		filepath = "{}/{}".format(tsv_path,filename)
		ssDict = parse_majiq_tsv(filepath)
		for ssid in ssDict:
			chrID = ssid.split("@")
			if (chrID[1] in ssName):
				if (chrID[0] not in ssName[chrID[1]]):
					ssName[chrID[1]].append(chrID[0])
			else:
				ssName[chrID[1]] = [ chrID[0] ]
			if (ssid in allSS):
				allSS[ssid].append(ssDict[ssid]) #allSS dictionary has summary of all ssID with lsvid with a list of psi-values for each time-point.
			else:
				l = ['NaN']*(k+1)
				l[k] = ssDict[ssid]
				allSS[ssid] = l
		otherIDs = set(allSS.keys()) - set(ssDict.keys())
		for oID in otherIDs:
			allSS[oID].append('NaN')
	
	f2 = open("../output/s05_output/LSVIDs_with_SS.csv", "w")
	for lsvid in ssName:
		ssids = ";".join(ssName[lsvid])
		f2.write("{}\t{}\n".format(lsvid, ssids))
	f2.close()

	varArr = list()
	cvArr = list()
	filteredSS = defaultdict(list)
	fcount = 0; pCount = 0;
	for sampleID in allSS:
		arrPSI = [ x for x in allSS[sampleID] if (x != 'NaN') ]
		if (len(arrPSI) >= 2):
			relPSI = [i for i,x in enumerate(arrPSI) if (float(x) >= 0.05)]
			if (len(relPSI) > 0):
				varPSI = np.var(arrPSI)
				cvPSI = sp.variation(arrPSI)
				varArr.append(varPSI)
				cvArr.append(cvPSI)
				allSS[sampleID].extend([varPSI, cvPSI])
				filteredSS[sampleID] = allSS[sampleID]
				#print(arrPSI,varPSI)	
				#if (cvPSI >= 0.5):	#Need to choose a proper cut-off
				#	ids = sampleID.split("@")
				#	filteredSS[ids[1]][ids[0]] = allSS[sampleID]
			else:
				fcount += 1
		else:
			pCount += 1
	
	print("Entries with less than .05 PSI: {}\nEntried with few occurrence: {}".format(fcount, pCount))
	newFilter = defaultdict(dict)
	hist = np.histogram(varArr)
	med = np.median(varArr)
	mean = np.mean(varArr)
	print(hist)
	print(med, mean, min(varArr), max(varArr))
	med1 = np.median(cvArr)
	mean1 = np.mean(cvArr)
	cvHist = np.histogram(cvArr)
	print(cvHist)
	print(med1, mean1, min(cvArr), max(cvArr))
	allSSList = allSS.keys()
	
	total = len(filteredSS.keys())
	fout = open(csv_out, "w")
	tsvStr = ",".join(tsvOrder)
	fout.write("LSV_ID,Chr_Coord,{},Variance,CV\n".format(tsvStr))
	fiveCount = 0
	cutoff = 0.25 * total
	print(cutoff, total)	
	for lsvid in sorted(filteredSS, key=lambda x: filteredSS[x][28], reverse=True):
		#if (filteredSS[lsvid][28] > mean1):
		if (fiveCount <= cutoff):
			ids = lsvid.split("@")
			arrPSI = filteredSS[lsvid]
			NaNs = [i for i,x in enumerate(arrPSI) if (x == "NaN")]
			for index in NaNs:
				arrPSI[index] = 0.0
			psistr = ",".join(map(str,arrPSI))
			fout.write("{},{},{}\n".format(ids[1],ids[0],psistr))
			#print(lsvid, filteredSS[lsvid], len(filteredSS[lsvid]))
			fiveCount += 1
		else:
			next
			#newFilter[lsvid][ssID] = filteredSS[lsvid][ssID]
	fout.close()

	#for lsvid in sorted(filteredSS, key=lambda x: filteredSS[x][28], reverse=True):
		#otherSS = [ x for x in allSSList if (re.search(r'(.*)@' + re.escape(lsvid), x))] #Syntax to combine variable content with regex
		#print(otherSS)
		#for oID in otherSS:
		#	sID = oID.split("@")
		#	newFilter[lsvid][sID[0]] = allSS[oID]
	#		print(lsvid, filteredSS[lsvid], len(filteredSS[lsvid]))
			#newFilter[lsvid][ssID] = filteredSS[lsvid][ssID]
	

