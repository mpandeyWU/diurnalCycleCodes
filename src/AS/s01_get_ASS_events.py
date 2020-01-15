#!/usr/bin/python

#################################################
#	Purpose: Summarize AS events from all timepoints
##################################################

#python s14_get_ASS_events.py -tsvpath ../output/Majiq_PSI_All/tsvFiles/ -rcfile ../output/181119_ReadCount_In_TP_format_Min2.csv
import re, sys, os
from collections import defaultdict
import argparse
import numpy as np
import time

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

def intron_retention_check(ssIDs):
	IR_flag = [ 0 ]
	coord = list()
	for ssid in ssIDs:
		ssA = get_coord(ssid)
		coord.append(ssA)
	diff = [(int(coord[i][1]) - int(coord[i][0])) for i in range(0,len(coord))]
	ir = [i for i,val in enumerate(diff) if (val == 1)]
	if (ir):
		IR_flag[0] = 1
		for i in ir:
			IR_flag.append(ssIDs[i])
	return IR_flag #First value is 0 or 1 based on if event is IR and second value is IR coordinate

def parse_majiq_tsv(filename):
	ssDict = defaultdict(dict)
	ssType = defaultdict(list)
	ssNames = defaultdict(list)
	with open(filename, "r") as fIN:
		for line in fIN:
			if not(line.startswith("#")):
				line = line.rstrip("\n")
				data = line.split("\t")
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
		flag = 0
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
					#Since Intron Retention events are not present in rcDict file, 
					#I verified them based on coordinates and included them in the analysis
					ssInfo = get_coord(ssid)
					if ((ssInfo[1] - ssInfo[0]) == 1):
						ssIndex[ssid] = nk
		if (flag >= 2):
			for ssid in ssIndex:
				filteredIndex[lsvid][ssid] = ssIndex[ssid] #index is the list of time-point index that qualify all filters: Used for binaryLSV
			count += 1
	print(count)
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
			#print(lsvid,ssid,ssLen,allLSV[lsvid][ssid])
			if (ssLen < 27):
				count += 1
	print(count)
	filteredIndex = apply_filters(allLSV, rcDict)
	
	count = 0
	#print(filteredIndex)
	for lsvid in filteredIndex:
		if (lsvid.startswith("Cre17.g745097")):
		#print(len(filteredIndex[lsvid].keys()))
		#if (len(filteredIndex[lsvid].keys()) > 1):
			for ssid in filteredIndex[lsvid]:
				print(lsvid,ssid,filteredIndex[lsvid][ssid])
	print(count)

	summary = defaultdict(list)
	binaryLSV = defaultdict(list)
	tpList = ",".join(tsvOrder)
	if not(os.path.isdir(outdir)):
		os.makedirs(outdir)	
	fout1 = open("{}/{}_Qualified_LSV_Summary.csv".format(outdir,date), "w")
	fout6 = open("{}/{}_SSID_Qualified_Summary.csv".format(outdir,date), "w")
	fout7 = open("{}/{}_LSV_Type_Information.csv".format(outdir,date), "w")
	fout2 = open("{}/{}_Binary_LSV.csv".format(outdir,date), "w")
	fout3 = open("{}/{}_summary_across_TP.csv".format(outdir,date), "w")
	fout4 = open("{}/{}_Group_based_binaryLSV.csv".format(outdir,date), "w")
	fout5 = open("{}/{}_LSV_Frequency.csv".format(outdir,date), "w")
	fout1.write("LSV_ID,SS_ID,{}\n".format(tpList))
	fout6.write("LSV_ID,SS_ID,{}\n".format(tpList))
	fout2.write("LSV_ID,{}\n".format(tpList))
	fout7.write("LSV_ID,Alt_5SS,Alt_3SS,ES\n")
	allTP = list(); qualLSV = list()
	count = 0
	for lsvid in allLSV:
		if(lsvid in filteredIndex):
			for ssid in allLSV[lsvid]:	
				ssArr = allLSV[lsvid][ssid]
				ssLine = ",".join(map(str,ssArr))
				#Write the PSI values for each LSV id in file_1
				fout1.write("{},{},{}\n".format(lsvid,ssid,ssLine))	
			#Intialize a numpy matrix of zeros with rows as all LSV ids and time-points as columns
			#Many AS events have >= 3 SS with one alternate site with PSI < 0.05; all SS with less than 0.05 PSI throughout, will not be included in binary 
			ssList = sorted(allLSV[lsvid].keys())
			#print(ssList)
			ssOccur = np.zeros((len(ssList),len(tsvOrder)))
			for k,ssid in enumerate(ssList):
				if (ssid in filteredIndex[lsvid]):		
					qualIndex = filteredIndex[lsvid][ssid]
					#if (13 in qualIndex):
					#	print(lsvid,ssid,qualIndex)
					ssArr = allLSV[lsvid][ssid]
					ssLine = ",".join(map(str,ssArr))
					fout6.write("{},{},{}\n".format(lsvid,ssid,ssLine))	
					#if (lsvid == "Cre03.g190311.v5.5:t:5950983-5954206"):
					#	print(lsvid,ssid,qualIndex)
					#Change the value from 0 to 1 for those indices
					for i in qualIndex:
						ssOccur[k][i] = 1
			if (lsvid == "Cre03.g179820.v5.5:t:4921423-4921559"):
				print(ssOccur)
			#For each LSV, sum the matrix across column for each time-point
			tp = np.sum(ssOccur, axis = 0)
			#If there are more than one SS existing for a time-point with defined PSI, mark that TP value as 1 else 0
			tpN = [1 if i > 1 else 0 for i in tp]
			#print(tpN[14])
			if (lsvid == "Cre03.g179820.v5.5:t:4921423-4921559"):
				print(tp, tpN)
			if (all(x == 0 for x in tpN)):
				continue
			else:
				allTP.append(tpN)
				tpNstr = ",".join(map(str,tpN))
				lsvtype_str = ",".join(lsvType[lsvid])
				fout2.write("{},{},{}\n".format(lsvid,tpNstr,lsvtype_str))
				fout7.write("{},{}\n".format(lsvid,lsvtype_str))
				summary[lsvid] = tpN
	allTP = np.asarray(allTP)
	#Count for each time-point
	tpCount = np.sum(allTP, axis = 0)
	tpDict = dict(zip(tsvOrder,tpCount))
	for tp in tsvOrder:
		print(tp,tpDict[tp])
		fout3.write("{},{}\n".format(tp,tpDict[tp]))
	#Frequency of occurrence of each LSV
	ssCount = np.sum(allTP, axis = 1)
	ssFD = dict(zip(filteredIndex.keys(),ssCount))
	for lsvid in ssFD:
		fout5.write("{},{}\n".format(lsvid,ssFD[lsvid]))

	G1 = range(0,10)
	G2 = range(10,18)
	G3 = range(18,26)
	tp1 = [ tsvOrder[i] for i in G1 ]
	tp2 = [ tsvOrder[i] for i in G2 ]
	tp3 = [ tsvOrder[i] for i in G3 ]
	print(tp1,tp2,tp3)
	ssFreq = defaultdict(list)
	count = 0
	for lsvid in summary:
		binaryArr = [0] * 3
		if (lsvid == "Cre08.g360450.v5.5:s:723475-723570"):
			print(summary[lsvid])
		phase1 = [summary[lsvid][i] for i in G1]
		phase2 = [summary[lsvid][i] for i in G2]
		phase3 = [summary[lsvid][i] for i in G3]
		fp1 = sum(phase1); fp2 = sum(phase2); fp3 = sum(phase3)
		if (fp1 > 0):
			binaryArr[0] = 1
		if (fp2):
			binaryArr[1] = 1
		if (fp3):
			binaryArr[2] = 1 
		ssFreq[lsvid] = binaryArr
		if (all(b == 0 for b in binaryArr)):
			print(lsvid,binaryArr)
			count += 1
	print("Zero hits: {}".format(count))	
	
	count = 0
	fout4.write("LSV_ID,Light_G1,S-M,Dark_G1\n")
	for lsvid in ssFreq:
		lCount = [i for i,val in enumerate(ssFreq[lsvid]) if (val > 0)]
		if (len(lCount) == len(ssFreq[lsvid])):
			count += 1
		pr = ",".join(map(str,ssFreq[lsvid]))
		fout4.write("{},{}\n".format(lsvid,pr))
	
	print(count)	

	
