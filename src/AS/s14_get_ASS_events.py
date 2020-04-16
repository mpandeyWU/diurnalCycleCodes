#!/usr/bin/python

#python s14_get_ASS_events.py -tsvpath ../output/Majiq_PSI_All/tsvFiles/ -rcfile ../output/181119_ReadCount_In_TP_format_Min2.csv
import re, sys, os
from collections import defaultdict
import argparse
import numpy as np

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
				#Index of entries that have PSI greater than 0.05
				index = [ i for i,k in enumerate(psi) if (k > 0.05)]
				if (len(index) > 1):
					for i in index:
						SSt = "{}:{}".format(chrom,sscoord[i])
						psi[i] = round(float(psi[i]), 4)
						#print(data[2],SSt,psi[i])
						ssDict[data[2]][SSt] = psi[i]
					ssType[data[2]] = [ data[6],data[7],data[8] ] 
	
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
					#if the ssID is in the rcDict and have read count of at least 2; then it is processed further
					if ((ssID in rcDict) & (rcDict[ssID] >= 2)):
						if (ssID in allLSV[lsvid]):
							allLSV[lsvid][ssID].append(ssDict[lsvid][ssID])
						else:
							l = ['N'] * (k+1)
							l[k] = ssDict[lsvid][ssID]
							allLSV[lsvid][ssID] = l
			else:
				allLSV[lsvid] = dict()
				lsvType[lsvid] = ssType[lsvid]
				for ssID in ssDict[lsvid]:
					if ((ssID in rcDict) & (rcDict[ssID] >= 2)):
						l = ['N'] * (k+1)
						l[k] = ssDict[lsvid][ssID]
						allLSV[lsvid][ssID] = l
			#print(lsvid,flag,len(ssDict[lsvid]))
		#All entries that are not filled with PSI value are filled with 'N'
		otherIDs = set(allLSV.keys()) - set(ssDict.keys())
		for oID in otherIDs:
			for ssID in allLSV[oID]:
				allLSV[oID][ssID].append('N')
	count = 0
	summary = defaultdict(list)
	binaryLSV = defaultdict(list)
	fout1 = open("../output/s14_output/181203_Qualified_LSV_Summary.csv", "w")
	fout2 = open("../output/s14_output/181203_Binary_LSV.csv", "w")
	fout3 = open("../output/s14_output/181203_summary_across_TP.csv", "w")
	fout4 = open("../output/s14_output/181203_Group_based_binaryLSV.csv", "w")
	fout5 = open("../output/s14_output/181203_LSV_Frequency.csv", "w")
	tpList = ",".join(tsvOrder)
	fout1.write("LSV_ID,SS_ID,{}\n".format(tpList))
	fout2.write("LSV_ID,{}\n".format(tpList))
	allTP = list(); qualLSV = list()
	count = 0
	for lsvid in allLSV:
		if(len(allLSV[lsvid]) > 1):
			qualLSV.append(lsvid)
			print(lsvid,lsvType[lsvid])
			#Intialize a numpy matrix of zeros with rows as all LSV ids and time-points as columns
			ssOccur = np.zeros((len(allLSV[lsvid]),len(tsvOrder)))
			for k,ssid in enumerate(allLSV[lsvid].keys()):
				ssArr = allLSV[lsvid][ssid]
				if (len(ssArr) == 27):
					count += 1
				#Get the indices/time-points that are not 'N' and have some value
				index = [i for i,val in enumerate(ssArr) if (val != 'N')]
				#Change the value from 0 to 1 for those indices
				for i in index:
					ssOccur[k][i] = 1
				ssLine = ",".join(map(str,ssArr))
				#Write the PSI values for each LSV id in file_1
				fout1.write("{},{},{}\n".format(lsvid,ssid,ssLine))
			#For each LSV, sum the matrix across column for each time-point
			tp = np.sum(ssOccur, axis = 0)
			#If there are more than one SS existing for a time-point with defined PSI, mark that TP value as 1 else 0
			tpN = [1 if i > 1 else 0 for i in tp]
			allTP.append(tpN)
			tpNstr = ",".join(map(str,tpN))
			fout2.write("{},{}\n".format(lsvid,tpNstr))
			summary[lsvid] = tpN
	print("ES temp: {}".format(count))
	allTP = np.asarray(allTP)
	#Count for each time-point
	tpCount = np.sum(allTP, axis = 0)
	tpDict = dict(zip(tsvOrder,tpCount))
	for tp in tsvOrder:
		fout3.write("{},{}\n".format(tp,tpDict[tp]))
	#Frequency of occurrence of each LSV
	ssCount = np.sum(allTP, axis = 1)
	ssFD = dict(zip(qualLSV,ssCount))
	for lsvid in qualLSV:
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
	print(count)	
	
	count = 0
	fout4.write("LSV_ID,G1_phase,S-M_phase,Post_SM_phase\n")
	for lsvid in ssFreq:
		lCount = [i for i,val in enumerate(ssFreq[lsvid]) if (val > 0)]
		if (len(lCount) == len(ssFreq[lsvid])):
			count += 1
		pr = ",".join(map(str,ssFreq[lsvid]))
		fout4.write("{},{}\n".format(lsvid,pr))
	
	print(count)	

	
