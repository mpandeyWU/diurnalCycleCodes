#!/usr/bin/python

###################################################
#	Author: Manishi Pandey
#	Date: December 5, 2018
#	Usage: python s02_3N_check_of_AS.py -infile ../output/altSS_s01_output/Dec6/1206_SSID_Qualified_Summary.csv -lsvtype ../output/altSS_s01_output/Dec6/1206_LSV_Type_Information.csv -outdir ../output/altSS_s02_output/Dec10/
#	Aim: Classify splicing events as 3N, Non-3N and ES. Further, classify ES events as 3N or Non-3N based on skipped exon length
#	Update: (Dec 6) Code updated to include ratio of GO and PFAM terms with new formatted output file
#	Update: (Dec 10) exon_skipping_check is updated to classify ES events as Alt_3N and Alt_Non3N based on added length of skipped exon
#	Update: (Dec 18) Added +1 to the length of Alt5 and Alt3
#####################################################

from __future__ import division
import re, sys, os
from collections import defaultdict
from collections import OrderedDict
import argparse
import numpy as np
import time

#Read qualified LSV file generated for all time-points
def parse_qualified_lsv(filename):
	lsvDict = defaultdict(dict)
	tagline = None
	with open(filename, "r") as fIN:
		for line in fIN:
			line = line.rstrip("\n")
			if not(line.startswith("LSV")):
				data = line.split(",")
				lsvDict[data[0]][data[1]] = data[2:]
			else:
				tagline = line
	return (lsvDict,tagline)

def get_coord(sscoord):
	ssInfo = re.search(r'(.*):(.*)-(.*)', sscoord)
	ssArr = [int(ssInfo.group(i)) for i in range(2,4) ]
	return(ssArr)

def parse_Binary_profile(filename):
	binary = defaultdict(list)
	with open(filename, "r") as fIN:
		for line in fIN:
			if not(line.startswith("LSV")):
				line = line.rstrip("\n")
				data = line.split(",")
				binary[data[0]] = data[1:]
	return binary

def threeN_check(ssIDs):
	threeN = 0
	count = 0
	ssID1 = ssIDs[0]
	ssID2 = ssIDs[1]
	ssInfo1 = get_coord(ssID1)
	ssInfo2 = get_coord(ssID2)	
	if (ssInfo1[0] == ssInfo2[0]) | (ssInfo1[1] == ssInfo2[1]):
		c1 = abs(int(ssInfo1[0]) - int(ssInfo2[0])) + 1
		c2 = abs(int(ssInfo1[1]) - int(ssInfo2[1])) + 1
		if (c1 == 1) & (c2%3 == 0):
			threeN = 1
		#	print(ssIDs,c1,c2)
			count += 1
		elif (c2 == 1) & (c1%3 == 0):
		#	print(ssIDs,c1,c2)
			threeN = 1
	print(count)
	return threeN

def parse_LSVtype_file(filename):
	lsvType = defaultdict(list)
	with open(filename, "r") as fIN:
		for line in fIN:
			if not(line.startswith("LSV")):
				binaryArr = [ 0, 0, 0]
				line = line.rstrip("\n")
				data = line.split(",")
				binaryArr = [ 0 if(x == "False") else 1 for x in data[1:] ]	
				lsvType[data[0]] = binaryArr
	return lsvType

def parse_gff_file(gff_file):
	exonInfo = defaultdict(dict)
	with open(gff_file, "r") as fIN:
		for line in fIN:
			if not(line.startswith("#")):
				line = line.rstrip("\n")
				data = line.split("\t")
				if (data[2] == "exon"):
					gene = re.search(r'ID=((.*).t(.*).v5.5.exon(.*));Parent=(.*)',data[8])
					geneName = gene.group(2)
					exonName = "{}:{}-{}".format(data[0],data[3],data[4])
					coord = map(int,[data[3], data[4]])
					exonInfo[geneName][exonName] = coord
	return exonInfo
				
#If either SS contains a known exon boundaries then it is marked as exon skipping event
def exon_skipping_check(ssIDs,exonDict):
	ES_flag = 0
	es3N = 0
	ES_SS = list()
	for ssid in ssIDs:
		ss = get_coord(ssid)
		exonLen = 0
		for exonid in exonDict:
			exon = exonDict[exonid]
			if (exon[0] > ss[0]) & (exon[1] < ss[1]):
				length = (exon[1] - exon[0] + 1)
				exonLen += length
				print(exon,ssid)
		ES_SS.append(exonLen)
	if (ES_SS[0] > 0):
		diff = 0
		if (len(ES_SS) > 1):
			diff = abs(ES_SS[1] - ES_SS[0])
		else:	
			diff = ES_SS[0]
		ES_flag = 1
		if (diff % 3 == 0):
			es3N = 1
			print("Frame preserving: {}".format(ssid))
		else:
			print("disrupting: {}".format(ssid))
		pass
	return (ES_flag, es3N)

def parse_pfam_GO_file(filename,lsvdict,outdir):
	fout = open("{}/GO_PFAM_frequency.txt".format(outdir), "w")
	fout.write("GO_PFAM_ID,Hits_Count,All_Count,Freq,Desc_1,Desc_2\n")
	pfamGO = defaultdict(dict)
	goDesc = defaultdict(list)
	goFreq = defaultdict(list)
	outDict = defaultdict(list)
	with open(filename, "r") as fIN:
		for line in fIN:
			line = line.rstrip("\n")
			if not(line.startswith("#")):
				data = line.split("\t")
				if (data[1].startswith("PF") | data[1].startswith("GO")):
					if (data[1] not in pfamGO):
						pfamGO[data[1]] = [ data[0] ]
					else:
						pfamGO[data[1]].append(data[0])
					desc = data[2].replace(",", ";")
					goDesc[data[1]] = [ data[2],data[4].rstrip('\r') ]
	#To get the total proteins associated with each GO term and PFAM ID
	for goid in pfamGO:
		geneList = list(set(pfamGO[goid]))
		count = len(geneList)
		goDesc[goid].append(count)
	for lsvid in lsvdict:
		gene = re.search(r'(.*).v5.5:(.*)',lsvid)
		geneName = gene.group(1)
		for goid in pfamGO:
			if (geneName in pfamGO[goid]):
				goFreq[goid].append(geneName)
				pass;
	for goid in goFreq:
		tLen = len(list(set(goFreq[goid])))
		allCount = goDesc[goid][2]
		if (tLen >= 3):
			freq = tLen/allCount
			outDict[goid] = [ tLen, allCount, freq ]
	for goid in sorted(outDict, key=lambda x:outDict[x][2], reverse = True):
		countStr = ",".join(list(map(str,outDict[goid])))
		if (goid.startswith("GO")):
			fout.write("{},{},{}\n".format(goid,countStr,goDesc[goid][0]))	
		else:
			fout.write("{},{},{}\n".format(goid,countStr,goDesc[goid][1]))	
	return (pfamGO, goDesc)

def intron_retention_check(ssIDs):
	irFlag = 0
	#print(ssIDs)
	for ssid in ssIDs:
		ssInfo = get_coord(ssid)
		diff = ssInfo[1] - ssInfo[0]
		#print(diff)
		if (diff == 1):
			#print(ssid)
			irFlag += 1
	return irFlag

def get_binaryAnnot():
	bdict = { "1,0,0" : "Light_G1",
		"0,1,0" : "SM",
		"0,0,1" : "Dark_G1",
		"1,1,0" : "Light_SM",
		"1,0,1" : "Light_Dark",
		"0,1,1" : "SM_Dark",
		"1,1,1" : "All",
		"0,0,0" : "N/A" }
	
	return bdict

if (__name__ == "__main__"):
	parser = argparse.ArgumentParser(prog = "Identify 3N vs Non-3N sites")
	parser.add_argument('-infile', dest='infile', help="Input file with qualified LSV")
	parser.add_argument('-lsvtype', dest='lsvtype', help="Input file with LSV type from MAJIQ")
	parser.add_argument('-outdir', dest='outdir', help="Input file with LSV type from MAJIQ")
	args = parser.parse_args()

	infile = args.infile
	lsvtypefile = args.lsvtype
	outdir = args.outdir
	date = time.strftime("%m%d")
	
	if not(os.path.isdir(outdir)):
		os.makedirs(outdir)
	
	ssAnnot = defaultdict(list)	
	bAnnot = get_binaryAnnot()
	lsvdict,tagline = parse_qualified_lsv(infile)
	
	binaryfile = "../output/altSS_s01_output/Dec6/1206_Group_based_binaryLSV.csv"
	binarydict = parse_Binary_profile(binaryfile)	
	
	lsvType = parse_LSVtype_file(lsvtypefile)
	
	gff_file = "/scratch/gslab/mpandey/lib/PhytozomeV12/Creinhardtii/annotation/Creinhardtii_281_v5.5.gene_exons.gff3"
	exonInfo = parse_gff_file(gff_file)

	pfam_file = "/scratch/gslab/mpandey/lib/GO_terms_Chlamydomonas.tsv"
	pfamGO, goDesc = parse_pfam_GO_file(pfam_file,lsvdict,outdir)

	ls = list(set(lsvdict.keys()) - set(lsvType.keys()))
	IRlsv = list()
	count = 0
	#print(len(lsvdict.keys()))
	#lsvdict have 5 IDs that do not exist in binary and lsvType because their RC == 2 reads qualification do not
	#align at the same time-point. So, they are filtered out in s01 program when creating binary file.
	for lsvid in lsvType: #lsvid are taken from lsvType file 
		#print(lsvid)
		gene = re.search(r'(.*).v5.5:(.*)',lsvid)
		geneName = gene.group(1)
		lsvT = lsvType[lsvid]
		#print(lsvT)
		ssIDs = lsvdict[lsvid].keys()
		irFlag = intron_retention_check(ssIDs)
		if (irFlag):
			IRlsv.append(lsvid)
		if (len(ssIDs) == 2):
			esFlag,flag_3N = exon_skipping_check(ssIDs,exonInfo[geneName])
			if (esFlag == 1):
				#ssAnnot['ES'].append(lsvid)
				if (flag_3N == 1):
					ssAnnot["Alt_3N"].append(lsvid)
				else:
					ssAnnot["Alt_Non3N"].append(lsvid)
				count += 1
			else:
				N3 = threeN_check(ssIDs)
				if (N3 == 1):
					ssAnnot["Alt_3N"].append(lsvid)
				else:	
					ssAnnot["Alt_Non3N"].append(lsvid)
		else:
			ssAnnot["Complex"].append(lsvid)
	for entry in ssAnnot:
		print(entry, len(ssAnnot[entry]))
	fout = open("{}/{}_Summary_All_LSV.csv".format(outdir,date), "w")
	fout.write("LSV_ID,Gene_ID,LSV_Type,Canonical_PSI,Alternate_PSI,Light_G1,S-M,Dark_G1,GO_Terms,Pfam_Terms\n")
	altSites = defaultdict(dict)
	for entry in ssAnnot:
		altN = ssAnnot[entry]
		for lsvid in altN:
			psiVal = list()
			gene = re.search(r'(.*).v5.5:(.*)', lsvid)
			geneName = gene.group(1)
			bdata = [ 0, 0, 0]
			if (lsvid in binarydict):	
				bdata = binarydict[lsvid]
			psiInfo = lsvdict[lsvid]
			ssName = sorted(psiInfo.keys())
			ssPSI = list()
			if (entry == "Alt_3N") | (entry == "Alt_Non3N"):# | (entry == "ES"):
				for ssid in ssName:
					#repPSI = [psi for psi in psiInfo[ssid] if (psi > '0.0001')]
					repPSI = [psi for psi in psiInfo[ssid] if (float(psi) > 0.05)]
					repPSI = map(float,repPSI)
					avgPSI = sum(repPSI)/len(repPSI)
					ssPSI.append(avgPSI)
				#canPSI = max(ssPSI)
				#altPSI = min(ssPSI)
				#psiVal = [ canPSI, altPSI ]
				canIndex, altIndex = 0,0
				if (ssPSI[0] > ssPSI[1]):
					canIndex = 0; altIndex = 1
				else:
					canIndex = 1; altIndex = 0
				psiVal = [ssPSI[canIndex], ssPSI[altIndex] ]
				ssA = [ ssName[canIndex], ssName[altIndex] ]
				altSites[lsvid][ssA[1]] = lsvdict[lsvid][ssA[1]]
				altSites[lsvid][ssA[1]].append(entry)
			else:
				psiVal = [ 'N', 'N' ]
			psiStr = ",".join(list(map(str,psiVal)))
			binaryStr = ",".join(map(str,bdata))
			fout.write("{},{},{},{},{}\n".format(lsvid,geneName,entry,psiStr,binaryStr))
				
	altout = open("{}/{}_AltSS_Only_Summary.csv".format(outdir,date), "w")
	altout.write("{},AltType\n".format(tagline))
	for lsvid in altSites:
		for ssid in altSites[lsvid]:
			altStr = ",".join(list(map(str,altSites[lsvid][ssid])))
			altout.write("{},{},{}\n".format(lsvid,ssid,altStr))
			
	#print(lsvid,entry,goOnt,pfamN,bdata,psiVal)
		
