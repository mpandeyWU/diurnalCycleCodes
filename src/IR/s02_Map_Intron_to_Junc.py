#!/usr/bin/python

import re, sys, os
from collections import defaultdict
from collections import OrderedDict
import argparse
import subprocess

#Function to return coordinates as integers from the standard format
def get_coord(sscoord):
	ssInfo = re.search(r'(.*):(.*)-(.*)', sscoord)
	ssArr = [int(ssInfo.group(i)) for i in range(2,4) ]
	#ssArr = [ ssInfo.group(1) ] + ssArr
	return(ssArr)

#Function to parse CIGAR and get the mapped length
def get_cigar_len(newcig):
	cigLen = 0
	while(len(newcig) > 0):
		if (re.match(r'(\d+)([MDP\=X])(.*)', newcig)):
			cigLen = cigLen + int(re.match(r'(\d+)([MDP\=X])(.*)', newcig).group(1))
			newcig = re.match(r'(\d+)([MDP\=X])(.*)', newcig).group(3)
		elif (re.match(r'(\d+)([NSIH])(.*)', newcig)):
			newcig = (re.match(r'(\d+)([NSIH])(.*)', newcig)).group(3)
		else:
			#print("Failed to parse CIGAR {}\n".format(newcig))
			newcig = ""
	return cigLen

#Function to parse CIGAR IDs with N and get the splice junction coordinate
def get_SJ_from_CIGAR(newcig,chrom,start):
	cigLen = 0; junctions = list();
	while(len(newcig) > 0):
		if (re.match(r'(\d+)([MDP\=X])(.*)', newcig)):
			cigLen = cigLen + int(re.match(r'(\d+)([MDP\=X])(.*)', newcig).group(1))
			newcig = re.match(r'(\d+)([MDP\=X])(.*)', newcig).group(3)
		elif (re.match(r'(\d+)([N])(.*)', newcig)):
			juncStart = start + cigLen
			cigLen = cigLen + int(re.match(r'(\d+)([N])(.*)', newcig).group(1))
			newcig = re.match(r'(\d+)([N])(.*)', newcig).group(3)
			juncEnd = start + cigLen
			juncID = "{}:{}-{}".format(chrom,juncStart,juncEnd)
			junctions.append(juncID)
			start = juncEnd
			cigLen = 0
		elif (re.match(r'(\d+)([SIH])(.*)', newcig)):
			newcig = (re.match(r'(\d+)([NSIH])(.*)', newcig)).group(3)
		else:
			#print("Failed to parse CIGAR {}\n".format(newcig))
			newcig = ""
	return junctions

#Function to map 
def map_junctions_to_intron(intronFile,intronJuncFile):
	intronJunc = defaultdict(list)
	inCoord = defaultdict(dict)
	with open(intronFile, "r") as fIN:
		for line in fIN:
			line = line.rstrip("\n")
			data = line.split("\t")
			chrID = "{}:{}-{}".format(data[0],data[1],data[2])
			inCoord[data[3]][chrID] = [ data[0],data[1],data[2] ]
	fIN.close()
	with open(intronJuncFile, "r") as fIN1:
		for line in fIN1:
			line = line.rstrip("\n")
			data = line.split("\t")
			juncID = "{}:{}-{}".format(data[0],data[1],data[2])
			for intronID in inCoord[data[3]]:
				ssInfo = inCoord[data[3]][intronID]
				#print(ssInfo,data)
				if ((data[1] == ssInfo[1]) | (data[1] == ssInfo[2])):
					#print("Qualified_IDs: {}\t{}".format(intronID,juncID))
					if (intronID in intronJunc):
						intronJunc[intronID].append(juncID)
					else: 
						intronJunc[intronID] = [ juncID ]
					next
	fIN1.close()
	return(intronJunc)

def reads_within_intron(filename):
	intronInfo = defaultdict(dict)
	intronCount = defaultdict(int)
	with open(filename, "r") as fIN:
		for line in fIN:
			line = line.rstrip("\n")
			data = line.split("\t")
			inID = "{}:{}-{}".format(data[-7],data[-6],data[-5])
			mapLoc = "{}:{}-{}".format(data[0],data[1],data[2])
			intronInfo[inID][mapLoc] = [ data[3],data[4],data[7] ]
	for intronID in intronInfo:
		intronCount[intronID] = 0
		for readID in intronInfo[intronID]:
			readInfo = intronInfo[intronID][readID]
			#print(readInfo)
			if (int(readInfo[1]) >= 20):
				mapLen = get_cigar_len(readInfo[2])
				if (mapLen > 50):
					intronCount[intronID] += 1
	return(intronCount)

def get_CIGAR_for_IR(newcig):
	cigMap = list(); cigIR = list()
	while(len(newcig) > 0):
		if (re.match(r'(\d+)([MP\=X])(.*)', newcig)):
			cigMap.append(int(re.match(r'(\d+)([MP\=X])(.*)', newcig).group(1)))
			newcig = re.match(r'(\d+)([MP\=X])(.*)', newcig).group(3)
		elif (re.match(r'(\d+)([N])(.*)', newcig)):
			cigIR.append(int(re.match(r'(\d+)([N])(.*)', newcig).group(1)))
			newcig = re.match(r'(\d+)([N])(.*)', newcig).group(3)
		elif (re.match(r'(\d+)([DSIH])(.*)', newcig)):
			newcig = (re.match(r'(\d+)([DSIH])(.*)', newcig)).group(3)
		else:
			#print("Failed to parse CIGAR {}\n".format(newcig))
			newcig = ""
	return (cigMap, cigIR)

def reads_on_intron_junctions(filename,intronJunc):
	splicedReads = defaultdict(int)
	IEjuncReads = defaultdict(int)
	count = 0
	with open(filename, "r") as fIN:
		for line in fIN:
			line = line.rstrip("\n")
			data = line.split("\t")
			if (int(data[10]) >= 10):
				if (re.search(r'(.*)M(\d+)N(.*)M', data[13])):
					count += 1
					junctions = get_SJ_from_CIGAR(data[13],data[6],int(data[7]))
					for juncID in junctions:
							splicedReads[juncID] += 1
				else:	
					mapLen = get_cigar_len(data[13])	
					if(mapLen > 50):
						convIE = "{}:{}-{}".format(data[0],int(data[1])+10,int(data[2])-10)
						IEjuncReads[convIE] += 1
	print(count)
	print(len(splicedReads.keys()))
	print(len(IEjuncReads.keys()))
	return(splicedReads, IEjuncReads)
	
if (__name__ == "__main__"):
        parser = argparse.ArgumentParser(prog="Get IR events for all important hits")
        parser.add_argument('-withinIRpath', dest="withinpath", help="Enter path to bed files that are mapped within a given intron")
        parser.add_argument('-juncIRpath', dest="juncpath", help="Enter path to bed files that are mapped within a given intron")
        parser.add_argument('-timepoints', dest="timepoints", help="Timepoints list file")
        parser.add_argument('-outdir', dest="outdir", help="File with readcounts at each timepoint")
        args = parser.parse_args()
	
	withinpath = args.withinpath
	juncpath = args.juncpath
	timepts_file = args.timepoints
	outdir = args.outdir

	#Intron Coordinate file and Intron-Exon junction coordinate file
	intronFile = "/scratch/gslab/mpandey/lib/updated_IntronCoord.bed"
	intronJuncFile = "../../output/IR_s01_output/a04_Chlamy_exonsIntronJunctions.bed"
	
	#
	fIN = open(timepts_file, "r")
	timepts = [ line.rstrip("\n") for line in fIN.readlines() ]
	print(timepts)
	intronJunc = map_junctions_to_intron(intronFile,intronJuncFile)
	
	for srrID in timepts:
		withinFile = "{}/{}.bed".format(withinpath,srrID)
		intronCount = reads_within_intron(withinFile)
		print(srrID,len(intronCount.keys()))
		junctionFile = "{}/{}.bed".format(juncpath,srrID)
		splicedReads, IEreads = reads_on_intron_junctions(junctionFile,intronJunc)
		outFile = "{}/{}.csv".format(outdir,srrID)
		fout = open(outFile, "w")
		for intronID in intronJunc:
			IEjunc = intronJunc[intronID]
			intronInfo = [0] * 4
			if (intronID in splicedReads):
				intronInfo[0] = splicedReads[intronID]
			if (intronID in intronCount):
				intronInfo[3] = intronCount[intronID]
			for i,IEid in enumerate(IEjunc):
				if (IEid in IEreads):
					#print(IEid,IEreads[IEid])
					intronInfo[i+1] = IEreads[IEid]
			infoStr = ",".join(map(str,intronInfo))
			fout.write("{},{}\n".format(intronID,infoStr))
		#Account for novel junctions and Intron reads	
		#	print(inID,intronCount[inID])
		#reads_on_intron_junctions(juncmapFile)
