#!/usr/bin/python

###############################################################
#	Author: Manishi Pandey
#	Command Line: python s03_merge_IRinfo.py -irPath ../output/IR_s02_output/Jan2_19/ -timepoints ../runs/sorted_bamFiles/timepoints.txt -outdir ../output/IR_s03_output
#	Date: January 3, 2019
###############################################################


from __future__ import division
import re, sys, os
from collections import defaultdict
from collections import OrderedDict
import argparse
import time

def parse_IRinfo_files(filename):
	IRinfo = defaultdict(list)
	with open(filename, "r")  as fIN:
		for line in fIN:
			line = line.rstrip("\n")
			data = line.split(",")
			IRinfo[data[0]] = data[1:]
	return IRinfo

def parse_IRfile(filename):
	irInfo = defaultdict(str)
	with open(filename, "r") as fIN:
		for line in fIN:
			line = line.rstrip("\n")
			data = line.split("\t")
			intronID = "{}:{}-{}".format(data[0],data[1],data[2])
			irInfo[intronID] = data[3]
	return(irInfo)

def parse_IR_coverage(filename):
	irCov = defaultdict(int)
	with open(filename, "r") as fIN:
		for line in fIN:
			line = line.rstrip("\n")
			data = line.split("\t")
			if (float(data[9]) >= 1.0):
				intronID = "{}:{}-{}".format(data[0],data[1],data[2])
				irCov[intronID] = int(data[6])
	return irCov
				

def read_timepoints(filename):
	timepts = OrderedDict()
	with open(filename, "r") as fIN:
		for line in fIN:
			line = line.rstrip("\n")
			data = line.split("\t")
			timepts[data[2]] = [ data[0],data[1] ]
	return timepts

def get_coord(sscoord):
	ssInfo = re.search(r'(.*):(.*)-(.*)', sscoord)
	ssArr = [int(ssInfo.group(i)) for i in range(2,4) ]
	return(ssArr)


if (__name__ == "__main__"):
	parser = argparse.ArgumentParser(prog="Get IR events for all important hits")
        parser.add_argument('-irPath', dest="irpath", help="Enter path to bed files that are mapped within a given intron")
	parser.add_argument('-covpath', dest="covpath", help="Enter the IR coverage path")
        parser.add_argument('-timepoints', dest="timepoints", help="Timepoints list file")
        parser.add_argument('-outdir', dest="outdir", help="File with readcounts at each timepoint")
        args = parser.parse_args()
	
	timepts_file = args.timepoints
	IRpath = args.irpath
	covpath = args.covpath
	outdir = args.outdir
	date = time.strftime("%y%m%d")

	IRfile = "/scratch/gslab/mpandey/lib/updated_IntronCoord.bed"
	geneIntrons = parse_IRfile(IRfile)
	timepts = read_timepoints(timepts_file)
	allVals = defaultdict(dict)
	pirInfo = defaultdict(list)
	for tp in timepts:
		srrIDs = timepts[tp]
		IRinfo = defaultdict(list)
		IRcov = defaultdict(list)
		for srrID in srrIDs:
			filename1 = "{}/{}.csv".format(IRpath,srrID)
			filename2 = "{}/{}.bed".format(covpath,srrID)
			#print(filename1)
			if (os.path.isfile(filename1)):
				tmpIR = parse_IRinfo_files(filename1)
				for iID in tmpIR:
					if (iID not in IRinfo):
						IRinfo[iID] = tmpIR[iID]
					else:
						IRinfo[iID].extend(tmpIR[iID])
			if (os.path.getsize(filename2) > 0):
				tmpcov = parse_IR_coverage(filename2)
				for iID in tmpcov:
					if (iID not in IRcov):
						IRcov[iID] = [ tmpcov[iID] ]
					else:
						IRcov[iID].append(tmpcov[iID])	
		count = 0; count1 = 0
		for intronID in sorted(IRinfo):
			tVals = map(int,IRinfo[intronID])
			coord = get_coord(intronID)
			intronLen = coord[1] - coord[0]
			IRcov_flag = 0
			if (intronID in IRcov):
				if (len(IRcov[intronID]) == 2):
					IRcov_flag = sum(IRcov[intronID])/2
					#print(IRcov_flag)
			if (tVals[0] > 0):
				count1 += 1
			avgIR = [ (tVals[i] + tVals[i+4])/2 for i in range(0,4) ]
			pir = 0.0001
			if (avgIR[0] >= 5):
				if (avgIR[1] >= 2) & (avgIR[2] >= 2) & (IRcov_flag >= 5):
					pir = (avgIR[1] + avgIR[2])/sum(avgIR[0:3])
					pir = round(pir, 4)
				else:
					pir = 0.01
				if (pir > 0.05):
					count += 1
			allVals[intronID][tp] = pir
		print("TP{},{},{}".format(tp,count,count1))
	#timep = timepts.keys()
	#print(timep)
	#timeStr = ",".join(timep)
	#if not(os.path.exists(outdir)):
	#	os.mkdir(outdir)
	#fout1 = open("{}/PIR_Summary_noDistinction.txt".format(outdir), "w")
	#fout2 = open("{}/{}_PIR_filtered_noDistinction.csv".format(outdir,date), "w")
	#fout1.write("intronID,GeneName,{}\n".format(timeStr))
	#fout2.write("intronID,GeneName,{}\n".format(timeStr))
	#for intronID in allVals:
	#	geneName = geneIntrons[intronID]
	#	#print(allVals[intronID][tp])
	#	pirArr = [ float(allVals[intronID][tp]) for tp in timep ]
	#	pirInfo[intronID] = pirArr
	#	if (all(pir < 0.05 for pir in pirArr)):
	#		next
	#	else:
	#		pirIndex = [ i for i,p in enumerate(pirArr) if (p != 0.0001) ]
	#		pMin = min([ pirArr[i] for i in pirIndex ])
	#		pMax = max([ pirArr[i] for i in pirIndex ])
	#		print(pirArr)
	#		print(pMax,pMin)
	#		if (pMin != pMax):
	#			#pirStd = [ round((p - pMin)/(pMax - pMin), 4) if (p != 0.0001) else "N" for p in pirArr ]
	#			pirStd = [ p if (p != 0.0001) else "N" for p in pirArr ]
	#			print(pirStd)
	#			noneInd = [ i for i,p in enumerate(pirStd) if (p == "N")]
	#			if (len(noneInd) < 0.75*len(pirStd)):
	#				pirStd_1 = [ 0.0001 if (i in noneInd) else p for i,p in enumerate(pirStd) ]
	#				pirStd_2 = [ 0.0001 if (float(p) == 0.0) else p for p in pirStd_1 ]
	#				pirStd_3 = [ 0.0001 if (float(p) == 0.01) else p for p in pirStd_2 ]
	#				print(pirStd_3)
	#				pirStd_str = ",".join(map(str,pirStd_3))
	#				fout2.write("{},{},{}\n".format(intronID,geneName,pirStd_str))
	#		pirStr = ",".join(map(str,pirArr))
	#		fout1.write("{},{},{}\n".format(intronID,geneName,pirStr))
	#}
