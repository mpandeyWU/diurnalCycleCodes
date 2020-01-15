#!/usr/bin/python

# Filter Criteria:
#1. SS should have readcount >= 2.
#2. SS should not be in any repeat region.
#3. SS should occur in more than one sample.

##Run the following command to generate the genome utilizing the output of this code:
#srun --mem=30G -c 4 /opt/apps/star/2.5.2b/bin/STAR --runMode genomeGenerate --runThreadN 4 --genomeDir ../lib/two-pass_genome/ --genomeFastaFiles /scratch/gslab/mpandey/lib/PhytozomeV12/Creinhardtii/assembly/Mod_Creinhardtii_281.fa --sjdbFileChrStartEnd ../lib/SJ_files/0904-SJs_for_STAR_annotation.tab 

import re, sys, os
import time
from collections import defaultdict
from collections import OrderedDict
import subprocess
import argparse

def run(cmd):
	p = subprocess.Popen(cmd, shell='True', stdout=subprocess.PIPE)
	p.wait()
	out, errs = p.communicate()
	return (out, errs)

#Function to read SJ.out.tab and store the data in the following dictionary:
#ssDict[SS_ID] = [ Array of read count associated with SS in each dataset ]

def read_SJout_file(filename):
	ssDict = defaultdict(int)
	with open(filename, "r") as fIN:
		for line in fIN:
			data = line.split("\t")
			ssID = "{}:{}-{}-{}".format(data[0],data[1],data[2],data[3])
			readCount = int(data[6])
			ssDict[ssID] = readCount
	return ssDict
		
#Mark all SS that occur in repeat regions and filter them out.
def get_SS_coord_with_ori(ssName):
	coord = re.search(r'(.*):(.*)-(.*)-(.*)', ssName)
	chrom = coord.group(1)
	start = int(coord.group(2))
	end = int(coord.group(3))
	os = coord.group(4)
	ori = '+'
	if (os == '2'):
		ori = '-'
	elif (os == '0'):
		ori = '+.'
	ssInfo = [chrom, start, end, ori]
	return ssInfo

def write_SS_bedfile(allSS_names, outfilename):
	fout = open(outfilename, "w")
	for ssID in allSS_names:
		ssInfo = get_SS_coord_with_ori(ssID)
		fout.write("{}\t{}\t{}\t{}\t1\t{}\n".format(ssInfo[0],ssInfo[1],ssInfo[2],ssID,ssInfo[3]))
	fout.close()

def process_intersect_output(cmd, outfile):
	interSum = defaultdict(dict)
	out, err = run(cmd)
	time.sleep(5)
	fIN = open(outfile, "r")
	outArr = fIN.readlines()
	for line in outArr:
		data = line.split("\t")
		#print(line)
		#print(line,data[3],data[9],data[6],data[7],data[8])
		interSum[data[3]][data[9]] = [ data[6],int(data[7]),int(data[8]) ] #interSum[SS_ID][repeatMask/ExonID] = [ Repeat/Exon coordinates ]
	return interSum
		
#Mark all SS that occur in more than one dataset.

if (__name__ == "__main__"):
	date = time.strftime("%m%d")
	parser = argparse.ArgumentParser(prog="Filtering first pass splice sites")
	parser.add_argument('-dir', dest="inputdir", help="Enter the input directory (STAR output directory)")
	parser.add_argument('-pass', dest="passType", help="Enter whether it is first-pass or second-pass")
	parser.add_argument('-outfile', dest="outfile", help="Enter output filename") 
	args = parser.parse_args()	
	inputdir = args.inputdir
	passType = args.passType
	outputfile = args.outfile
	
	alldir = sorted(os.listdir(inputdir))
	allSS = defaultdict(list)
	#print(alldir)
	for i,seqID in enumerate(alldir):
		print(seqID)
		sj_pairedFile = "{}/{}/SJ.out.tab".format(inputdir,seqID)
		ssDict1 = read_SJout_file(sj_pairedFile)
		if (passType == "first-pass"):
			sj_nonpairedFile = "{}/{}/Unpaired_Run/SJ.out.tab".format(inputdir,seqID)
			ssDict2 = read_SJout_file(sj_nonpairedFile)
			for ssID in ssDict2:
				if (ssID in ssDict1):
					ssDict1[ssID] += ssDict2[ssID]
				else:
					ssDict1[ssID] = ssDict2[ssID]
		for ssID in ssDict1:
			if(ssID not in allSS):
				l = [0] * (i+1)
				l[i] = ssDict1[ssID]
				allSS[ssID] = l
			else:
				allSS[ssID].append(ssDict1[ssID])	
		otherIDs = set(allSS.keys()) - set(ssDict1.keys())
		for oID in otherIDs:
			allSS[oID].append(0)
		print(len(ssDict1),len(otherIDs),len(allSS))
	allSS_names = allSS.keys()
	bedfileName = "../tmp/{}-allSS_bedfile.bed".format(date)
	write_SS_bedfile(allSS_names, bedfileName)
	repeatBedfile = "/scratch/gslab/mpandey/lib/processedFiles/repeatmasked.bed"
	tmpOut = "../tmp/{}-intersect_RM_SS.bed".format(date)
	rmSect = "intersectBed -a {} -b {} -wa -wb > {}".format(bedfileName, repeatBedfile, tmpOut)	
	print(rmSect)
	rmDict = process_intersect_output(rmSect,tmpOut) #Function to process intersect output. Details at function
	print(len(rmDict.keys()))
	#RepeatMasked intersect output is further to exclude those SS that jumps over the repeat region
	# and do not intersect at the 5' or 3' splice site.
	#Concept: If repeat masked coordinates are inclusive into the SS boundaries, then the SS hops over the repeats
	#and are not artefacts.
	quarantSS = []
	for ssID in rmDict:
		for repeat in rmDict[ssID]:
			repCoord = rmDict[ssID][repeat]
			#print(ssID)
			ssCoord = get_SS_coord_with_ori(ssID)
			if (ssCoord[0] == repCoord[0]):
				#if not((ssCoord[1] > int(repCoord[1])) & (ssCoord[2] < int(repCoord[2]))):
				if not((ssCoord[1] < int(repCoord[1])) & (ssCoord[2] > int(repCoord[2]))):
					if (ssID not in quarantSS):
						quarantSS.append(ssID)
	print(quarantSS)
	repRegion = len(quarantSS)
	print("IDs in repeat regions: {}".format(repRegion))
	count_relSS = 0
	ssForSTAR_file = "../lib/SJ_files/{}-SJs_for_STAR_annotation.tab".format(date)
	#outfile = "../output/s01_output/{}-firstPass_SS.tab".format(date)
	fout2 = open(outputfile, "w")
	if (passType == "first-pass"):
		fout1 = open(ssForSTAR_file, "w")
		for ssID in allSS:
			if (ssID not in quarantSS):
				if (sum([i >= 2 for i in allSS[ssID]]) >= 2):
					ssInfo = get_SS_coord_with_ori(ssID)
					fout1.write("{}\t{}\t{}\t{}\n".format(ssInfo[0],ssInfo[1],ssInfo[2],ssInfo[3]))
					readCount = "\t".join(map(str,allSS[ssID]))
					fout2.write("{}:{}-{}\t{}\t{}\n".format(ssInfo[0],ssInfo[1],ssInfo[2],ssInfo[3],readCount))	
					count_relSS += 1
				#else:
				#	if (sum([i >= 2 for i in allSS[ssID]]) >= 1):
				#		tmp = [ i for i in allSS[ssID] if i >= 2 ]
				#		newTmp = sum([ i >= 2 for i in allSS[ssID] ])
				#		readCount = "\t".join(map(str,allSS[ssID]))
				#		print(ssID,readCount,tmp,newTmp)
		fout1.close()
		fout2.close()
	elif (passType == "second-pass"):
		outfile1 = "../output/s02_output/{}-constitutive_splicing.bed".format(date)
		#outfile2 = "../output/s02_output/{}-alternative_splicing.txt".format(date)
		f1 = open(outfile1, "w")
		#f2 = open(outfile2, "w")
		for ssID in allSS:
			if (ssID not in quarantSS):
				if (sum([i >= 1 for i in allSS[ssID]]) >= 2):
					ssInfo = get_SS_coord_with_ori(ssID)
					readCount = "\t".join(map(str,allSS[ssID]))
					fout2.write("{}:{}-{}\t{}\t{}\n".format(ssInfo[0],ssInfo[1],ssInfo[2],ssInfo[3],readCount))
			 		count_relSS += 1
					if (sum([i >= 1 for i in allSS[ssID]]) >= 33):
						temp = [i for i in allSS[ssID] if i >= 1]
						avgCount = sum(temp)/len(temp)
						f1.write("{}\t{}\t{}\t{}\t1\t{}\n".format(ssInfo[0],ssInfo[1],ssInfo[2],avgCount,ssInfo[3]))
	print(count_relSS)
