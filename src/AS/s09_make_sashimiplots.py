#!/usr/bin/python

import re, sys, os
from collections import defaultdict
import pandas as pd
import argparse
import numpy as np
import time

def parse_gff3_file():
	gff = defaultdict(str)
	gff3_file = "../../lib/PhytozomeV12/Creinhardtii/annotation/Creinhardtii_281_v5.5.gene.gff3"
	with open(gff3_file, "r") as fIN:
		for line in fIN:
			if not(line.startswith("#")):
				line = line.rstrip("\n")
				data = line.split("\t")
				if (data[2] == "gene"):
					coord = "{}:{}-{}".format(data[0],data[3],data[4])
					gene = re.search(r'ID=(.*).v5.5(.*)',data[8])
					geneID = gene.group(1)
					gff[geneID] = coord
	return gff

def parse_timepoints():
	tpDict = defaultdict(dict)
	timepoints_file = "../data/timepoints.txt"
	common_path = "/scratch/gslab/mpandey/Chlamy_SS_model/runs/sorted_bamFiles/"	
	#fout = open(outFilename, "w")
	with open(timepoints_file, "r") as fIN:
		for line in fIN:
			line = line.rstrip("\n")
			data = line.split("\t")
			f2_path = "{}/{}".format(common_path,data[1])
			for i,val in enumerate(list(["A", "B"])):
				smID = "TP{}{}".format(data[2],val)
				fpath = "{}/{}.bam".format(common_path,data[i])	
				tpDict["TP_{}".format(data[2])][smID] = fpath
	return tpDict

def parse_PSI_values(filename):
	df1 = pd.read_csv(filename, sep=",")
	df1 = df1.set_index("ClusterID")
	#print(df1)
	tp_df = df1.drop(["clusterNum"], axis = 1)
	#print(tp_df) 
	clusterMax = tp_df.idxmax(axis = 1)
	tp_df["min_TP"] = tp_df[tp_df > 0].idxmin(axis = 1)
	tp_df["max_TP"] = clusterMax
	sub_df = tp_df[["min_TP", "max_TP"]]
	dfDict = sub_df.to_dict('index')
	return dfDict

def updated_parser_for_PSI_values(filename):
	tpName = None
	psiTP = defaultdict(str)
	with open(filename, "r") as fIN:
		for line in fIN:
			line = line.rstrip("\n")
			line = line.replace('\r', "")
			if not(line.startswith("#")):
				data = line.split(",")
				a = list(map(float,data[2:]))
				filteredPSI = [ (i,psi) for i,psi in enumerate(a) if (psi > 0.0001)]
				sorted_PSI = sorted(filteredPSI, key = lambda x:x[1], reverse = True)
				last = len(sorted_PSI) - 1
				minPSI = sorted_PSI[last]
				maxPSI = sorted_PSI[0]
				val = (float(maxPSI[1]) - float(minPSI[1]))
				#print("{} PSI diff: {}".format(data[0],val))
				if (tpName != None) & ((float(maxPSI[1]) - float(minPSI[1])) >= 0.2):
					minTP = tpName[minPSI[0]].replace("\"", "")
					maxTP = tpName[maxPSI[0]].replace("\"", "")
					psiTP[data[0]] = [ minTP, maxTP ]
					print(data[0],minTP,minPSI,maxTP,maxPSI)
				#print(minval, maxval)
			else:
				data = line.split(",")
				tpName = data[2:]
	return psiTP
				
			
if (__name__ == "__main__"):
	parser = argparse.ArgumentParser(prog = "Enter clustering input for ggsashimi")
	parser.add_argument('-psifile', dest="psifile", help="Input PSI file")
	parser.add_argument('-outdir', dest="outdir", help="Output directory")
	args = parser.parse_args()

	psifile = args.psifile
	outdir = args.outdir
	date = time.strftime("%m%d")

	if not(os.path.isdir(outdir)):
		os.mkdir(outdir)
		os.mkdir("{}/Output".format(outdir))

	tpDict = parse_timepoints()
	geneList = list()	
	#inputType = "Dark"
	#filename = "../output/Clustering_Output/181025_PSI_{}_Set1.csv".format(inputType)
	#psiTP = parse_PSI_values(filename)
	#print(psiTP)
	psiTP = updated_parser_for_PSI_values(psifile)
	print((psiTP))
	geneTP = defaultdict(dict)
	for lsvid in psiTP:
		genetmp = re.search(r'"(.*).v5.5(.*)', lsvid)
		geneID = genetmp.group(1)
		if (geneID not in geneList):
			geneList.append(geneID)
		geneTP[geneID][lsvid] = psiTP[lsvid]
		#print(lsvid[0], psiTP[tmpid])
	
	gffIDs = parse_gff3_file()

	#cmd_file = "../plots/r02_{}_set1_cmdline.sh".format(inputType)
	cmd_file = "../plots/{}_r03_cmdline.sh".format(date)
	summary_file = "{}/summary_of_all_genes.txt".format(outdir)
	sumOut = open(summary_file, "w")
	cmdOut = open(cmd_file, "w")
	sumGene = list()
	for geneID in geneList:
		coord = gffIDs[geneID]
		tsvName = "{}/{}_input.tsv".format(outdir,geneID)
		fout = open(tsvName, "w")
		for lsvid in geneTP[geneID]:
			minTP = geneTP[geneID][lsvid][0]
			maxTP = geneTP[geneID][lsvid][1]
			print(minTP,maxTP)
			minSRR = list(); maxSRR = list()
			for tID in tpDict[minTP]:
				minSRR.append(tpDict[minTP][tID])
				fout.write("{}\t{}\tMin\tOrange\n".format(tID, tpDict[minTP][tID]))
			for tID in tpDict[maxTP]:
				maxSRR.append(tpDict[maxTP][tID])
				fout.write("{}\t{}\tMax\tBlue\n".format(tID, tpDict[maxTP][tID]))
			if not(geneID in sumGene):
				minStr = ";".join(minSRR)
				maxStr = ";".join(maxSRR)
				sumOut.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(geneID,coord,minTP,maxTP,minStr,maxStr))
				sumGene.append(geneID)
		fout.close()
		cmdOut.write("./Mod_sashimi-plot_v01.py -b {} -c {} -o {}/Output/{} -O 3 -C 4 -P pallet.txt -L 1 -A mean\n".format(tsvName,coord,outdir,geneID))		
	#
	#cmdOut.write("./Mod_sashimi-plot_v01.py -b ./timepoints_file.tsv -c {} -o ./Updated_Light_Set1/{} -O 3 -M 2 -C 4 -P pallet.txt -L 1 -A mean\n".format(coord,geneID))		
