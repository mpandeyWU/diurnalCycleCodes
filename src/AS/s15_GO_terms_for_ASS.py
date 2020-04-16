#!/usr/bin/python

import re, sys, os
from collections import defaultdict
import argparse
import subprocess
import time

def run(cmd):
	p = subprocess.Popen(cmd, shell='True', stdout = subprocess.PIPE)
	p.wait()
	out, errs = p.communicate()
	print(errs)
	return(out)
	
def parse_binary_file(filename):
	uniqSS = defaultdict(list)
	uniqSS = {'G1' : list(), 'SM' : list(), 'Post_SM' :list() }
	with open(filename, "r") as fIN:
		for line in fIN:
			line = line.rstrip("\n")
			data = line.split(",", 1)
			if (data[1] == "1,0,0"):
				uniqSS["G1"].append(data[0])
			elif (data[1] == "0,1,0"):
				uniqSS["SM"].append(data[0])
			elif (data[1] == "0,0,1"):
				uniqSS["Post_SM"].append(data[0])
	return uniqSS

def process_GO_mapping(filename):
	goEnrich = defaultdict(list)
	goDesc = defaultdict(list)
	with open(filename, "r") as fIN:
		for line in fIN:
			line = line.rstrip("\n")
			data = line.split("\t")
			if(data[1].startswith("PF")):
				if (data[1] not in goDesc):
					goDesc[data[1]] = [ data[4],data[3] ]
				if (data[1] in goEnrich):
					if (data[0] not in goEnrich[data[1]]):
						goEnrich[data[1]].append(data[0])
				else:
					goEnrich[data[1]] = [ data[0] ]
				
	return (goEnrich, goDesc)

if (__name__ == "__main__"):
	parser = argparse.ArgumentParser(prog="GO terms for unique hits in the group")
	parser.add_argument('-bfile', dest="binary", help="File that have binary information on LSV presence/absence")
	parser.add_argument('-outdir', dest="outdir", help="File that have binary information on LSV presence/absence")
	args = parser.parse_args()

	binaryfile = args.binary
	outdir = args.outdir

	uniqSS = parse_binary_file(binaryfile)
	uniqFilelist = list()
	for entry in uniqSS:
		outfile = "{}/{}_file.txt".format(outdir,entry)
		fout = open(outfile, "w")
		geneList = list()
		for lsvid in uniqSS[entry]:
			geneName = lsvid.split(".v5.5")
			if (geneName[0] not in geneList):
				geneList.append(geneName[0])
		geneStr = "\n".join(geneList)
		fout.write("{}".format(geneStr))
		fout.close()
		uniqFilelist.append(outfile)
	
	GO_file = "/scratch/gslab/mpandey/lib/GO_terms_Chlamydomonas.tsv"
	
	outfileList = list()
	for filename in uniqFilelist:
		outfile = "{}.gomap".format(filename)
		cmd_awk = "awk -F \"\\t\" \'NR==FNR {{a[$1]; next}} $1 in a {{print $0}}\' {} {} > {}".format(filename,GO_file,outfile)
		print(cmd_awk)
		out = run(cmd_awk)
		outfileList.append(outfile)

	for outfile in outfileList:
		goTerms, goDesc = process_GO_mapping(outfile)
		print(outfile)
		for goid in sorted(goTerms, key=lambda k: len(goTerms[k]), reverse=True): 
			print(goid,len(goTerms[goid]),goDesc[goid])

		

	
			
