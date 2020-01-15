#!/usr/bin/python

import re, sys, os
import argparse
from collections import defaultdict
import time

if (__name__ == "__main__"):
	date = time.strftime("%m%d")
	parser = argparse.ArgumentParser(prog="Make MAJIQ and VOILA run files")
	parser.add_argument('-tp', dest="timepoints", help="Enter timepoints to analyze (separated by comma)")
	parser.add_argument('-mdir', dest="mdir", help="Path to the directory where *.majiq files are located")
	parser.add_argument('-outmajiq', dest="outmajiq", help="MAJIQ PSI output directory")
	parser.add_argument('-outvoila', dest="outvoila", help="VOILA PSI output directory")
	parser.add_argument('--delta', dest="cmdType", const=1; default=0, help="VOILA PSI output directory")
	
	args = parser.parse_args()	
	tp = args.timepoints
	majiqdir = args.mdir
	outmajiq = args.outmajiq
	outvoila = args.outvoila
	cmdType = args.cmdType
	
	timepoints = list()
	if (tp != "all"):
		timepoints = tp.split(",")
	else:
		timepoints = list(map(str, range(1,25)))
		print(timepoints)
		timepoints.extend(["11.5", "12.5", "13.5", "14.5"])
	print(timepoints)

	fout1 = open("./r11_majiq_psi_run.sh", "w")
	fout2 = open("./r12_voila_psi_run.sh","w")
	idfile = "../data/timepoints.txt"
	idInfo = defaultdict(list)
	with open(idfile, "r") as fIN:
		for line in fIN:
			data = line.split("\t")
			print(data)
			idInfo[data[2]] = [ data[0],data[1] ]
	
	majiq_cmd = "srun majiq psi -j 4 -n TP_{} -o {} {} {}"
	voila_cmd = "srun voila psi -s {}/splicegraph.sql -o {} {}"
	if (cmdType):
		majiq_cmd = "srun majiq deltapsi "
	for tp in timepoints:
		if (tp in idInfo):
			f1Path = "{}/{}.majiq".format(majiqdir,idInfo[tp][0])
			f2Path = "{}/{}.majiq".format(majiqdir,idInfo[tp][1])
			outdir = "{}/TP_{}".format(outmajiq,tp)
			majiq_psi = majiq_cmd.format(tp,outdir,f1Path,f2Path)
			fout1.write("{}\n".format(majiq_psi))
			voilafile = "{}/TP_{}.psi.voila".format(outdir,tp)
			voila_psi = voila_cmd.format(majiqdir,outdir,voilafile)
			fout2.write("{}\n".format(voila_psi))
		
