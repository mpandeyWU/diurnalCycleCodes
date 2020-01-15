#!/usr/bin/python

import re, sys, os
import subprocess
import logging
import argparse

def run(cmd):
	p = subprocess.Popen(cmd, shell='True', stdout=subprocess.PIPE)
	p.wait()
	out, errs = p.communicate()
	return (out, errs)


if (__name__ == "__main__"):
	parser = argparse.ArgumentParser(prog="Get Intron Coordinates from GTF file")
	parser.add_argument('-gtf_file', dest = "gtf_file", help="Enter GTF file path of Organism (in this case, C. reinhardtii)")
	parser.add_argument('-outdir', dest = "outdir", help="Enter output directory path")
	args = parser.parse_args()

	gtf_file = args.gtf_file
	outdir = args.outdir
	
	if not(os.path.exists(outdir)):
		os.mkdir(outdir)

	logging.basicConfig(filename = "{}/app.log".format(outdir), filemode = "w", format='%(name)s - %(levelname)s - %(message)s')
	logging.info("Step 1: Get exon coordinates from the GTF file")
	s01_cmd = "awk \'($3==\"exon\")\' {} | gtf2bed - | cut -f1-6 - > {}/a01_Chlamy_exons.bed".format(gtf_file,outdir)
	out1, err1 = run(s01_cmd)
	logging.warning("{}".format(err1))

	logging.info("Step 2: Convert transcripts to merged exons")
	s02_cmd = "awk -f r01_transcripts2MergedExons.awk {}/a01_Chlamy_exons.bed > {}/a02_Chlamy_mergedExons.bed".format(outdir,outdir)
	out2, err2 = run(s02_cmd)
	logging.warning("{}".format(err2))
	
	logging.info("Step 3: Make an Exon - Intron list")
	s03_cmd = "awk -f r02_mergedExons_to_IntronList.awk {}/a02_Chlamy_mergedExons.bed > {}/a03_Chlamy_ExonsandIntrons.bed".format(outdir,outdir)
	out3, err3 = run(s03_cmd)
	logging.warning("{}".format(err3))

	logging.info("Step 4: Convert them into junctions")
	s04_cmd = "awk -f r03_exonIntronList_to_junctionList.awk {}/a03_Chlamy_ExonsandIntrons.bed > {}/a04_Chlamy_exonsIntronJunctions.bed".format(outdir,outdir)
	out4, err4 = run(s04_cmd)
	logging.warning("{}".format(err4))
	
	logging.info("Step 5: Pad the exon-intron junction by 10bp around the junction")
	s05_cmd = "bedops --everything --range 10 {}/a04_Chlamy_exonsIntronJunctions.bed > {}/a05_Chlamy_exonsIntronJunctions.pad10.bed".format(outdir,outdir)
	out5, err5 = run(s05_cmd)
	logging.warning("{}".format(err5))


