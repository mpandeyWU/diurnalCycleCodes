#!/usr/bin/python

import re, sys, os
from collections import defaultdict

filename = "../output/altSS_s01_output/May18/0518_Group_based_binaryLSV.csv"
bDict = defaultdict(list)
with open(filename, "r") as fIN:
	for line in fIN:
		if not(line.startswith("#")):
			line = line.rstrip("\n")
			data = line.split(",")
			binaryTag = "".join(data[1:4])
			if (binaryTag in bDict):
				bDict[binaryTag].append(data[0])
			else:
				bDict[binaryTag] = [ data[0] ]
for tag in bDict:
	if (tag == "000"):
		print(bDict[tag])
	print(tag,len(bDict[tag]))
