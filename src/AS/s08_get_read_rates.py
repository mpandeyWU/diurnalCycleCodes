#!/usr/bin/python

import re, sys, os
sys.path.append('/home/manishipandey/.virtualenvs/majiq-test/lib/python3.5/site-packages/')
from voila.api import SpliceGraph

fout = open("../output/1015_Majiq_All_RC.txt.tmp", "w")
with SpliceGraph('../runs/1006-MajiqBuild_All/splicegraph.sql') as sq:
	for gene in sq.genes:
		for junction in [ js for js in gene.junctions if not js.intron_retention ]:
			for r in junction.reads:
				tmpID = re.search(r'Cre(.*)\.g(.*).v5.5', r.junction_gene_id)
				chrNum = (int(tmpID.group(1)))
				chrID = "chromosome_{}".format(str(chrNum))
				print(gene)
				fout.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(chrID,r.junction_start,r.junction_end,r.junction_gene_id,r.experiment_name,r.reads))
