#!/usr/bin/bash


for name in `cat ./timepoints.txt`
do
	#cmd="convert2bed -i bam -o bed < ../runs/sorted_bamFiles/$name.bam > ../runs/bedfiles/$name.bed"
	#cmd="bedtools intersect -a ../output/IR_s01_output/a05_Chlamy_exonsIntronJunctions.pad10.bed -b ../runs/bedfiles/$name.bed -f 1.0 -wo | awk -F \"\\t\" '!(\$14 ~ /N/) {print \$0;}' > ../runs/IR_juncMap/$name.bed"
	#cmd="bedtools intersect -a ../runs/bedfiles/$name.bed -b /scratch/gslab/mpandey/lib/updated_IntronCoord.bed -f 1.0 -wo > ../runs/IR_withinMap/$name.bed"
	#cmd="bedtools coverage -a /scratch/gslab/mpandey/lib/updated_IntronCoord.bed -b ../runs/sorted_bamFiles/$name.bam -F 0.10 -split > ../runs/IntronCoverage/$name.bed"
	cmd="samtools view -bq 10 ../runs/sorted_bamFiles/$name.bam | bedtools coverage -a /scratch/gslab/mpandey/lib/updated_IntronCoord.bed -b stdin -F 0.10 -split > ../runs/IntronCoverage/$name.bed"
	echo $cmd
done
