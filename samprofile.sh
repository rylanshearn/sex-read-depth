#!/bin/bash
# ---------------------------------------------------------------------------
# samprofile.sh generates polymorphism profile file from a .sam alignment

# Copyright 2018, Rylan Shearn

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License at <http://www.gnu.org/licenses/> for
# more details.

# usage: samprofile.sh /path/to/Referece.fasta /path/to/alignedreads.sam memory(gb)

# Revision history:
# 2018.08.04
# ---------------------------------------------------------------------------

reference=$1
input=$2
memory=$3
sample=$(basename $input .sam) #get name of input file
outputdir=$(dirname "${input}") #get output directory from input

#sorting prep in samtools (15 mins)
samtools view -t "$reference" -F 4 -h -S -b -o "$outputdir"/"$sample".bam "$input"

#sorting in samtools (20 mins)
samtools sort -m "$memory"G -o "$outputdir"/"$sample"_sorted.bam "$outputdir"/"$sample".bam

#create .bai index file (instant)
samtools index "$outputdir"/"$sample"_sorted.bam

#generate profile file with mpileup and sam2pro
samtools view -b "$outputdir"/"$sample"_sorted.bam | samtools mpileup - | sam2pro -c 5 > "$outputdir"/"$sample".pro
