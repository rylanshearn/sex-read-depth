#!/bin/bash
# ---------------------------------------------------------------------------
# alignbowtiecal.sh aligns male reads to bowtie2 indexed reference genome

# Copyright 2018, Rylan Shearn

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License at <http://www.gnu.org/licenses/> for
# more details.

# usage: samtoolssort.sh /path/to/Referece.fasta /path/to/alignedreads.sam memory(gb)

# Revision history:
# 2016.02.26
# 2016.02.19 added output directory spec
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
