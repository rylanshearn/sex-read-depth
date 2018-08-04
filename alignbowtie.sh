#!/bin/bash
# ---------------------------------------------------------------------------
# alignbowtiecal.sh aligns male reads to bowtie2 indexed reference genome

# Copyright 2018, Rylan Shearn

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License at <http://www.gnu.org/licenses/> for
# more details.

# usage: alignbowtiecal.sh RefereceIndexedname /path/to/forwardreads.fastq.gz /path/to/reversereads.fastq.gz /path/to/output.sam <number of threads>

# relies on ReferenceIndexedname (indexed reference file) to be in the same directory as the reads files
# outputs to same directory as reads files

# Revision history:
# 2016.02.16
# ---------------------------------------------------------------------------

reference=$1
forward=$2
reverse=$3
output=$4
nodes=$5
alength=$6

outputdir=$(dirname "${forward}") #get output directory from reference

bowtie2 -p $nodes -N 0 -L $alength -x "$outputdir"/"$reference" -1 $forward -2 $reverse -S "$outputdir"/"$output".sam

