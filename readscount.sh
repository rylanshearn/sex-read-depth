#!/bin/bash
# ---------------------------------------------------------------------------
# readscount.sh counts the number of reads for two paired end fastq files 

# Copyright 2018, Rylan Shearn

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License at <http://www.gnu.org/licenses/> for
# more details.

# usage: readscount.sh input1.fastq.gz input2.fastq.gz samplename samplenamereadsnum.txt

# Revision history:
# 2016.03.15
# ---------------------------------------------------------------------------

###########################################################
## Count number of reads for sample
###########################################################

input1=$1
input2=$2
samplename=$3

outputdir=$(dirname "${input1}") #get output directory from input

#number of lines in first sample reads files
num1=$(wc -l < "$input1")
num2=$(wc -l < "$input2")

#number of reads for first sample
rnum=$((($num1 + $num2)/4))

#output
echo "$rnum" > "$outputdir"/"$samplename"readsnum.txt
