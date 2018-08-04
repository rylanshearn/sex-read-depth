#!/bin/bash
# ---------------------------------------------------------------------------
# samtoolsdepth.sh takes two sorted bam files and the number of reads from a male and female sample, and calculates the read depth ratio between them

# Copyright 2018, Rylan Shearn

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License at <http://www.gnu.org/licenses/> for
# more details.

# usage: samtoolsdepth.sh /path/to/sample1_sorted.bam /path/to/sample2_sorted.bam sample1name sample2name sample1readsnum.txt sample2readsnum.txt

# Revision history:
# 2016.03.16
# ---------------------------------------------------------------------------

###########################################################
## read depth and reshaping data
###########################################################

s1bam=$1
s2bam=$2
s1=$3
s2=$4
s1readsnum=$5
s2readsnum=$6

outputdir=$(dirname $s1bam)
r1num=$(cat "$s1readsnum") #get number of reads for sample 1
r2num=$(cat "$s2readsnum") #get number of reads for sample 2

#pipe samtools depth to awk (adds zeros to front, then fills in gaps with zeros and corrects for number of reads)
#1 - prints header
#2 - 0s will continue being printed until coord value >= 1
#3 - default 0s line printed if the prev coord not 'current coord - 1' and from same scaffold
#4 - if coverage=0 for 2nd sample (male) but not 1st, a line is printed that gives ratio=0, this is a workaround for the fact that you cannot divide by 0
#5 - line printed with ratio calculation if only 2nd sample, or both samples coverage not equal to 0
samtools depth "$s1bam" "$s2bam" | awk  -v OFS="\t" 'BEGIN {print "mCoord","chr","coord","samp'''$s1'''","samp'''$s2'''","corr'''$s1'''","corr'''$s2'''","ratio";} NR==1{for(i=1;i<$2;i++) print ++cont,$1,i,0,0,0,0,"NA"} NR>1 && $2!=exp_idx{for (i=exp_idx;i<$2;i++){printf("%d\t%s\t%d\t0\t0\t0\t0\tNA\n",++cont,exp_coord,i)}} {if ($3!=0 && $4==0) print ++cont,$0,$3/'''$r1num''',$4/'''$r2num''',"NA";exp_coord=$1;exp_idx=$2+1;} {if ($3!=0 && $4!=0 || $3==0 && $4!=0) print ++cont,$0,$3/'''$r1num''',$4/'''$r2num''',(($3/'''$r1num''')/($4/'''$r2num'''));exp_coord=$1;exp_idx=$2+1}' > "$outputdir"/"$s1""$s2"coverageCorr.txt

