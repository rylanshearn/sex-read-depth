#!/bin/bash
# ---------------------------------------------------------------------------
# profilesum.sh takes two snp profiles as output from sam2pro, and combines them in a single summary file

# Copyright 2018, Rylan Shearn

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License at <http://www.gnu.org/licenses/> for
# more details.

# usage: profilesum.sh /path/to/profile1.pro /path/to/profile2.pro sample1name sample2name

# Revision history:
# 2017.03.14
# ---------------------------------------------------------------------------

###########################################################
## reshaping data from snp profile
###########################################################

s1pro=$1
s2pro=$2
s1=$3
s2=$4

outputdir=$(dirname $s1pro)

#generate summary of sample 1 profile
#first awk script takes average+1 of four snps
#second awk script inserts zeros for missing coordinates and outputs
awk 'BEGIN {FS=OFS="\t"; print "coord\t'''$s1'''snpM+1"} NR>1{sum=0; n=0; for(i=2;i<=NF;i++) {sum+=$i; ++n} print $1,1+(sum/4)}' $s1pro | awk -v OFS="\t" 'BEGIN {print "coord","'''$s1'''snpM+1";prev_pos=0;} NR>1{ if(prev_pos+1!=int($0)) {for(i=prev_pos+1;i<int($1);++i) {printf("%d\t0\n",i);}} print;prev_pos=int($0);}' > "$outputdir"/"$s1"sum.pro

#generate summary of sample 2 profile
#first awk script takes average+1 of four snps
#second awk script inserts zeros for missing coordinates and outputs
awk 'BEGIN {FS=OFS="\t"; print "coord\t'''$s2'''snpM+1"} NR>1{sum=0; n=0; for(i=2;i<=NF;i++) {sum+=$i; ++n} print $1,1+(sum/4)}' $s2pro | awk -v OFS="\t" 'BEGIN {print "coord","'''$s2'''snpM+1";prev_pos=0;} NR>1{ if(prev_pos+1!=int($0)) {for(i=prev_pos+1;i<int($1);++i) {printf("%d\t0\n",i);}} print;prev_pos=int($0);}' > "$outputdir"/"$s2"sum.pro

#paste both samples together, then run a sliding window over them
paste "$outputdir"/"$s1"sum.pro <(cut -f2 "$outputdir"/"$s2"sum.pro) | awk -v OFS="\t" 'BEGIN{window=150000;slide=10000} { if(NR==1) print "coord",$2,$3} {mod=NR%window; if(NR<=window){count++}else{sum-=array[mod];sum2-=array2[mod]}sum+=$2;sum2+=$3;array[mod]=$2;array2[mod]=$3;} (NR%slide)==0{print NR,sum/count,sum2/count}' > "$outputdir"/"$s1""$s2"sum.pro
