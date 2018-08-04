#!/bin/bash
# ---------------------------------------------------------------------------
# bowiteindex.sh counts the number of reads for two paired end fastq files 

# Copyright 2018, Rylan Shearn

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License at <http://www.gnu.org/licenses/> for
# more details.

# usage: bowtieindex.sh reference.fasta referenceindexed

# Revision history:
# 2016.03.15
# ---------------------------------------------------------------------------

###########################################################
## Index fasta reference with bowtie2
###########################################################

reference=$1
indexed=$2

#get output directory
outputdir=$(dirname "${reference}") #get output directory from reference

#index the reference
bowtie2-build "$reference" "$outputdir"/"$indexed"
