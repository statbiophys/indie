#!/bin/bash

#
#  preprocessing.sh
#
#  ---------------------------------------------------------------------------
#
#  Copyright (C) 2019-2022 Natanael Spisak, Cosimo Lupo
#
#  This source code is distributed as part of the 'indie' software.
#  'indie' (INference on Deletion and InsErtions) is a versatile software
#  for evaluation and inference of indel and point substitution hypermutations
#  in high-throughput Ig antibody sequencing data, as well as
#  for the generation of synthetic repertoires with custom models.
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  You should have received a copy of the GNU General Public License
#  along with this program. If not, see <https://www.gnu.org/licenses/>.
#
#  ---------------------------------------------------------------------------
#
#  For any issue or question, please send an email to <cosimo.lupo89@gmail.com>.
#

########## settings ##########

# data location
data_dir = "./"

# primers
primers_dir = "./"
CPrimers = $primers_dir"CPrimers.fasta"
VPrimers = $primers_dir"VPrimers.fasta"

#UMI counts
Umithr = 3

########## main ##########

FilterSeq.py quality -s $data_dir"R1.fastq" -q 30 --outname $data_dir"all_R1" --log $data_dir"FS1.log"
FilterSeq.py quality -s $data_dir"R2.fastq" -q 30 --outname $data_dir"all_R2" --log $data_dir"FS2.log"

MaskPrimers.py align -s $data_dir"all_R1_quality-pass.fastq" --maxlen 50 -p $CPrimers --mode cut --barcode --outname $data_dir"all_R1" --log $data_dir"MP1.log"
MaskPrimers.py align -s $data_dir"all_R2_quality-pass.fastq" --maxlen 30 -p $VPrimers --mode mask --outname $data_dir"all_R2" --log $data_dir"MP2.log"

PairSeq.py -1 $data_dir"all_R1_primers-pass.fastq" -2 $data_dir"all_R2_primers-pass.fastq" --1f BARCODE --coord sra

ClusterSets.py set -s $data_dir"new_all_R1_primers-pass_pair-pass.fastq" -f BARCODE -k INDEX_SEQ --ident 0.9 --exec usearch

ParseHeaders.py merge -s $data_dir"new_all_R1_primers-pass_pair-pass_cluster-pass.fastq" -f BARCODE INDEX_SEQ -k INDEX_MERGE

PairSeq.py -1 $data_dir"new_all_R1_primers-pass_pair-pass_cluster-pass_reheader.fastq" -2 $data_dir"new_all_R2_primers-pass_pair-pass.fastq" --1f INDEX_MERGE --coord sra

BuildConsensus.py -s $data_dir"new_all_R1_primers-pass_pair-pass_cluster-pass_reheader_pair-pass.fastq" --bf INDEX_MERGE --pf PRIMER --prcons 0.6 --maxerror 0.1 --outname $data_dir"all_R1" --log $data_dir"BC1.log"
BuildConsensus.py -s $data_dir"new_all_R2_primers-pass_pair-pass_pair-pass.fastq" --bf INDEX_MERGE --pf PRIMER --maxerror 0.1 --outname $data_dir"all_R2" --log $data_dir"BC2.log"

PairSeq.py -1 $data_dir"all_R1_consensus-pass.fastq" -2 $data_dir"all_R2_consensus-pass.fastq" --1f INDEX_MERGE --coord presto

AssemblePairs.py align -1 $data_dir"all_R2_consensus-pass_pair-pass.fastq" -2 $data_dir"all_R1_consensus-pass_pair-pass.fastq" --coord presto --rc tail --1f CONSCOUNT --2f CONSCOUNT PRCONS --outname $data_dir"all" --log $data_dir"AP.log"

ParseHeaders.py collapse -s $data_dir"all_assemble-pass.fastq" -f CONSCOUNT --act min

CollapseSeq.py -s $data_dir"all_assemble-pass_reheader.fastq" -n 0 --inner --uf PRCONS --cf CONSCOUNT --act sum --outname $data_dir"all"

SplitSeq.py group -s $data_dir"all_collapse-unique.fastq" -f CONSCOUNT --num $Umithr --outname $data_dir"all"

ParseHeaders.py table -s $data_dir"all_atleast-"$Umithr".fastq" -f ID PRCONS CONSCOUNT DUPCOUNT
ParseHeaders.py table -s $data_dir"all_under-"$Umithr".fastq" -f ID PRCONS CONSCOUNT DUPCOUNT

seqtk seq -a $data_dir"all_atleast-"$Umithr".fastq" > $data_dir"all_atleast-"$Umithr".fasta"
seqtk seq -a $data_dir"all_under-"$Umithr".fastq" > $data_dir"all_under-"$Umithr".fasta"
