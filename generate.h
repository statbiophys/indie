/*
 *
 *  generate.h
 *
 *  ---------------------------------------------------------------------------
 *
 *  Copyright (C) 2019-2022 Cosimo Lupo
 *
 *  This source code is distributed as part of the 'indie' software.
 *  'indie' (INference on Deletion and InsErtions) is a versatile software
 *  for evaluation and inference of indel and point substitution hypermutations
 
 *  in high-throughput Ig antibody sequencing data, as well as
 *  for the generation of synthetic repertoires with custom models.
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  You should have received a copy of the GNU General Public License
 *  along with this program. If not, see <https://www.gnu.org/licenses/>.
 *
 *  ---------------------------------------------------------------------------
 *
 *  For any issue or question, please send an email to <cosimo.lupo89@gmail.com>.
 *
 */

#ifndef GENERATE_H
#define GENERATE_H

const string nuclAlphabet[4] = {"A", "C", "G", "T"};
const string notA[3] = {"C", "G", "T"};
const string notC[3] = {"A", "G", "T"};
const string notG[3] = {"A", "C", "T"};
const string notT[3] = {"A", "C", "G"};

void init_synth_seqs(int N_seqs_synth, struct synth_seq *pt, unordered_map <int, string> genomic_V_map_inverse, unordered_map <string, int> genomic_V_map);
void generate_synth_seqs_type_1(int N_seqs_synth, struct synth_seq *pt);
void generate_synth_seqs_type_2(int N_seqs_synth, struct synth_seq *pt);
void write_generate_files(int N_seqs_synth, struct synth_seq *pt, string batch_name);

#endif
