/*
 *
 *  funcs.h
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

#ifndef FUNCS_H
#define FUNCS_H

#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <time.h>
#include <math.h>
#include <cmath>
#include <algorithm>
#include <random>
#include <omp.h>
#include <getopt.h>
#include <chrono>
#include <unordered_map>
#include "spline.h"

#define frand (double(rand())/RAND_MAX)
#define sign(x) ((x > 0) ? 1 : ((x < 0) ? -1 : 0))
#define AG 0
#define GAMMA 0
#define SPLINE 1

using namespace std;
using namespace tk;

const int max_gap_bound = 60;

extern int gap_bound;
extern int del_params_begin, ins_params_begin;
extern int N_seqs;
extern vector <pair<const int, const string>> genomic_V;
extern vector <pair<const int, const string>> indexed_seqList;
extern vector <double> params;
extern double gamma_del_array[max_gap_bound], gamma_ins_array[max_gap_bound];
extern double log_gamma_del_array[max_gap_bound], log_gamma_ins_array[max_gap_bound];
extern vector <double> prior;
extern int N_prior;
extern double d_prior, prior_min, prior_max;
extern double *new_prior;
extern double Gauss_KDE_sigma;
extern double frac_mutated;
extern string batch_name;
extern minstd_rand0 generator;
extern double frac_post_update;
extern int N_sampled_mu_errs;
extern int W;

struct synth_seq{
	int seq_ID;  // contains the integer ID of the sequence
	string seq;  // contains the sequence (in nt)
	int V_best_idx;  // integer index for the best V annotation
	string V_best;  // name for the best V annotation
	int V_start;  // starting position of the alignment along the V germline
	int V_end;  // ending position (not included) of the alignment along the V germline
	double mu_err;  // contains the random, sequence-specific point mutation rate
	int N_err;  // contains the number of independent point mutations to be realized
	int N_del;  // contains the number of independent deletions to be realized
	int N_ins;  // contains the number of independent insertions to be realized
	vector <int> vec_del;
	vector <int> vec_ins;
	string error_pos;
	string error_list;
	string insertion_pos;
	string insertion_list;
	string deletion_pos;
	string deletion_list;
};

struct active_cell{
	int row;  // row index - usually labelled as i
	int col;  // column index - usually labelled as j
	vector <int> parents;  // list of parent cells (the one with i-1,j-1, or the ones on the left, or the ones above)
};

struct seq_align{
	int seq_ID;  // contains the integer ID of the sequence
	string seq;  // contains the sequence (in nt)
	string V_best;  // name for the best V annotation
	int V_start;  // starting position of the alignment along the V germline
	int V_end;  // ending position (not included) of the alignment along the V germline
	int V_best_idx;  // integer index for the best V annotation
	int N_err;  // contains the total number of mismatches as given by the greedy alignment
	int N_del;  // contains the number of independent deletions as given by the greedy alignment
	int N_ins;  // contains the number of independent insertions as given by the greedy alignment
	vector <int> vec_del;  // contains the length of each deletion event
	vector <int> vec_ins;  // contains the length of each insertion event
	double greedy_logLik;  // contains the likelihood of the greedy alignment
	string greedy_map;  // encodes the greedy alignment through a sequence of 'm', 'e', 'd' and 'i'
	int N_active_cells;  // number of active cells in the pruned alignment matrix
	struct active_cell *active_cells;  // list of the active cells in the pruned alignment matrix
	int N_parents;  // number of active cells out of the diagonal in the alignment matrix
	bool mutated;  // tells if the sequence has undergone affinity maturation (true, mu_err>0) or not (false, mu_err=0)
	double p_mutated;  // probability in [0,1] that the sequence has undergone affinity maturation
	double mu_err;  // sequence-specific mutation rate
	double Z_prime;  // normalization factor from the update of the continuous part of the posterior
	#if GAMMA
		double alpha;
		double beta;
	#endif
	#if SPLINE
		double mode;
		double mean;
		double sigma;
		vector <double> spline_bins;
		vector <double> spline_posterior;
		vector <double> sampled_mu_errs;
	#endif
};

std::ostream& bold_on(std::ostream& os);
std::ostream& bold_off(std::ostream& os);
bool sort_using_smaller_than(double u, double v);
bool sort_using_greater_than(double u, double v);
vector <string> split_string(string line_str, char separator);
vector <pair<const int, const string>> read_fasta(string filename);
vector <pair<const int , const string>> read_indexed_csv(string filename);
vector <double> educated_init_params(vector <double> p);
void update_gammas();
vector <double> project_by_hand(vector <double> p);
vector <double> simplex_projection(vector <double> p);
void check_model(vector <double> p);
void write_model_file(const string filename, double frac_mutated, vector <double> p);
struct seq_align greedy_likelihood(const string read);
struct seq_align reduced_greedy_likelihood(struct seq_align *pt);
double gammapdf(double x, double alpha, double beta);
double gausspdf(double x, double mean, double sigma);
double gausscdf(double x, double mean, double sigma);
double gausscdf_complement(double x, double mean, double sigma);

#endif
