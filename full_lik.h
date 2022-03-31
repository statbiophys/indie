/*
 *
 *  full_lik.h
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

#ifndef FULL_LIK_H
#define FULL_LIK_H

extern int max_T;
extern double alpha, eps, M;
extern double prune_thr;
extern double p_mutated;

double full_likelihood(struct seq_align *pt);
void prune_by_new_method(struct seq_align *pt);
double partial_likelihood_new_method(struct seq_align *pt);
double negLogLikelihood(struct seq_align *pt, string method);
#if GAMMA
	void init_gamma_posterior(struct seq_align *pt);
	double update_gamma_posterior(struct seq_align *pt, int t);
#endif
#if SPLINE
	void init_spline_posterior(struct seq_align *pt);
	void update_spline_posterior(struct seq_align *pt, int t);
#endif
vector <double> smooth_prior(vector <double> prior, double prior_mean, string method);
void tune_KDE_sigma(string batch_name, int attempt_number, vector <double> prior_mean_history, vector <double> prior_std_history, int t_w);
void explore_scenarios(struct seq_align *pt);

#endif
