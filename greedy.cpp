/*
 *
 *  greedy.cpp
 *
 *  ---------------------------------------------------------------------------
 *
 *  Copyright (C) 2019-2021 Cosimo Lupo
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

#include "funcs.h"
#include "greedy.h"

struct seq_align greedy_likelihood(const string read){
	int genomic_V_size, i, max_i, new_i, j, max_j, new_j, k, max_k, n, Lx, Ly, left_bound;
	double S, max_S, logLik;
	double mu_s_err, log_mu_s_err, log_cost_s_match, log_cost_s_error;
	string germline;
	
	struct seq_align greedy_seq;
	
	genomic_V_size = (int)genomic_V.size();
	greedy_seq.seq = read;
	Ly = read.size() + 1;
	
	mu_s_err = params[0];    // initialized by the average value
	log_mu_s_err = log(mu_s_err);
	log_cost_s_match = log((1. - mu_s_err) * (1. - mu_s_err*params[1] - mu_s_err*params[2]));
	log_cost_s_error = log(mu_s_err / 3  * (1. - mu_s_err*params[1] - mu_s_err*params[2]));
	
	greedy_seq.greedy_logLik = 0.;
	
	for(n=0;n<genomic_V_size;n++){
		germline = genomic_V[n].second;
		Lx = germline.size() + 1;
		
		// Allocate the memory
		double** logP = new double*[Lx];
		int** traceback_row = new int*[Lx];
		int** traceback_col = new int*[Lx];
		for(i=0;i<Lx;i++){
			logP[i] = new double[Ly];
			traceback_row[i] = new int[Ly];
			traceback_col[i] = new int[Ly];
		}
		
		// Boundary conditions
		logP[0][0] = 0;
		traceback_row[0][0] = 0;
		traceback_col[0][0] = 0;
		for(i=1;i<Lx;i++){
			left_bound = max(0,i-gap_bound);
			for(k=left_bound;k<i;k++){
				if(k==left_bound){
					S = logP[k][0] + log_mu_s_err + log_gamma_del_array[i-k-1];
					max_S = S;
					max_k = k;
				}else{
					S = logP[k][0] + log_mu_s_err + log_gamma_del_array[i-k-1];
					if(S > max_S){
						max_S = S;
						max_k = k;
					}
				}
			}
			logP[i][0] = max_S;
			traceback_row[i][0] = max_k;
			traceback_col[i][0] = 0;
		}
		for(j=1;j<Ly;j++){
			S = 0.;
			left_bound = max(0,j-gap_bound);
			for(k=left_bound;k<j;k++){
				if(k==left_bound){
					S = logP[0][k] + log_mu_s_err + log_gamma_ins_array[j-k-1];
					max_S = S;
					max_k = k;
				}else{
					S = logP[0][k] + log_mu_s_err + log_gamma_ins_array[j-k-1];
					if(S > max_S){
						max_S = S;
						max_k = k;
					}
				}
			}
			logP[0][j] = max_S;
			traceback_row[0][j] = 0;
			traceback_col[0][j] = max_k;
		}
		
		// Fill the matrix
		for(i=1;i<Lx;i++){
			for(j=1;j<Ly;j++){
				// match/mismatch
				logP[i][j] = logP[i-1][j-1] + (germline.at(i-1)==read.at(j-1) ? log_cost_s_match : log_cost_s_error);
				traceback_row[i][j] = i-1;
				traceback_col[i][j] = j-1;
				// deletion
				left_bound = max(0,i-gap_bound);
				for(k=left_bound;k<i;k++){
					S = logP[k][j] + log_mu_s_err + log_gamma_del_array[i-k-1];
					if(S > logP[i][j]){
						logP[i][j] = S;
						traceback_row[i][j] = k;
						traceback_col[i][j] = j;
					}
				}
				// insertion
				left_bound = max(0,j-gap_bound);
				for(k=left_bound;k<j;k++){
					S = logP[i][k] + log_mu_s_err + log_gamma_ins_array[j-k-1];
					if(S > logP[i][j]){
						logP[i][j] = S;
						traceback_row[i][j] = i;
						traceback_col[i][j] = k;
					}
				}
			}
		}
		
		logLik = logP[Lx-1][Ly-1];
		
		if(n==0){
			greedy_seq.greedy_logLik = logLik;
			greedy_seq.V_best_idx = n;
			// reconstruct the alignment
			greedy_seq.N_err = 0;
			greedy_seq.N_del = 0;
			greedy_seq.N_ins = 0;
			greedy_seq.vec_del.clear();
			greedy_seq.vec_ins.clear();
			greedy_seq.greedy_map = "";
			max_i = Lx-1;
			max_j = Ly-1;
			while(max_i>0 and max_j>0){
				new_i = traceback_row[max_i][max_j];
				new_j = traceback_col[max_i][max_j];
				if(new_i==(max_i-1) and new_j==(max_j-1)){
					// match/mismatch
					if(germline.at(max_i-1) == read.at(max_j-1)){
						greedy_seq.greedy_map = greedy_seq.greedy_map + "m";
					}else if(germline.at(max_i-1) != read.at(max_j-1)){
						greedy_seq.N_err += 1;
						greedy_seq.greedy_map = greedy_seq.greedy_map + "e";
					}
					max_i = new_i;
					max_j = new_j;
				}else if(new_i<max_i and new_j==max_j){
					// deletion
					greedy_seq.N_del += 1;
					greedy_seq.vec_del.push_back(max_i-new_i);
					for(k=0;k<max_i-new_i;k++){
						greedy_seq.greedy_map = greedy_seq.greedy_map + "d";
					}
					max_i = new_i;
					max_j = new_j;
				}else if(new_i==max_i and new_j<max_j){
					// insertion
					greedy_seq.N_ins += 1;
					greedy_seq.vec_ins.push_back(max_j-new_j);
					for(k=0;k<max_j-new_j;k++){
						greedy_seq.greedy_map = greedy_seq.greedy_map + "i";
					}
					max_i = new_i;
					max_j = new_j;
				}
			}
		}else{
			if(logLik > greedy_seq.greedy_logLik){
				greedy_seq.greedy_logLik = logLik;
				greedy_seq.V_best_idx = n;
				// reconstruct the alignment
				greedy_seq.N_err = 0;
				greedy_seq.N_del = 0;
				greedy_seq.N_ins = 0;
				greedy_seq.vec_del.clear();
				greedy_seq.vec_ins.clear();
				max_i = Lx-1;
				max_j = Ly-1;
				while(max_i>0 and max_j>0){
					new_i = traceback_row[max_i][max_j];
					new_j = traceback_col[max_i][max_j];
					if(new_i==(max_i-1) and new_j==(max_j-1)){
						// match/mismatch
						if(germline.at(max_i-1) == read.at(max_j-1)){
							greedy_seq.greedy_map = greedy_seq.greedy_map + "m";
						}else if(germline.at(max_i-1) != read.at(max_j-1)){
							greedy_seq.N_err += 1;
							greedy_seq.greedy_map = greedy_seq.greedy_map + "e";
						}
						max_i = new_i;
						max_j = new_j;
					}else if(new_i<max_i and new_j==max_j){
						// deletion
						greedy_seq.N_del += 1;
						greedy_seq.vec_del.push_back(max_i-new_i);
						for(k=0;k<max_i-new_i;k++){
							greedy_seq.greedy_map = greedy_seq.greedy_map + "d";
						}
						max_i = new_i;
						max_j = new_j;
					}else if(new_i==max_i and new_j<max_j){
						// insertion
						greedy_seq.N_ins += 1;
						greedy_seq.vec_ins.push_back(max_j-new_j);
						for(k=0;k<max_j-new_j;k++){
							greedy_seq.greedy_map = greedy_seq.greedy_map + "i";
						}
						max_i = new_i;
						max_j = new_j;
					}
				}
			}
		}
		
		// Free the memory
		for(i=0;i<Lx;i++){
			delete [] logP[i];
			delete [] traceback_row[i];
			delete [] traceback_col[i];
		}
		delete [] logP;
		delete [] traceback_row;
		delete [] traceback_col;
	}
	
	// Reverse the map for the greedy alignment
	reverse(greedy_seq.greedy_map.begin(), greedy_seq.greedy_map.end());
	
	//return maxLogLik;
	return greedy_seq;
}

struct seq_align reduced_greedy_likelihood(struct seq_align *pt){
	int i, max_i, new_i, j, max_j, new_j, k, max_k, Lx, Ly, left_bound;
	double S, max_S, logLik;
	double mu_s_err, log_mu_s_err, log_cost_s_match, log_cost_s_error;
	string germline, read;
	
	struct seq_align greedy_seq;
	
	greedy_seq = *pt;
	read = pt -> seq;
	Ly = read.size() + 1;
	
	mu_s_err = params[0];    // initialized by the average value
	log_mu_s_err = log(mu_s_err);
	log_cost_s_match = log((1. - mu_s_err) * (1. - mu_s_err*params[1] - mu_s_err*params[2]));
	log_cost_s_error = log(mu_s_err / 3  * (1. - mu_s_err*params[1] - mu_s_err*params[2]));
	
	greedy_seq.greedy_logLik = 0.;
	
	germline = genomic_V[pt -> V_best_idx].second.substr(pt -> V_start,(pt -> V_end)-(pt -> V_start));
	Lx = germline.size() + 1;
	
	// Allocate the memory
	double** logP = new double*[Lx];
	int** traceback_row = new int*[Lx];
	int** traceback_col = new int*[Lx];
	for(i=0;i<Lx;i++){
		logP[i] = new double[Ly];
		traceback_row[i] = new int[Ly];
		traceback_col[i] = new int[Ly];
	}
	
	// Boundary conditions
	logP[0][0] = 0;
	traceback_row[0][0] = 0;
	traceback_col[0][0] = 0;
	for(i=1;i<Lx;i++){
		left_bound = max(0,i-gap_bound);
		for(k=left_bound;k<i;k++){
			if(k==left_bound){
				S = logP[k][0] + log_mu_s_err + log_gamma_del_array[i-k-1];
				max_S = S;
				max_k = k;
			}else{
				S = logP[k][0] + log_mu_s_err + log_gamma_del_array[i-k-1];
				if(S > max_S){
					max_S = S;
					max_k = k;
				}
			}
		}
		logP[i][0] = max_S;
		traceback_row[i][0] = max_k;
		traceback_col[i][0] = 0;
	}
	for(j=1;j<Ly;j++){
		S = 0.;
		left_bound = max(0,j-gap_bound);
		for(k=left_bound;k<j;k++){
			if(k==left_bound){
				S = logP[0][k] + log_mu_s_err + log_gamma_ins_array[j-k-1];
				max_S = S;
				max_k = k;
			}else{
				S = logP[0][k] + log_mu_s_err + log_gamma_ins_array[j-k-1];
				if(S > max_S){
					max_S = S;
					max_k = k;
				}
			}
		}
		logP[0][j] = max_S;
		traceback_row[0][j] = 0;
		traceback_col[0][j] = max_k;
	}
	
	// Fill the matrix
	for(i=1;i<Lx;i++){
		for(j=1;j<Ly;j++){
			// match/mismatch
			logP[i][j] = logP[i-1][j-1] + (germline.at(i-1)==read.at(j-1) ? log_cost_s_match : log_cost_s_error);
			traceback_row[i][j] = i-1;
			traceback_col[i][j] = j-1;
			// deletion
			left_bound = max(0,i-gap_bound);
			for(k=left_bound;k<i;k++){
				S = logP[k][j] + log_mu_s_err + log_gamma_del_array[i-k-1];
				if(S > logP[i][j]){
					logP[i][j] = S;
					traceback_row[i][j] = k;
					traceback_col[i][j] = j;
				}
			}
			// insertion
			left_bound = max(0,j-gap_bound);
			for(k=left_bound;k<j;k++){
				S = logP[i][k] + log_mu_s_err + log_gamma_ins_array[j-k-1];
				if(S > logP[i][j]){
					logP[i][j] = S;
					traceback_row[i][j] = i;
					traceback_col[i][j] = k;
				}
			}
		}
	}
	
	logLik = logP[Lx-1][Ly-1];
	
	greedy_seq.greedy_logLik = logLik;
	// reconstruct the alignment
	greedy_seq.N_err = 0;
	greedy_seq.N_del = 0;
	greedy_seq.N_ins = 0;
	greedy_seq.vec_del.clear();
	greedy_seq.vec_ins.clear();
	greedy_seq.greedy_map = "";
	max_i = Lx-1;
	max_j = Ly-1;
	while(max_i>0 and max_j>0){
		new_i = traceback_row[max_i][max_j];
		new_j = traceback_col[max_i][max_j];
		if(new_i==(max_i-1) and new_j==(max_j-1)){
			// match/mismatch
			if(germline.at(max_i-1) == read.at(max_j-1)){
				greedy_seq.greedy_map = greedy_seq.greedy_map + "m";
			}else if(germline.at(max_i-1) != read.at(max_j-1)){
				greedy_seq.N_err += 1;
				greedy_seq.greedy_map = greedy_seq.greedy_map + "e";
			}
			max_i = new_i;
			max_j = new_j;
		}else if(new_i<max_i and new_j==max_j){
			// deletion
			greedy_seq.N_del += 1;
			greedy_seq.vec_del.push_back(max_i-new_i);
			for(k=0;k<max_i-new_i;k++){
				greedy_seq.greedy_map = greedy_seq.greedy_map + "d";
			}
			max_i = new_i;
			max_j = new_j;
		}else if(new_i==max_i and new_j<max_j){
			// insertion
			greedy_seq.N_ins += 1;
			greedy_seq.vec_ins.push_back(max_j-new_j);
			for(k=0;k<max_j-new_j;k++){
				greedy_seq.greedy_map = greedy_seq.greedy_map + "i";
			}
			max_i = new_i;
			max_j = new_j;
		}
	}
	
	// Free the memory
	for(i=0;i<Lx;i++){
		delete [] logP[i];
		delete [] traceback_row[i];
		delete [] traceback_col[i];
	}
	delete [] logP;
	delete [] traceback_row;
	delete [] traceback_col;
	
	// Reverse the map for the greedy alignment
	reverse(greedy_seq.greedy_map.begin(), greedy_seq.greedy_map.end());
	
	//return maxLogLik;
	return greedy_seq;
}
