/*
 *
 *  full_lik.cpp
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
#include "full_lik.h"

double full_likelihood(struct seq_align *pt){
	
	int genomic_V_size, i, j, k, n, Lx, Ly;
	double S, lik;
	double mu_s_err, cost_s_match, cost_s_error;
	string germline, read;
	
	genomic_V_size = (int)genomic_V.size();
	read = pt -> seq;
	Ly = read.size() + 1;
	
	mu_s_err = pt -> mu_err;
	cost_s_match = (1. - mu_s_err) * (1. - mu_s_err*params[1] - mu_s_err*params[2]);
	cost_s_error = mu_s_err / 3  * (1. - mu_s_err*params[1] - mu_s_err*params[2]);
	
	lik = 0.;
	
	for(n=0;n<genomic_V_size;n++){
		germline = genomic_V[n].second;
		Lx = germline.size() + 1;
		
		// Allocate the memory
		double** P = new double*[Lx];
		for(i=0;i<Lx;i++){
			P[i] = new double[Ly];
		}
		
		// Fill the first column
		P[0][0] = 1;
		for(i=1;i<gap_bound;i++){
			S = 0.;
			// first place from above is 0
			for(k=0;k<i;k++){
				// actual length of the deletion is i-k
				S += P[k][0] * mu_s_err * gamma_del_array[i-k-1];
			}
			P[i][0] = S;
		}
		for(i=gap_bound;i<Lx;i++){
			S = 0.;
			// first place from above is i-gap_bound
			for(k=i-gap_bound;k<i;k++){
				// actual length of the deletion is i-k
				S += P[k][0] * mu_s_err * gamma_del_array[i-k-1];
			}
			P[i][0] = S;
		}
		
		// Fill the first row
		for(j=1;j<gap_bound;j++){
			S = 0.;
			// first place from left is 0
			for(k=0;k<j;k++){
				// actual length of the insertion is j-k
				S += P[0][k] * mu_s_err * gamma_ins_array[j-k-1];
			}
			P[0][j] = S;
		}
		for(j=gap_bound;j<Ly;j++){
			S = 0.;
			// first place from left is j-gap_bound
			for(k=j-gap_bound;k<j;k++){
				// actual length of the insertion is j-k
				S += P[0][k] * mu_s_err * gamma_ins_array[j-k-1];
			}
			P[0][j] = S;
		}
		
		// Fill the matrix
		for(i=1;i<gap_bound;i++){
			for(j=1;j<gap_bound;j++){
				S = 0.;
				S += P[i-1][j-1] * (germline.at(i-1)==read.at(j-1) ? cost_s_match : cost_s_error);
				// first place from above is 0
				for(k=0;k<i;k++){
					// actual length of the deletion is i-k
					S += P[k][j] * mu_s_err * gamma_del_array[i-k-1];
				}
				// first place from left is 0
				for(k=0;k<j;k++){
					// actual length of the insertion is j-k
					S += P[i][k] * mu_s_err * gamma_ins_array[j-k-1];
				}
				P[i][j] = S;
			}
			for(j=gap_bound;j<Ly;j++){
				S = 0.;
				S += P[i-1][j-1] * (germline.at(i-1)==read.at(j-1) ? cost_s_match : cost_s_error);
				// first place from above is 0
				for(k=0;k<i;k++){
					// actual length of the deletion is i-k
					S += P[k][j] * mu_s_err * gamma_del_array[i-k-1];
				}
				// first place from left is j-gap_bound
				for(k=j-gap_bound;k<j;k++){
					// actual length of the insertion is j-k
					S += P[i][k] * mu_s_err * gamma_ins_array[j-k-1];
				}
				P[i][j] = S;
			}
		}
		for(i=gap_bound;i<Lx;i++){
			for(j=1;j<gap_bound;j++){
				S = 0.;
				S += P[i-1][j-1] * (germline.at(i-1)==read.at(j-1) ? cost_s_match : cost_s_error);
				// first place from above is i-gap_bound
				for(k=i-gap_bound;k<i;k++){
					// actual length of the deletion is i-k
					S += P[k][j] * mu_s_err * gamma_del_array[i-k-1];
				}
				// first place from left is 0
				for(k=0;k<j;k++){
					// actual length of the insertion is j-k
					S += P[i][k] * mu_s_err * gamma_ins_array[j-k-1];
				}
				P[i][j] = S;
			}
			for(j=gap_bound;j<Ly;j++){
				S = 0.;
				S += P[i-1][j-1] * (germline.at(i-1)==read.at(j-1) ? cost_s_match : cost_s_error);
				// first place from above is i-gap_bound
				for(k=i-gap_bound;k<i;k++){
					// actual length of the deletion is i-k
					S += P[k][j] * mu_s_err * gamma_del_array[i-k-1];
				}
				// first place from left is j-gap_bound
				for(k=j-gap_bound;k<j;k++){
					// actual length of the insertion is j-k
					S += P[i][k] * mu_s_err * gamma_ins_array[j-k-1];
				}
				P[i][j] = S;
			}
		}
		
		lik += P[Lx-1][Ly-1];
		
		// Free the memory
		
		for(i=0;i<Lx;i++){
			delete [] P[i];
		}
		delete [] P;
	}
	
	return lik;
}

void prune_by_new_method(struct seq_align *pt){
	
	int i, j, k, Lx, Ly, nn, N_cells, N_parents;
	double S, max_P;
	double mu_s_err, cost_s_match, cost_s_error;
	string germline, read;
	pair <int, int> cell;
	vector <pair <int, int>> active_cells_list;
	
	read = pt -> seq;
	Ly = read.size() + 1;
	//germline = genomic_V[pt -> V_best_idx].second;
	germline = genomic_V[pt -> V_best_idx].second.substr(pt -> V_start,(pt -> V_end)-(pt -> V_start));
	Lx = germline.size() + 1;
	
	mu_s_err = pt -> mu_err;
	cost_s_match = (1. - mu_s_err) * (1. - mu_s_err*params[1] - mu_s_err*params[2]);
	cost_s_error = mu_s_err / 3  * (1. - mu_s_err*params[1] - mu_s_err*params[2]);
	
	// Allocate the memory for the fw matrix
	double** fw_P = new double*[Lx];
	for(i=0;i<Lx;i++){
		fw_P[i] = new double[Ly];
	}
	
	fw_P[0][0] = 1;
	// Fill the first column
	for(i=1;i<Lx;i++){
		S = 0.;
		for(k=max(0,i-gap_bound);k<i;k++){
			// actual length of the deletion is i-k
			S += fw_P[k][0] * mu_s_err * gamma_del_array[i-k-1];
		}
		fw_P[i][0] = S;
	}
	// Fill the first row
	for(j=1;j<Ly;j++){
		S = 0.;
		for(k=max(0,j-gap_bound);k<j;k++){
			// actual length of the insertion is j-k
			S += fw_P[0][k] * mu_s_err * gamma_ins_array[j-k-1];
		}
		fw_P[0][j] = S;
	}
	// Fill the rest of the matrix
	for(i=1;i<Lx;i++){
		for(j=1;j<Ly;j++){
			S = 0.;
			S += fw_P[i-1][j-1] * (germline.at(i-1)==read.at(j-1) ? cost_s_match : cost_s_error);
			for(k=max(0,i-gap_bound);k<i;k++){
				// actual length of the deletion is i-k
				S += fw_P[k][j] * mu_s_err * gamma_del_array[i-k-1];
			}
			for(k=max(0,j-gap_bound);k<j;k++){
				// actual length of the insertion is j-k
				S += fw_P[i][k] * mu_s_err * gamma_ins_array[j-k-1];
			}
			fw_P[i][j] = S;
		}
	}
	
	// Reverse both germline and read
	reverse(germline.begin(), germline.end());
	reverse(read.begin(), read.end());
	
	// Allocate the memory for the bw matrix
	double** bw_P = new double*[Lx];
	for(i=0;i<Lx;i++){
		bw_P[i] = new double[Ly];
	}
	
	bw_P[0][0] = 1;
	// Fill the first column
	for(i=1;i<Lx;i++){
		S = 0.;
		for(k=max(0,i-gap_bound);k<i;k++){
			// actual length of the deletion is i-k
			S += bw_P[k][0] * mu_s_err * gamma_del_array[i-k-1];
		}
		bw_P[i][0] = S;
	}
	// Fill the first row
	for(j=1;j<Ly;j++){
		S = 0.;
		for(k=max(0,j-gap_bound);k<j;k++){
			// actual length of the insertion is j-k
			S += bw_P[0][k] * mu_s_err * gamma_ins_array[j-k-1];
		}
		bw_P[0][j] = S;
	}
	// Fill the rest of the matrix
	for(i=1;i<Lx;i++){
		for(j=1;j<Ly;j++){
			S = 0.;
			S += bw_P[i-1][j-1] * (germline.at(i-1)==read.at(j-1) ? cost_s_match : cost_s_error);
			for(k=max(0,i-gap_bound);k<i;k++){
				// actual length of the deletion is i-k
				S += bw_P[k][j] * mu_s_err * gamma_del_array[i-k-1];
			}
			for(k=max(0,j-gap_bound);k<j;k++){
				// actual length of the insertion is j-k
				S += bw_P[i][k] * mu_s_err * gamma_ins_array[j-k-1];
			}
			bw_P[i][j] = S;
		}
	}
	
	// Reverse back both germline and read
	reverse(germline.begin(), germline.end());
	reverse(read.begin(), read.end());
	
	// Allocate the memory for the full matrix
	double** P = new double*[Lx];
	for(i=0;i<Lx;i++){
		P[i] = new double[Ly];
	}
	
	// Compute (and rescale) the full matrix from fw and bw matrices
	max_P = fw_P[0][0] * bw_P[Lx-1][Ly-1];
	for(i=0;i<Lx;i++){
		for(j=0;j<Ly;j++){
			P[i][j] = fw_P[i][j] * bw_P[Lx-i-1][Ly-j-1];
			if(P[i][j] > max_P){
				max_P = P[i][j];
			}
		}
	}
	for(i=0;i<Lx;i++){
		for(j=0;j<Ly;j++){
			P[i][j] = P[i][j] / max_P;
		}
	}
	
	// Identify and store the cells above threshold
	active_cells_list.clear();
	for(i=0;i<Lx;i++){
		for(j=0;j<Ly;j++){
			if( P[i][j] > prune_thr ){
				cell.first = i;
				cell.second = j;
				active_cells_list.push_back(cell);
			}
		}
	}
	
	// Check over first and last active cells
	cell = active_cells_list.at(0);
	if(cell.first!=0 or cell.second!=0){
		cout << "\nERROR! active_cell_list vector does not have the correct initial cell!\n" << endl;
		exit(EXIT_FAILURE);
	}
	cell = active_cells_list.at((int)active_cells_list.size()-1);
	if(cell.first!=(Lx-1) or cell.second!=(Ly-1)){
		cout << "\nERROR! active_cell_list vector does not have the correct final cell!\n" << endl;
		exit(EXIT_FAILURE);
	}
	
	// new way of storing active cells
	// now each cell also contains its parent cells (listed as integers)
	// we are not storing the parent cell on the diagonal (corresponding to matches/mismatches)
	N_cells = (int)active_cells_list.size();
	struct active_cell *new_active_cells_list;
	new_active_cells_list = new active_cell[N_cells];
	N_parents = 0;
	
	new_active_cells_list[0].row = 0;
	new_active_cells_list[0].col = 0;
	
	for(nn=1;nn<N_cells;nn++){
		i = active_cells_list.at(nn).first;
		j = active_cells_list.at(nn).second;
		new_active_cells_list[nn].row = i;
		new_active_cells_list[nn].col = j;
		for(k=0;k<N_cells;k++){
			if( active_cells_list.at(k).first<=i or active_cells_list.at(k).second<=j ){
				/*if( active_cells_list.at(k).first==i-1 and active_cells_list.at(k).second==j-1 ){
				 // match/mismatch
				 new_active_cells_list[nn].parents.push_back(k);
				 }else */
				if( active_cells_list.at(k).first==i and active_cells_list.at(k).second<j and active_cells_list.at(k).second>=j-gap_bound){
					// insertion
					new_active_cells_list[nn].parents.push_back(k);
					N_parents++;
				}else if( active_cells_list.at(k).first<i and active_cells_list.at(k).first>=i-gap_bound and active_cells_list.at(k).second==j){
					// deletion
					new_active_cells_list[nn].parents.push_back(k);
					N_parents++;
				}
			}
		}
	}
	
	//pt -> active_cells = active_cells_list;
	pt -> N_active_cells = N_cells;
	pt -> active_cells = new_active_cells_list;
	pt -> N_parents = N_parents;
	
	// Free all the matrices
	for(i=0;i<Lx;i++){
		delete [] fw_P[i];
		delete [] bw_P[i];
		delete [] P[i];
	}
	delete [] fw_P;
	delete [] bw_P;
	delete [] P;
	
	return;
}

double partial_likelihood_new_method(struct seq_align *pt){
	
	int i, parent_i, j, parent_j, Lx, Ly, nn, N_cells, pp;
	double prunedLik;
	double mu_s_err, cost_s_match, cost_s_error;
	string germline, read;
	//pair <int, int> cell;
	//vector <pair <int, int>> active_cells_list;
	struct active_cell cell;
	struct active_cell parent_cell;
	
	read = pt -> seq;
	Ly = read.size() + 1;
	
	mu_s_err = pt -> mu_err;
	cost_s_match = (1. - mu_s_err) * (1. - mu_s_err*params[1] - mu_s_err*params[2]);
	cost_s_error = mu_s_err / 3  * (1. - mu_s_err*params[1] - mu_s_err*params[2]);
	
	//germline = genomic_V[pt -> V_best_idx].second;
	germline = genomic_V[pt -> V_best_idx].second.substr(pt -> V_start,(pt -> V_end)-(pt -> V_start));
	Lx = germline.size() + 1;
	//active_cells_list = pt -> active_cells;
	N_cells = pt -> N_active_cells;
	
	// Allocate and initialize the memory for the pruned matrix
	double** pruned_P = new double*[Lx];
	for(i=0;i<Lx;i++){
		pruned_P[i] = new double[Ly];
	}
	for(i=0;i<Lx;i++){
		for(j=0;j<Ly;j++){
			pruned_P[i][j] = 0.;
		}
	}
	
	pruned_P[0][0] = 1;
	// First active cell is (0,0)
	// so nn index has to start from 1
	for(nn=1;nn<N_cells;nn++){
		cell = (pt -> active_cells)[nn];
		i = cell.row;
		j = cell.col;
		if(i>0 and j>0){
			pruned_P[i][j] += pruned_P[i-1][j-1] * (germline.at(i-1)==read.at(j-1) ? cost_s_match : cost_s_error);
		}
		for(pp=0;pp<(int)(cell.parents.size());pp++){
			parent_cell = (pt -> active_cells)[cell.parents.at(pp)];
			parent_i = parent_cell.row;
			parent_j = parent_cell.col;
			/*
			 if(i==(parent_i+1) and j==(parent_j+1)){
			 //pruned_P[i][j] += pruned_P[i-1][j-1] * (germline.at(i-1)==read.at(j-1) ? cost_s_match : cost_s_error);
			 }else */
			if(j==parent_j){
				// actual length of the deletion is i-parent_i
				pruned_P[i][j] += pruned_P[parent_i][j] * mu_s_err * gamma_del_array[i-parent_i-1];
			}else if(i==parent_i){
				// actual length of the insertion is j-parent_j
				pruned_P[i][j] += pruned_P[i][parent_j] * mu_s_err * gamma_ins_array[j-parent_j-1];
			}
		}
	}
	
	prunedLik = pruned_P[Lx-1][Ly-1];
	
	// Free all the matrices
	for(i=0;i<Lx;i++){
		delete [] pruned_P[i];
	}
	delete [] pruned_P;
	
	return prunedLik;
	
}

double negLogLikelihood(struct seq_align *pt, string method){
	int n;
	double negLogLik = 0.;
	
	if(method=="full_lik"){
		// computes the full alignment likelihood
		#pragma omp parallel for reduction (+:negLogLik)
		for(n=0;n<N_seqs;n++){
			//negLogLik -= log(full_likelihood(indexed_seqList[n].second));
			negLogLik -= log(full_likelihood(pt+n));
		}
	}else if(method=="pruned"){
		// computes the pruned alignment likelihood (thanks to the forward+backward trick)
		#pragma omp parallel for reduction (+:negLogLik)
		for(n=0;n<N_seqs;n++){
			if((pt+n) -> mu_err > 0. or (pt+n) -> mutated == false){
				negLogLik -= log(partial_likelihood_new_method(pt+n));
			}else{
				negLogLik -= 0.;
			}
		}
	}
	
	return negLogLik;
}

#if GAMMA
	void init_gamma_posterior(struct seq_align *pt){
		
		int i, j, k, max_i, Lx, N_err, N_parents;
		double mode, int_mu_shift, alpha, beta;
		double mu_1, mu_2, mu_3, f_1, f_2, f_3, max_f;
		double a, b, c1, c2, c3;
		double m, curv;
		double x, p_mut, p_unmut;
		string where;
		vector <int> index_vec;
		vector <double> mu_vec, f_vec;
		
		N_err = pt -> N_err;
		N_parents = pt -> N_parents;
		Lx = (pt -> V_end) - (pt -> V_start);
		
		if(pt -> mutated == false){
			// Posterior can be approximated by a Gamma with alpha = 1
			
			// I have to estimate the decay of the likelihood
			mu_1 = prior_min + d_prior;
			pt -> mu_err = mu_1;
			f_1 = log(partial_likelihood_new_method(pt));
			
			mu_2 = prior_min + 2*d_prior;
			pt -> mu_err = mu_2;
			f_2 = log(partial_likelihood_new_method(pt));
			
			// we recover the parameters for the new gamma distribution
			pt -> alpha = 1.;
			pt -> beta = -(f_2-f_1)/(mu_2-mu_1);
			
			if(pt -> beta <= 0.){
				cout << endl;
				cout << "We have a problem with init of posterior for seq " << pt -> seq_ID << " (mu_hat_err=0)." << endl;
				cout << "We got an unexpected value for beta: " << pt -> beta << endl;
				exit(EXIT_FAILURE);
			}
			
			// we set the sequence-specific value for p_mutated
			x = d_prior;
			i = (int)((x-prior_min)/d_prior);
			pt -> mu_err = x;
			p_mut = p_mutated * partial_likelihood_new_method(pt) * prior.at(i) / gammapdf(x, pt -> alpha, pt -> beta);
			pt -> mu_err = 0.;
			p_unmut = (1. - p_mutated) * partial_likelihood_new_method(pt);
			p_mut = p_mut / (p_mut + p_unmut);
			pt -> p_mutated = p_mut;
			
		}else if(pt -> mutated == true){
			// Posterior can be approximated by a Gamma with alpha > 1
			
			mode = pt -> mu_err;  // includes also deletions and insertions from greedy alignment
			
			i = (int)((mode-prior_min)/d_prior);
			int_mu_shift = min(50,max(4,(int)(0.10*i)));
			
			index_vec.clear();
			mu_vec.clear();
			f_vec.clear();
			
			k = 2;
			// we choose 2*k+1 values of mu_s
			// then, we look for the maximum
			
			if(i-k*int_mu_shift<2){
				index_vec.push_back(2);
				for(j=1;j<2*k+1;j++){
					index_vec.push_back(2+j*int_mu_shift);
				}
			}else if(i+k*int_mu_shift>N_prior-1){
				for(j=0;j<2*k+1;j++){
					index_vec.push_back(N_prior-1-(2*k-j)*int_mu_shift);
				}
			}else{
				for(j=0;j<2*k+1;j++){
					index_vec.push_back(i+(j-2)*int_mu_shift);
				}
			}
			
			for(i=0;i<2*k+1;i++){
				if(index_vec.at(i)<0 or index_vec.at(i)>N_prior-1){
					cout << "Problem with the prior index i=" << index_vec.at(i) << endl;
					exit(EXIT_FAILURE);
				}
				mu_vec.push_back(prior_min + index_vec.at(i)*d_prior);
				pt -> mu_err = mu_vec.at(i);
				f_vec.push_back(log(partial_likelihood_new_method(pt)));
				if(i==0){
					max_i = 0;
					max_f = f_vec.at(0);
				}else{
					if(f_vec.at(i) > max_f){
						max_i = i;
						max_f = f_vec.at(i);
					}
				}
			}
			
			if(max_i>0 and max_i<2*k){
				mu_1 = mu_vec.at(max_i-1);
				f_1 = f_vec.at(max_i-1);
				mu_2 = mu_vec.at(max_i);
				f_2 = f_vec.at(max_i);
				mu_3 = mu_vec.at(max_i+1);
				f_3 = f_vec.at(max_i+1);
			}else if(max_i==0){
				where = "left";
				j = index_vec.at(0);
				mu_1 = mu_vec.at(0);
				f_1 = f_vec.at(0);
				mu_2 = mu_vec.at(1);
				f_2 = f_vec.at(1);
				mu_3 = mu_vec.at(2);
				f_3 = f_vec.at(2);
				while(where=="left"){
					j -= int_mu_shift;
					if(j>0){
						mu_3 = mu_2;
						f_3 = f_2;
						mu_2 = mu_1;
						f_2 = f_1;
						mu_1 = prior_min + j*d_prior;
						pt -> mu_err = mu_1;
						f_1 = log(partial_likelihood_new_method(pt));
						if(f_1<f_2){
							where = "stop";
						}
					}else{
						where = "stop";
					}
				}
			}else if(max_i==2*k){
				where = "right";
				j = index_vec.at(2*k);
				mu_1 = mu_vec.at(2*k-2);
				f_1 = f_vec.at(2*k-2);
				mu_2 = mu_vec.at(2*k-1);
				f_2 = f_vec.at(2*k-1);
				mu_3 = mu_vec.at(2*k);
				f_3 = f_vec.at(2*k);
				while(where=="right"){
					j += int_mu_shift;
					if(j<N_prior){
						mu_1 = mu_2;
						f_1 = f_2;
						mu_2 = mu_3;
						f_2 = f_3;
						mu_3 = prior_min + j*d_prior;
						pt -> mu_err = mu_3;
						f_3 = log(partial_likelihood_new_method(pt));
						if(f_3<f_2){
							where = "stop";
						}
					}else{
						where = "stop";
					}
				}
			}
			
			// we interpolate the parabola in order to find the maximum
			c1 = f_1/(mu_1-mu_2)/(mu_1-mu_3);
			c2 = f_2/(mu_2-mu_1)/(mu_2-mu_3);
			c3 = f_3/(mu_3-mu_1)/(mu_3-mu_2);
			a = c1 + c2 + c3;
			b = - c1*(mu_2+mu_3) - c2*(mu_1+mu_3) - c3*(mu_1+mu_2);
			//c = c1*mu_2*mu_3 + c2*mu_1*mu_3 + c3*mu_1*mu_2;
			m = -0.5*b/a;
			curv = 2*a;
			
			// we use directly the gamma interpolation
			alpha = 1. + ((mu_3-mu_2)*(f_2-f_1)-(mu_2-mu_1)*(f_3-f_2)) / ((mu_3-mu_2)*(log(mu_2)-log(mu_1))-(mu_2-mu_1)*(log(mu_3)-log(mu_2)));
			if(alpha > 1.01){
				beta = ((alpha-1.)*(log(mu_3)-log(mu_2))-(f_3-f_2)) / (mu_3-mu_2);
			}else{
				cout << "A posterior with mode>0 is trying to be initialized with alpha=1." << endl;
				alpha = 1.;
				beta = - (f_3-f_2) / (mu_3-mu_2);
			}
			
			// we recover the parameters for the new gamma distribution
			pt -> alpha = alpha;
			pt -> beta = beta;
			
			if(pt -> alpha < 1. or pt -> beta <= 0.){
				cout << endl;
				cout << "We have a problem with init of posterior for seq " << pt -> seq_ID << " (mu_hat_err>0)." << endl;
				cout << "Lx=" << Lx << ", N_err=" << pt -> N_err << ", N_del=" << pt -> N_del << ", N_ins=" << pt -> N_ins << endl;
				cout << "mode=" << mode << endl;
				cout << "mu_1=" << mu_1 << ", mu_2=" << mu_2 << ", mu_3=" << mu_3 << endl;
				cout << "f_1=" << f_1 << ", f_2=" << f_2 << ", f_3=" << f_3 << endl;
				cout << "m=" << m << ", curv=" << curv << endl;
				cout << "new_alpha=" << pt -> alpha << ", new_beta=" << pt -> beta << endl << endl;
				exit(EXIT_FAILURE);
			}
			
			// we set the sequence-specific value for p_mutated
			x = ((pt -> alpha) - 1.) / (pt -> beta);
			i = (int)((x-prior_min)/d_prior);
			pt -> mu_err = x;
			p_mut = p_mutated * partial_likelihood_new_method(pt) * prior.at(i) / gammapdf(x, pt -> alpha, pt -> beta);
			pt -> mu_err = 0.;
			p_unmut = (1. - p_mutated) * partial_likelihood_new_method(pt);
			p_mut = p_mut / (p_mut + p_unmut);
			pt -> p_mutated = p_mut;
			
		}
		
		return;
	}
#endif

#if GAMMA
	double update_gamma_posterior(struct seq_align *pt, int t){
		
		int i, j, k, max_i, Lx;
		double eta, omega, mode, alpha, old_alpha, beta, old_beta, int_mu_shift;
		double x, p_mut, p_unmut;
		double mu_1, mu_2, mu_3, lik, lik_1, lik_2, prior_1, prior_2, max_f, f_1, f_2, f_3;
		double a, b, c1, c2, c3;
		double m, curv;
		string where;
		vector <int> index_vec;
		vector <double> mu_vec, loglik_vec, logprior_vec, f_vec;
		
		old_alpha = pt -> alpha;
		old_beta = pt -> beta;
		Lx = (pt -> V_end) - (pt -> V_start);
		
		omega = 2. / (1 + exp(-(double)t/20.)) - 1.;
		eta = 0.9;  // eta=1 means no update, eta=0 means full update
		
		if(pt -> mutated == false){
			// Posterior can be approximated by a Gamma with alpha = 1
			
			// I have to estimate the decay of the likelihood times the prior
			mu_1 = prior_min + d_prior;
			pt -> mu_err = mu_1;
			lik_1 = log(partial_likelihood_new_method(pt));
			prior_1 = prior.at(1);
			f_1 = lik_1 + omega*log(prior_1);
			
			mu_2 = prior_min + 2*d_prior;
			pt -> mu_err = mu_2;
			lik_2 = log(partial_likelihood_new_method(pt));
			prior_2 = prior.at(2);
			f_2 = lik_2 + omega*log(prior_2);
			
			// we recover the parameters for the new gamma distribution
			pt -> alpha = 1.;
			pt -> beta = eta*old_beta + (1.-eta)*(-(f_2-f_1)/d_prior);
			
			if(pt -> beta <= 0.){
				cout << endl << "We have a problem with update of posterior for seq " << pt -> seq_ID << " (alpha=1 staying alpha=1)." << endl;
				cout << "loglik_1=" << lik_1 << ", loglik_2=" << lik_2 << endl;
				cout << "logprior_1=" << log(prior_1) << ", logprior_2=" << log(prior_2) << endl;
				cout << "f_1=" << f_1 << ", f_2=" << f_2 << endl;
				cout << "We got an unexpected value for beta: " << pt -> beta << endl << endl;
				exit(EXIT_FAILURE);
			}
			
			// we set the sequence-specific value for p_mutated
			x = d_prior;
			i = (int)((x-prior_min)/d_prior);
			pt -> mu_err = x;
			p_mut = p_mutated * partial_likelihood_new_method(pt) * prior.at(i) / gammapdf(x, pt -> alpha, pt -> beta);
			pt -> mu_err = 0.;
			p_unmut = (1. - p_mutated) * partial_likelihood_new_method(pt);
			p_mut = p_mut / (p_mut + p_unmut);
			pt -> p_mutated = p_mut;
			
		}else if(pt -> mutated == true){
			// Posterior can be approximated by a Gamma with alpha > 1
			
			mode = (old_alpha - 1.) / (old_beta);
			
			i = (int)((mode-prior_min)/d_prior);
			int_mu_shift = min(50,max(4,(int)(0.10*i)));
			
			index_vec.clear();
			mu_vec.clear();
			loglik_vec.clear();
			logprior_vec.clear();
			f_vec.clear();
			
			k = 2;
			// we choose 2*k+1 values of mu_s
			// then, we look for the maximum
			
			if(i-k*int_mu_shift<2){
				index_vec.push_back(2);
				for(j=1;j<2*k+1;j++){
					index_vec.push_back(2+j*int_mu_shift);
				}
			}else if(i+k*int_mu_shift>N_prior-1){
				for(j=0;j<2*k+1;j++){
					index_vec.push_back(N_prior-1-(2*k-j)*int_mu_shift);
				}
			}else{
				for(j=0;j<2*k+1;j++){
					index_vec.push_back(i+(j-2)*int_mu_shift);
				}
			}
			
			for(i=0;i<2*k+1;i++){
				if(index_vec.at(i)<0 or index_vec.at(i)>N_prior-1){
					cout << "Problem with the prior index i=" << index_vec.at(i) << endl;
					exit(EXIT_FAILURE);
				}
				mu_vec.push_back(prior_min + index_vec.at(i)*d_prior);
				pt -> mu_err = mu_vec.at(i);
				loglik_vec.push_back(log(partial_likelihood_new_method(pt)));
				logprior_vec.push_back(log(prior.at(index_vec.at(i))));
				f_vec.push_back(loglik_vec.at(i) + omega*logprior_vec.at(i));
				if(i==0){
					max_i = 0;
					max_f = f_vec.at(0);
				}else{
					if(f_vec.at(i) > max_f){
						max_i = i;
						max_f = f_vec.at(i);
					}
				}
			}
			
			if(max_i>0 and max_i<2*k){
				mu_1 = mu_vec.at(max_i-1);
				f_1 = f_vec.at(max_i-1);
				mu_2 = mu_vec.at(max_i);
				f_2 = f_vec.at(max_i);
				mu_3 = mu_vec.at(max_i+1);
				f_3 = f_vec.at(max_i+1);
			}else if(max_i==0){
				where = "left";
				j = index_vec.at(0);
				mu_1 = mu_vec.at(0);
				f_1 = f_vec.at(0);
				mu_2 = mu_vec.at(1);
				f_2 = f_vec.at(1);
				mu_3 = mu_vec.at(2);
				f_3 = f_vec.at(2);
				while(where=="left"){
					j -= int_mu_shift;
					if(j>0){
						mu_3 = mu_2;
						f_3 = f_2;
						mu_2 = mu_1;
						f_2 = f_1;
						mu_1 = prior_min + j*d_prior;
						pt -> mu_err = mu_1;
						f_1 = log(partial_likelihood_new_method(pt)) + omega*log(prior.at(j));
						if(f_1<f_2){
							where = "stop";
						}
					}else{
						where = "stop";
					}
				}
			}else if(max_i==2*k){
				where = "right";
				j = index_vec.at(2*k);
				mu_1 = mu_vec.at(2*k-2);
				f_1 = f_vec.at(2*k-2);
				mu_2 = mu_vec.at(2*k-1);
				f_2 = f_vec.at(2*k-1);
				mu_3 = mu_vec.at(2*k);
				f_3 = f_vec.at(2*k);
				while(where=="right"){
					j += int_mu_shift;
					if(j<N_prior){
						mu_1 = mu_2;
						f_1 = f_2;
						mu_2 = mu_3;
						f_2 = f_3;
						mu_3 = prior_min + j*d_prior;
						pt -> mu_err = mu_3;
						f_3 = log(partial_likelihood_new_method(pt)) + omega*log(prior.at(j));
						if(f_3<f_2){
							where = "stop";
						}
					}else{
						where = "stop";
					}
				}
			}
			
			// we interpolate the parabola
			c1 = f_1/(mu_1-mu_2)/(mu_1-mu_3);
			c2 = f_2/(mu_2-mu_1)/(mu_2-mu_3);
			c3 = f_3/(mu_3-mu_1)/(mu_3-mu_2);
			a = c1 + c2 + c3;
			b = - c1*(mu_2+mu_3) - c2*(mu_1+mu_3) - c3*(mu_1+mu_2);
			//c = c1*mu_2*mu_3 + c2*mu_1*mu_3 + c3*mu_1*mu_2;
			m = -0.5*b/a;
			curv = 2*a;
			
			// we use directly the gamma interpolation
			alpha = 1. + ((mu_3-mu_2)*(f_2-f_1)-(mu_2-mu_1)*(f_3-f_2)) / ((mu_3-mu_2)*(log(mu_2)-log(mu_1))-(mu_2-mu_1)*(log(mu_3)-log(mu_2)));
			if(alpha > 1.01){
				beta = ((alpha-1.)*(log(mu_3)-log(mu_2))-(f_3-f_2)) / (mu_3-mu_2);
			}else{
				alpha = 1.;
				beta = - (f_3-f_2) / (mu_3-mu_2);
			}
			
			// we recover the parameters for the new gamma distribution
			pt -> alpha = eta*old_alpha + (1.-eta)*alpha;
			pt -> beta = eta*old_beta + (1.-eta)*beta;
			
			if(pt -> alpha < 1. or pt -> beta <= 0.){
				cout << endl << "We have a problem with update of posterior for seq " << pt -> seq_ID << " (alpha>1 staying alpha>1)." << endl;
				cout << "N_err=" << pt -> N_err << ", N_del=" << pt -> N_del << ", N_ins=" << pt -> N_ins << endl;
				cout << "mode=" << mode << ", old_alpha=" << old_alpha << ", old_beta=" << old_beta << endl;
				cout << "mu_1=" << mu_1 << ", mu_2=" << mu_2 << ", mu_3=" << mu_3 << endl;
				//cout << "lik_1=" << lik_1 << ", lik_2=" << lik_2 << ", lik_3=" << lik_3 << endl;
				//cout << "prior_1=" << prior_1 << ", prior_2=" << prior_2 << ", prior_3=" << prior_3 << endl;
				cout << "f_1=" << f_1 << ", f_2=" << f_2 << ", f_3=" << f_3 << endl;
				cout << "m=" << m << ", curv=" << curv << endl;
				cout << "new_alpha=" << pt -> alpha << ", new_beta=" << pt -> beta << endl << endl;
				exit(EXIT_FAILURE);
			}
			
			// we set the sequence-specific value for p_mutated
			x = ((pt -> alpha) - 1.) / (pt -> beta);
			i = (int)((x-prior_min)/d_prior);
			pt -> mu_err = x;
			p_mut = p_mutated * partial_likelihood_new_method(pt) * prior.at(i) / gammapdf(x, pt -> alpha, pt -> beta);
			pt -> mu_err = 0.;
			p_unmut = (1. - p_mutated) * partial_likelihood_new_method(pt);
			p_mut = p_mut / (p_mut + p_unmut);
			pt -> p_mutated = p_mut;
			
			// return a quantity that is not interesting per se, but useful to check if things are going on properly
			return f_1+f_2+f_3;
			
		}
		
	}
#endif

#if SPLINE
	void init_spline_posterior(struct seq_align *pt){
		
		int i, Lx;
		double mode, sigma, x, d_left, d_right;
		vector <double> bins, posterior;
		
		Lx = (pt -> V_end) - (pt -> V_start);
		
		if(pt -> mutated == false){
			
			// Posterior is supposed to have a max in zero, or very close to it,
			// and then to rapidly decay
			
			bins.clear();
			mode = 0.;
			sigma = 1./Lx;
			d_right = 0.667*sigma;
			bins.push_back(mode);
			for(i=1;i<=W;i++){
			x = mode + i*d_right * (1. + 0.5*(i-1)*0.25);
				if(x>prior_min and x<prior_max){
					bins.push_back(x);
				}
			}
			
			for(i=0;i<=5;i++){
				x = 0.2*i;
				bins.push_back(x);
			}
			
			pt -> p_mutated = 0.5;
			
		}else{
			
			// Posterior is supposed to have a max at some mu_err larger than zero,
			// so both left and right tails have to be taken into account
			
			mode = (double)(pt->N_err+pt->N_del+pt->N_ins)/Lx;
			sigma = sqrt(mode*(1.-mode)/Lx);
			d_left = 0.333*sigma;
			d_right = 0.667*sigma;
			
			bins.push_back(mode);
			for(i=1;i<=W;i++){
				x = mode - i*d_left * (1+0.5*(i-1)*0.5);
				if(x>prior_min and x<prior_max){
					bins.push_back(x);
				}
				x = mode + i*d_right * (1+0.5*(i-1)*0.5);
				if(x>prior_min and x<prior_max){
					bins.push_back(x);
				}
			}
			
			for(i=0;i<=5;i++){
				x = 0.2*i;
				bins.push_back(x);
			}
			
			pt -> p_mutated = 1.0;
			
		}
		
		pt -> mode = mode;
		pt -> mean = mode;
		pt -> sigma = sigma;
		(pt -> sampled_mu_errs).clear();
		(pt -> sampled_mu_errs).push_back(pt -> mu_err);
		
		stable_sort(bins.begin(),bins.end(),sort_using_smaller_than);
		bins.erase(unique(bins.begin(),bins.end()),bins.end());
		
		pt -> spline_bins = bins;
		posterior.clear();
		for(i=0;i<(int)bins.size();i++){
			// Has to be properly initialized with info only coming from Lik
			// -> can actually be done by calling update_spline_posterior() with t=0 !
			posterior.push_back(1./(int)bins.size());
		}
		pt -> spline_posterior = posterior;
		
		return;
		
	}
#endif

#if SPLINE
	void update_spline_posterior(struct seq_align *pt, int t){
		
		int i, k;
		double omega;
		double mu, x, lik, pr, f, S, Sp, mode_x, mode_y, mean, mean2, sigma, d_left, d_right, p_unmut, p_mut;
		
		vector <double> bins, posterior, post_pdf, post_cdf;
		
		bins = pt -> spline_bins;
		
		omega = 2. / (1 + exp(-(double)t/20.)) - 1.;
		
		posterior.clear();
		posterior.resize((int)bins.size(),0.);
		for(i=0;i<(int)bins.size();i++){
			mu = bins.at(i);
			pt -> mu_err = mu;
			if(mu==0.){
				if(pt -> mutated == true){
					lik = 0.;
				}else{
					lik = partial_likelihood_new_method(pt);
				}
				pr = prior.at(0);
			}else if(mu==1.){
				lik = 0.;
				pr = prior.at(N_prior-1);
			}else if(mu>0. and mu<1.){
				lik = partial_likelihood_new_method(pt);
				k = (int)((mu-prior_min)/d_prior);
				pr = prior.at(k) + (mu-prior_min-k*d_prior) * (prior.at(k+1)-prior.at(k))/d_prior;
			}
			f = (1.-omega)*lik + omega*lik*pr;
			posterior.at(i) = f;
		}
		
		post_pdf.clear();
		tk::spline s;
		//s.set_boundary(tk::spline::second_deriv, 0.0, tk::spline::first_deriv, 0.0);
		s.set_points(bins,posterior);
		s.make_monotonic();
		S = 0.;
		post_cdf.push_back(0.);
		for(k=0;k<N_prior;k++){
			x = max(0.,s(prior_min+k*d_prior));
			post_pdf.push_back(x);
			if(k==0){
				mode_x = 0.;
				mode_y = post_pdf.at(0);
			}else{
				if(post_pdf.at(k)>mode_y){
					mode_x = prior_min+k*d_prior;
					mode_y = post_pdf.at(k);
				}
				post_cdf.push_back(post_cdf.at(k-1)+x);
			}
			S += x;
		}
		S -= 0.5*post_pdf.at(0);
		S *= d_prior;
		Sp = post_cdf.at(N_prior-1);
		mean = 0.;
		mean2 = 0.;
		for(k=0;k<N_prior;k++){
			post_pdf.at(k) /= S;
			post_cdf.at(k) /= Sp;
			x = prior_min+k*d_prior;
			mean += x * post_pdf.at(k);
			mean2 += x * x * post_pdf.at(k);
		}
		mean *= d_prior;
		mean2 *= d_prior;
		sigma = sqrt(mean2-mean*mean);
		for(i=0;i<(int)bins.size();i++){
			posterior.at(i) /= S;
		}
		
		pt -> mode = mode_x;
		pt -> mean = mean;
		pt -> sigma = sigma;
		pt -> spline_posterior = posterior;
		pt -> Z_prime = S;
		
		// sample from continuous posterior
		(pt -> sampled_mu_errs).clear();
		for(i=0;i<N_sampled_mu_errs;i++){
			discrete_distribution<int> post_distr(post_pdf.begin(), post_pdf.end());
			do{
				k = post_distr(generator);
				mu = prior_min + k*d_prior;
				mu += (frand-0.5)*d_prior;
			}while(mu<prior_min or mu>prior_max);
			(pt -> sampled_mu_errs).push_back(mu);
		}
		
		// refresh bins
		if(t>0 and t%5==0 and pt -> mutated == true){
			mode_x = pt -> mode;
			sigma = pt -> sigma;
			d_left = 0.333*sigma;
			d_right = 0.667*sigma;
			
			bins.clear();
			bins.push_back(mode_x);
			for(i=1;i<=W;i++){
				x = mode_x - i*d_left * (1+0.5*(i-1)*0.5);
				if(x>prior_min and x<prior_max){
					bins.push_back(x);
				}
				x = mode_x + i*d_right * (1+0.5*(i-1)*0.5);
				if(x>prior_min and x<prior_max){
					bins.push_back(x);
				}
			}
			
			for(i=0;i<=5;i++){
				x = 0.2*i;
				bins.push_back(x);
			}
			
			stable_sort(bins.begin(),bins.end(),sort_using_smaller_than);
			bins.erase(unique(bins.begin(),bins.end()),bins.end());
			
			pt -> spline_bins = bins;
			
		}
		
		// update p_mutated
		if(pt -> mutated == false){
			pt -> mu_err = 0.;
			lik = partial_likelihood_new_method(pt);
			p_unmut = (1.-p_mutated) * lik;
			p_mut = p_mutated * (pt -> Z_prime);
			pt -> p_mutated = p_mut / (p_unmut + p_mut);
		}else{
			pt -> p_mutated = 1.0;
		}
		
		// update the prior
		for(k=0;k<N_prior;k++){
			new_prior[k] += (pt->p_mutated) * post_pdf.at(k);
		}
		
		return;
		
	}
#endif

vector <double> smooth_prior(vector <double> prior, double prior_mean, string method){
	int i, j, k, delta_j, delta_k, N_y;
	double S, s, sigma, x, mean, y, min_y, d_y, jac;
	vector <double> prior_temp, prior_smoothed, prior_y;
	//vector <double> prior_y_smoothed;
	double *prior_y_array;
	double *prior_y_smoothed_array;
	double *norm_array;
	
	if(method=="Gaussian_KDE"){
		// Use Gaussian KDE in the "linear" space
		
		double xx = Gauss_KDE_sigma;
		double ss = Gauss_KDE_sigma;
		
		delta_j = (int)(5.*ss / d_prior);
		
		double prior_array[N_prior], prior_smoothed_array[N_prior], norm_array[N_prior], sigma[N_prior];
		
		for(i=0;i<N_prior;i++){
			x = prior_min + i*d_prior;
			prior_array[i] = prior.at(i);
			prior_smoothed_array[i] = 0.;
			norm_array[i] = 0.;
			if(x<=xx){
				sigma[i] = ss * exp(-0.5*(x-xx)*(x-xx)/ss/ss);  // Gaussian decay
				//sigma[i] = 0.5*ss + x*(0.5*ss/xx);  // Linear decay
			}else if(x>xx and x<prior_max-xx){
				sigma[i] = ss;
			}else if(x>=prior_max-xx){
				sigma[i] = ss * exp(-0.5*(x-prior_max+xx)*(x-prior_max+xx)/ss/ss);  // Gaussian decay
				//sigma[i] = ss - (x-prior_max+xx)*(0.5*ss/xx);  // Linear decay
			}
			sigma[i] = ss;
		}
		
		#pragma omp parallel for private(j) shared(sigma) reduction(+:prior_smoothed_array,norm_array)
		for (i=0;i<N_prior;i++){
			for(j=max(0,i-delta_j);j<min(N_prior,i+delta_j);j++){
				prior_smoothed_array[i] += prior_array[j] * gausspdf(prior_min + i*d_prior,prior_min+j*d_prior,sigma[i]);
				norm_array[i] += gausspdf(prior_min + i*d_prior,prior_min+j*d_prior,sigma[i]);
			}
		}
		
		prior_temp.clear();
		S = 0.;
		for(i=0;i<N_prior;i++){
			prior_smoothed_array[i] /= norm_array[i];
			S += prior_smoothed_array[i];
		}
		S -= 0.5*prior_smoothed_array[0];
		S *= d_prior;
		
		
		/////
		// here a fix for close-to-zero region
		for(i=0;i<N_prior;i++){
			prior_smoothed_array[i] /= S;
		}
		tk::spline s;
		vector <double> nodes {0.003,0.004,0.005}; // {3.*ss,4.*ss,5.*ss}
		vector <double> values;
		values.clear();
		for(i=0;i<(int)nodes.size();i++){
			x = nodes.at(i);
			values.push_back(prior_smoothed_array[(int)((x-prior_min)/d_prior)]);
		}
		//s.set_boundary(tk::spline::second_deriv, 0.0, tk::spline::first_deriv, 0.0);
		s.set_points(nodes,values);
		s.make_monotonic();
		int k_max = (int)((nodes.at(0)-prior_min)/d_prior);
		S = 0.;
		for(k=0;k<N_prior;k++){
			if(k<k_max){
				prior_smoothed_array[k] = max(0.,s(prior_min+k*d_prior));
			}
			S += prior_smoothed_array[k];
		}
		S -= 0.5*prior_smoothed_array[0];
		S *= d_prior;
		/////
		
		
		prior_temp.clear();
		for(i=0;i<N_prior;i++){
			prior_temp.push_back(prior_smoothed_array[i]/S);
		}
		
	}else if(method=="Gaussian_KDE_log"){
		// Use Gaussian KDE in the log space
		// y = log(x)
		
		sigma = Gauss_KDE_sigma;
		
		min_y = -8.;
		d_y = 0.005;
		
		// Build f(y)
		prior_y.clear();
		y = min_y;
		x = exp(y);
		i = max(0,(int)((x-prior_min)/d_prior));
		while(i < N_prior){
			jac = 1./x;
			prior_y.push_back(prior.at(i)/jac);
			y += d_y;
			x = exp(y);
			i = max(0,(int)((x-prior_min)/d_prior));
		}
		N_y = (int)(prior_y.size());
		
		prior_y_array = new double[N_y];
		for(j=0;j<N_y;j++){
			prior_y_array[j] = prior_y.at(j);
		}
		
		// Smooth it
		prior_y_smoothed_array = new double[N_y];
		norm_array = new double[N_y];
		delta_k = (int)(5.*sigma/d_y);
		
		#pragma omp parallel for private(j,k) shared(prior_y_array,sigma) reduction(+:prior_y_smoothed_array[0:N_y],norm_array[0:N_y])
		for(j=0;j<N_y;j++){
			for(k=max(0,j-delta_k);k<min(N_y,j+delta_k);k++){
				prior_y_smoothed_array[j] += prior_y_array[k] * gausspdf(min_y+j*d_y,min_y+k*d_y,sigma);
				norm_array[j] += gausspdf(min_y+j*d_y,min_y+k*d_y,sigma);
			}
		}
		
		S = 0.;
		for(j=0;j<N_y;j++){
			prior_y_smoothed_array[j] /= norm_array[j];
			S += prior_y_smoothed_array[j];
		}
		S *= d_y;
		
		// Back to linear space
		prior_temp.clear();
		S = 0.;
		for(i=0;i<N_prior;i++){
			if(i==0){
				x = prior_min + d_prior; // placeholder to avoid issues when x = 0
			}else{
				x = prior_min + i*d_prior;
			}
			y = log(x);
			j = (int)((y-min_y+0.0001*d_y)/d_y);
			jac = 1./x;
			s = prior_y_smoothed_array[j] + (y - (min_y+j*d_y)) * (prior_y_smoothed_array[j+1]-prior_y_smoothed_array[j]) / d_y; // interpolation over prior_y_smoothed
			s *= jac;
			s = max(s,1e-30);
			prior_temp.push_back(s);
			S += prior_temp.at(i);
		}
		// Here we fix the first component
		S -= prior_temp.at(0);
		prior_temp.at(0) = max(0.,2.*prior_temp.at(1) - prior_temp.at(2));
		S += prior_temp.at(0);
		S *= d_prior;
		for(i=0;i<N_prior;i++){
			prior_temp.at(i) /= S;
		}
		
	}else if(method=="Gaussian_KDE_log_rescaled"){
		// Use Gaussian KDE in the rescaled log space
		// y = log(x/mean)
		
		mean = prior_mean;
		sigma = Gauss_KDE_sigma;
		
		d_y = 0.005;
		min_y = (int)(log((prior_min+d_prior)/mean)-1.) - 4.;
		
		// Build f(y)
		prior_y.clear();
		y = min_y;
		x = mean*exp(y);
		i = max(0,(int)((x-prior_min)/d_prior));
		while(i < N_prior){
			jac = 1./x;
			prior_y.push_back(prior.at(i)/jac);
			y += d_y;
			x = mean*exp(y);
			i = max(0,(int)((x-prior_min)/d_prior));
		}
		N_y = (int)(prior_y.size());
		
		prior_y_array = new double[N_y];
		for(j=0;j<N_y;j++){
			prior_y_array[j] = prior_y.at(j);
		}
		
		// Smooth it
		prior_y_smoothed_array = new double[N_y];
		norm_array = new double[N_y];
		delta_k = (int)(5.*sigma/d_y);
		
#pragma omp parallel for private(j,k) shared(prior_y_array,sigma) reduction(+:prior_y_smoothed_array[0:N_y],norm_array[0:N_y])
		for(j=0;j<N_y;j++){
			for(k=max(0,j-delta_k);k<min(N_y,j+delta_k);k++){
				prior_y_smoothed_array[j] += prior_y_array[k] * gausspdf(min_y+j*d_y,min_y+k*d_y,sigma);
				norm_array[j] += gausspdf(min_y+j*d_y,min_y+k*d_y,sigma);
			}
		}
		
		S = 0.;
		for(j=0;j<N_y;j++){
			prior_y_smoothed_array[j] /= norm_array[j];
			S += prior_y_smoothed_array[j];
		}
		S *= d_y;
		
		// Back to linear space
		prior_temp.clear();
		S = 0.;
		for(i=0;i<N_prior;i++){
			if(i==0){
				x = prior_min + d_prior; // placeholder to avoid issues when x = 0
			}else{
				x = prior_min + i*d_prior;
			}
			y = log(x/mean);
			j = (int)((y-min_y+0.0001*d_y)/d_y);
			jac = 1./x;
			s = prior_y_smoothed_array[j] + (y - (min_y+j*d_y)) * (prior_y_smoothed_array[j+1]-prior_y_smoothed_array[j]) / d_y; // interpolation over prior_y_smoothed
			s *= jac;
			s = max(s,1e-30);
			prior_temp.push_back(s);
			S += prior_temp.at(i);
		}
		// Here we fix the first component
		S -= prior_temp.at(0);
		prior_temp.at(0) = max(0.,2.*prior_temp.at(1) - prior_temp.at(2));
		S += prior_temp.at(0);
		S *= d_prior;
		for(i=0;i<N_prior;i++){
			prior_temp.at(i) /= S;
		}
		
	}else if(method=="Gaussian_KDE_soft_log"){
		// Use Gaussian KDE in the soft-log space
		// y = log(exp(x/mean)-1)
		
		mean = prior_mean;
		sigma = Gauss_KDE_sigma / mean; // works properly from y~y_min to y~y_max
		
		// d_y given by d_x goes from ~log(i+1)-log(i) for x~x_min to ~d_prior/mean for x~x_max
		// so 0.001 should always be ok, as d_x=0.0005 and mean<<1
		d_y = 0.001;
		min_y = (int)(log(exp((prior_min+d_prior)/mean)-1.)) - 4.;
		
		// Build f(y)
		prior_y.clear();
		y = min_y;
		x = mean*log(exp(y)+1.);
		i = (int)((x-prior_min)/d_prior);
		while(i < N_prior){
			jac = exp(x/mean)/mean/(exp(x/mean)-1.);
			prior_y.push_back(prior.at(i)/jac);
			y += d_y;
			x = mean*log(exp(y)+1.);
			i = (int)((x-prior_min)/d_prior);
		}
		N_y = (int)(prior_y.size());
		
		prior_y_array = new double[N_y];
		for(j=0;j<N_y;j++){
			prior_y_array[j] = prior_y.at(j);
		}
		
		// Smooth it
		prior_y_smoothed_array = new double[N_y];
		norm_array = new double[N_y];
		delta_k = (int)(5.*sigma/d_y);
		
		#pragma omp parallel for private(j,k) shared(prior_y_array,sigma) reduction(+:prior_y_smoothed_array[0:N_y],norm_array[0:N_y])
		for(j=0;j<N_y;j++){
			for(k=max(0,j-delta_k);k<min(N_y,j+delta_k);k++){
				prior_y_smoothed_array[j] += prior_y_array[k] * gausspdf(min_y+j*d_y,min_y+k*d_y,sigma);
				norm_array[j] += gausspdf(min_y+j*d_y,min_y+k*d_y,sigma);
			}
		}
		
		S = 0.;
		for(j=0;j<N_y;j++){
			prior_y_smoothed_array[j] /= norm_array[j];
			S += prior_y_smoothed_array[j];
		}
		S *= d_y;
		
		// Back to linear space
		prior_temp.clear();
		S = 0.;
		for(i=0;i<N_prior;i++){
			if(i==0){
				x = prior_min + d_prior; // placeholder to avoid issues when x = 0
			}else{
				x = prior_min + i*d_prior;
			}
			y = log(exp(x/mean)-1.);
			j = (int)((y-min_y+0.0001*d_y)/d_y);
			jac = exp(x/mean)/mean/(exp(x/mean)-1.);
			//s = prior_y_smoothed.at(j);
			s = prior_y_smoothed_array[j] + (y - (min_y+j*d_y)) * (prior_y_smoothed_array[j+1]-prior_y_smoothed_array[j]) / d_y; // interpolation over prior_y_smoothed
			s *= jac;
			s = max(s,1e-30);
			prior_temp.push_back(s);
			S += prior_temp.at(i);
		}
		// Here we fix the first component
		S -= prior_temp.at(0);
		prior_temp.at(0) = max(0.,2.*prior_temp.at(1) - prior_temp.at(2));
		S += prior_temp.at(0);
		S *= d_prior;
		for(i=0;i<N_prior;i++){
			prior_temp.at(i) /= S;
		}
		
	}else if(method=="none"){
		prior_temp = prior;
	}
	
	return prior_temp;
}

void tune_KDE_sigma(string batch_name, int attempt_number, vector <double> prior_mean_history, vector <double> prior_std_history, int t_w){
	int i, i_min, i_max, j, dt, t_max;
	double prior_min_history, prior_max_history;
	string filename, line;
	ifstream infile;
	vector <string> vec_str;
	vector <double> prior_history;
	bool minimum, maximum, increase, reduce;
	
	t_max = (int)prior_mean_history.size() - 1;
	dt = (int)(0.8*t_w);
	
	// check if prior_std is becoming smaller and smaller
	reduce = false;
	minimum = false;
	maximum = false;
	i_min = t_max-dt;
	i_max = t_max-dt;
	prior_min_history = prior_std_history.at(i_min);
	prior_max_history = prior_std_history.at(i_max);
	for(i=t_max-dt+1;i<=t_max-1;i++){
		if(prior_std_history.at(i)<prior_min_history){
			i_min = i;
			prior_min_history = prior_std_history.at(i_min);
		}
		if(prior_std_history.at(i)<prior_std_history.at(i-1) and prior_std_history.at(i)<prior_std_history.at(i+1)){
			minimum = true;
		}
		if(prior_std_history.at(i)>prior_max_history){
			i_max = i;
			prior_max_history = prior_std_history.at(i_max);
		}
		if(prior_std_history.at(i)>prior_std_history.at(i-1) and prior_std_history.at(i)>prior_std_history.at(i+1)){
			maximum = true;
		}
	}
	if(minimum==false and maximum==false and i_min==t_max-1){
		if(prior_min_history < 0.95*prior_max_history){
			cout << "Gauss_KDE_sigma was too large: " << Gauss_KDE_sigma  << " vs prior_std=" << prior_std_history.at(t_max) << endl;
			Gauss_KDE_sigma *= prior_min_history/prior_max_history;
			//Gauss_KDE_sigma /= 1.2;
			cout << "It has been lowered to: " << Gauss_KDE_sigma << endl;
			reduce = true;
		}
	}
	
	// otherwise, check if a larger Gauss_KDE_sigma is needed
	if(reduce==false){
		increase = false;
		for(j=-2;j<=2;j++){
			i = (int)((prior_mean_history.at(t_max)+j*prior_std_history.at(t_max)-prior_min)/d_prior);
			if(i>0 and i<N_prior-1){
				
				filename = batch_name + "inference/aux_prior_file_attempt_" + to_string(attempt_number) + ".txt";
				infile.open(filename);
				if(!infile){
					throw runtime_error("File not found: " + filename);
				}
				
				prior_history.clear();
				getline(infile,line);  // we get rid of headers
				getline(infile,line);  // we get rid of mu spacings
				while(getline(infile,line)){
					vec_str = split_string(line, ';');
					prior_history.push_back(stod(vec_str.at(i+1)));
				}
				infile.close();
				
				minimum = false;
				maximum = false;
				i_min = t_max-dt;
				i_max = t_max-dt;
				prior_min_history = prior_history.at(i_min);
				prior_max_history = prior_history.at(i_max);
				for(i=t_max-dt+1;i<=t_max-1;i++){
					if(prior_history.at(i)<prior_min_history){
						i_min = i;
						prior_min_history = prior_history.at(i_min);
					}
					if(prior_history.at(i)<prior_history.at(i-1) and prior_history.at(i)<prior_history.at(i+1)){
						minimum = true;
					}
					if(prior_history.at(i)>prior_max_history){
						i_max = i;
						prior_max_history = prior_history.at(i_max);
					}
					if(prior_history.at(i)>prior_history.at(i-1) and prior_history.at(i)>prior_history.at(i+1)){
						maximum = true;
					}
				}
				if((prior_max_history-prior_min_history)>0.10*0.5*(prior_max_history+prior_min_history) and (minimum==true or maximum==true)){
					increase = true;
				}
			}
		}
		if(increase==true){
			cout << "Gauss_KDE_sigma was too small: " << Gauss_KDE_sigma  << " vs prior_std=" << prior_std_history.at(t_max) << endl;
			Gauss_KDE_sigma *= 1.1;
			cout << "It has been increased to: " << Gauss_KDE_sigma << endl;
		}
	}
	
	return;
}

void explore_scenarios(struct seq_align *pt){
	
	int child_i, parent_i, child_j, parent_j, k, nn, N_cells, N_children, p, N_scenarios;
	
	double logLik;
	double mu_s_err, log_cost_s_match, log_cost_s_error;
	string germline, read, alignment;
	
	read = pt -> seq;
	germline = genomic_V[pt -> V_best_idx].second.substr(pt -> V_start,(pt -> V_end)-(pt -> V_start));
	
	mu_s_err = pt -> mu_err;
	log_cost_s_match = log((1. - mu_s_err) * (1. - mu_s_err*params[1] - mu_s_err*params[2]));
	log_cost_s_error = log(mu_s_err / 3  * (1. - mu_s_err*params[1] - mu_s_err*params[2]));
	
	N_cells = (pt -> N_active_cells);
	cout << "N_cells: " << N_cells << endl;
	
	struct active_cell *new_active_cells_list;
	new_active_cells_list = new active_cell[N_cells];
	new_active_cells_list = pt -> active_cells;
	
	cout << "Building the tree..." << endl;
	
	// variables for the bulding of the tree
	struct tree_node{
		int n;  // is the integer that labels the cell into new_active_cells_list
		int parent;  // is the integer for the parent node in the tree
	};
	struct tree_node node, child_node, parent_node;
	vector <struct tree_node> tree;
	vector <int> stack;  // set of leaves up to the previous iteration (integers refer to numbering in tree vector)
	vector <int> wk_list;  // new set of leaves (integers refer to numbering in tree vector)
	vector <int> final_leaves;  // set of final leaves (integers refer to numbering in tree vector)
	
	// Root of the tree
	node.n = N_cells-1;
	node.parent = -1;  // root has no parent
	tree.push_back(node);
	stack.push_back(0);
	
	// Building the tree
	while((int)stack.size()>0){
		wk_list.clear();
		for(k=0;k<(int)stack.size();k++){
			N_children = (int)new_active_cells_list[tree.at(stack.at(k)).n].parents.size();
			for(p=0;p<N_children;p++){
				node.n = new_active_cells_list[tree.at(stack.at(k)).n].parents.at(p);
				node.parent = stack.at(k);
				tree.push_back(node);
				if(node.n==0){
					final_leaves.push_back((int)tree.size()-1);
				}else{
					wk_list.push_back((int)tree.size()-1);
				}
			}
		}
		stack.clear();
		stack = wk_list;
	}
	
	//cout << "final_leaves[0]: " << tree.at(final_leaves.at(0)).n << ", " << tree.at(final_leaves.at(0)).parent << endl;
	
	// Reconstructing all the paths
	N_scenarios = (int)final_leaves.size();
	cout << "N_scenarios: " << N_scenarios << endl;
	cout << "tree_size: " << (int)tree.size() << endl;
	
	for(nn=0;nn<min(N_scenarios,5);nn++){
		
		alignment = "";
		logLik = 0.;
		child_node = tree.at(final_leaves.at(nn));
		while(child_node.parent!=-1){
			child_i = new_active_cells_list[child_node.n].row;
			child_j = new_active_cells_list[child_node.n].col;
			parent_node = tree.at(child_node.parent);
			parent_i = new_active_cells_list[parent_node.n].row;
			parent_j = new_active_cells_list[parent_node.n].col;
			if( parent_i==child_i+1 and parent_j==child_j+1 ){
				// match/mismatch
				logLik += (germline.at(child_i)==read.at(child_j) ? log_cost_s_match : log_cost_s_error);
				alignment += (germline.at(child_i)==read.at(child_j) ? "m" : "e");
			}else if( parent_i>child_i and parent_j==child_j ){
				// deletion (actual deletion lenght is child_i-parent_i)
				logLik += log_gamma_del_array[parent_i-child_i-1];
				alignment += "{d," + to_string(parent_i-child_i) + "}";
			}else if( parent_i==child_i and parent_j>child_j ){
				// insertion (actual insertion lenght is child_j-parent_j)
				logLik += log_gamma_ins_array[parent_j-child_j-1];
				alignment += "{i," + to_string(parent_j-child_j) + "}";
			}
			child_node = parent_node;
		}
		
		cout << "n_scen: " << nn << "; align: " << alignment << "; logLik: " << logLik << endl;
	}
	
}
