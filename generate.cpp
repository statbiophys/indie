/*
 *
 *  generate.cpp
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

#include "funcs.h"
#include "generate.h"

void init_synth_seqs(int N_seqs_synth, struct synth_seq *pt, unordered_map <int, string> genomic_V_map_inverse, unordered_map <string, int> genomic_V_map){
	int i, n, Lx, initial_V_cut;
	double mu_err, mode, var, alpha, beta;
	string prior_mode, V_usage_mode, line_str, prior_filename, V_usage_filename;
	vector <string> temp_vec_str;
	vector <double> prior_binning, prior_density, V_usage;
	vector <string> V_names;
	ifstream infile;
	
	prior_mode = "explicit";  // to be chosen in ['explicit','gauss','gamma']
	prior_filename = "./models/prior.csv"; // default file is "./models/prior.csv"
	V_usage_mode = "uniform";  // to be chosen in ['uniform','data']
	V_usage_filename = "./models/V_usage.csv"; // default file is "./models/V_usage.csv"
	initial_V_cut = 0;  // number of nt to cut on 5' side of V templates when generating synth seqs
	
	// The distribution from which draw the sequence-specific mutation rate can be either loaded explicitly or given analitically
	if(prior_mode=="explicit"){
		infile.open(prior_filename);
		if(!infile){
			throw runtime_error("File not found: " + prior_filename);
		}
		system(&("cp " + prior_filename + " " + batch_name + "generate/generative_mu_err_distr.txt")[0]);
		getline(infile,line_str);
		temp_vec_str.clear();
		temp_vec_str = split_string(line_str, ';');
		for(i=0;i<(int)temp_vec_str.size();i++){
			prior_binning.push_back(stod(temp_vec_str.at(i)));
		}
		getline(infile,line_str);
		temp_vec_str.clear();
		temp_vec_str = split_string(line_str, ';');
		for(i=0;i<(int)temp_vec_str.size();i++){
			prior_density.push_back(stod(temp_vec_str.at(i)));
		}
		infile.close();
	}else if(prior_mode=="gauss"){
		mode = params[0];
		var = 0.01;
	}else if(prior_mode=="gamma"){
		mode = params[0];
		var = 0.01;
		beta = (mode + sqrt(mode*mode+4*var)) / (2.*var);
		alpha = 1. + mode*beta;
	}
	
	// non-uniform V usage
	if(V_usage_mode=="data"){
		infile.open(V_usage_filename);
		if(!infile){
			throw runtime_error("File not found: " + V_usage_filename);
		}
		system(&("cp " + V_usage_filename + " " + batch_name + "generate/V_usage.csv")[0]);
		getline(infile,line_str);
		temp_vec_str.clear();
		temp_vec_str = split_string(line_str, ';');
		for(i=0;i<(int)temp_vec_str.size();i++){
			V_names.push_back(temp_vec_str.at(i));
		}
		getline(infile,line_str);
		temp_vec_str.clear();
		temp_vec_str = split_string(line_str, ';');
		for(i=0;i<(int)temp_vec_str.size();i++){
			V_usage.push_back(stod(temp_vec_str.at(i)));
		}
		infile.close();
	}
	
	for(n=0;n<N_seqs_synth;n++){
		// choose the germline template
		if(V_usage_mode=="uniform"){
			(pt+n) -> V_best_idx = rand() % (int)genomic_V.size();
			(pt+n) -> V_best = genomic_V_map_inverse[(pt+n) -> V_best_idx];
		}else if(V_usage_mode=="data"){
			discrete_distribution<int> V_distr(V_usage.begin(), V_usage.end());
			i = V_distr(generator);
			(pt+n) -> V_best = V_names.at(i) + "*01";
			(pt+n) -> V_best_idx = genomic_V_map[(pt+n) -> V_best];
		}
		(pt+n) -> V_start = initial_V_cut;
		(pt+n) -> V_end = (int)(genomic_V[(pt+n) -> V_best_idx].second.size());
		
		// draw the point mutation rate
		if(frand<frac_mutated){
			do{
				if(prior_mode=="explicit"){
					discrete_distribution<int> prior_distr(prior_density.begin(), prior_density.end());
					i = prior_distr(generator);
					mu_err = prior_binning.at(i) + (frand-0.5)*(prior_binning.at(1)-prior_binning.at(0));
				}else if(prior_mode=="gauss"){
					normal_distribution<double> gauss_distr(mode,sqrt(var));
					mu_err = gauss_distr(generator);
				}else if(prior_mode=="gamma"){
					gamma_distribution<double> gamma_distr(alpha,1./beta);
					mu_err = gamma_distr(generator);
				}
			}while(mu_err<prior_min or mu_err>prior_max);
		}else{
			mu_err = 0.;
		}
		(pt+n) -> mu_err = mu_err;
		
		// Draw the numbers of point mutations, deletions and insertions to be realized
		// Useful for type_2 generation of point mutations and indels
		Lx = (int)(genomic_V[(pt+n) -> V_best_idx].second.size()) - initial_V_cut;
		
		poisson_distribution<int> poiss_distr(mu_err*Lx);
		(pt+n) -> N_err = poiss_distr(generator);
		
		poisson_distribution<int>(mu_err*params[1]*Lx);
		(pt+n) -> N_del = poiss_distr(generator);
		
		poisson_distribution<int>(mu_err*params[2]*Lx);
		(pt+n) -> N_ins = poiss_distr(generator);
	}
	
	return;
	
}

void generate_synth_seqs_type_1(int N_seqs_synth, struct synth_seq *pt){
	
	int V_start, V_end, del_flag, del_len, ins_len, i, j, k, n;
	string germline, err, batch_name;
	vector<double> del_profile_cumul;
	vector<double> ins_profile_cumul;
	double mu_err, mu_del, mu_ins, p;
	
	// Build cumulative of length profiles
	for(i=0;i<gap_bound;i++){
		if(i==0){
			del_profile_cumul.push_back(params[del_params_begin+i]);
			ins_profile_cumul.push_back(params[ins_params_begin+i]);
		}else{
			del_profile_cumul.push_back(params[del_params_begin+i]+del_profile_cumul[i-1]);
			ins_profile_cumul.push_back(params[ins_params_begin+i]+ins_profile_cumul[i-1]);
		}
	}
	
	// Build vectors for deletion/insertion length profiles
	vector <double> p_del, p_ins;
	p_del.clear();
	p_ins.clear();
	for(i=0;i<gap_bound;i++){
		p_del.push_back(params[del_params_begin+i]);
		p_ins.push_back(params[ins_params_begin+i]);
	}
	discrete_distribution<int> del_distr(p_del.begin(), p_del.end());
	discrete_distribution<int> ins_distr(p_ins.begin(), p_ins.end());
	
	for(n=0;n<N_seqs_synth;n++){
		
		(pt+n) -> N_err = 0;
		(pt+n) -> N_del = 0;
		(pt+n) -> N_ins = 0;
		
		mu_err = (pt+n) -> mu_err;
		mu_del = mu_err * params[1];
		mu_ins = mu_err * params[2];
		
		germline = genomic_V[(pt+n) -> V_best_idx].second;
		V_start = (pt+n) -> V_start;
		V_end = (pt+n) -> V_end;
		
		(pt+n) -> error_pos = "[";
		(pt+n) -> error_list = "[";
		(pt+n) -> insertion_pos = "[";
		(pt+n) -> insertion_list = "[";
		(pt+n) -> deletion_pos = "[";
		(pt+n) -> deletion_list = "[";
		
		// insert mutations and hyper-indels on the germline
		(pt+n) -> seq = "";
		del_flag = 0;
		for(j=V_start;j<V_end;j++){
			if(del_flag==0){
				p = frand;
				if(p < mu_del){
					// Deletion
					del_len = 1 + del_distr(generator);
					// check if the number of nt to be deleted goes beyond the end of the sequence
					// this introduces a bias in favour of short deletions in the original deletion profile!!!
					del_len = min((int)germline.length()-j,del_len);
					del_flag = del_len - 1;
					if((pt+n) -> deletion_pos.length()>1){
						(pt+n) -> deletion_pos.append(",");
						(pt+n) -> deletion_list.append(",");
					}
					(pt+n) -> deletion_pos.append(to_string(j));
					(pt+n) -> deletion_list.append(to_string(del_len));
				}else if(p >= mu_del and p < mu_del+mu_ins){
					// Insertion
					ins_len = 1 + ins_distr(generator);
					err = "";
					for(k=0;k<ins_len;k++){
						err += nuclAlphabet[rand() % 4];
					}
					(pt+n) -> seq += err;
					if((pt+n) -> insertion_pos.length()>1){
						(pt+n) -> insertion_pos.append(",");
						(pt+n) -> insertion_list.append(",");
					}
					(pt+n) -> insertion_pos.append(to_string(j));
					(pt+n) -> insertion_list.append(to_string(ins_len)+"->");
					(pt+n) -> insertion_list.append(err);
					//(pt+n) -> seq += germline.at(j);
					j--;  // needed to compensate the j++ in the for loop, as insertions do not increase the index along the germline
				}else{
					// Aligned symbols
					if(frand < mu_err){
						// Point mutation
						if(germline.at(j)=='A'){
							err = notA[rand() % 3];
						}else if(germline.at(j)=='C'){
							err = notC[rand() % 3];
						}else if(germline.at(j)=='G'){
							err = notG[rand() % 3];
						}else if(germline.at(j)=='T'){
							err = notT[rand() % 3];
						}
						(pt+n) -> seq += err;
						if((pt+n) -> error_pos.length()>1){
							(pt+n) -> error_pos.append(",");
							(pt+n) -> error_list.append(",");
						}
						(pt+n) -> error_pos.append(to_string(j));
						(pt+n) -> error_list.append(1,germline.at(j));
						(pt+n) -> error_list.append("->"+err);
					}else{
						// Match
						(pt+n) -> seq += germline.at(j);
					}
				}
			}else{
				del_flag -= 1;
			}
		}
		
		(pt+n) -> error_pos.append("]");
		(pt+n) -> error_list.append("]");
		(pt+n) -> deletion_pos.append("]");
		(pt+n) -> deletion_list.append("]");
		(pt+n) -> insertion_pos.append("]");
		(pt+n) -> insertion_list.append("]");
		
	}
	
	return;
	
}

void generate_synth_seqs_type_2(int N_seqs_synth, struct synth_seq *pt){
	
	// TODO: different method for the generation of synthetic sequences, work in progress...
	
	int V_start, V_end, del_len, ins_len, i, j, k, L, n, RandIndex, N_err, N_del, N_ins;
	string germline, seq, err;
	
	// Build vectors for deletion/insertion length profiles
	vector <double> p_del, p_ins;
	p_del.clear();
	p_ins.clear();
	for(i=0;i<gap_bound;i++){
		p_del.push_back(params[del_params_begin+i]);
		p_ins.push_back(params[ins_params_begin+i]);
	}
	discrete_distribution<int> del_distr(p_del.begin(), p_del.end());
	discrete_distribution<int> ins_distr(p_ins.begin(), p_ins.end());
	
	for(n=0;n<N_seqs_synth;n++){
		
		N_err = (pt+n) -> N_err;
		N_del = (pt+n) -> N_del;
		N_ins = (pt+n) -> N_ins;
		
		germline = genomic_V[(pt+n) -> V_best_idx].second;
		V_start = (pt+n) -> V_start;
		V_end = (pt+n) -> V_end;
		
		(pt+n) -> error_pos = "[";
		(pt+n) -> error_list = "[";
		(pt+n) -> insertion_pos = "[";
		(pt+n) -> insertion_list = "[";
		(pt+n) -> deletion_pos = "[";
		(pt+n) -> deletion_list = "[";
		
		seq = germline;
		
		// insert mutations and hyper-indels on the germline
		vector <string> errors;
		for(i=0;i<N_err;i++){
			errors.push_back("err");
		}
		for(i=0;i<N_del;i++){
			errors.push_back("del");
		}
		for(i=0;i<N_ins;i++){
			errors.push_back("ins");
		}
		random_shuffle(errors.begin(),errors.end());
		
		// TODO: in what follows, index j has to be defined so to be coherent with the case V_start>0
		
		while((int)(errors.size())>0){
			RandIndex = rand() % (int)(errors.size());
			L = (int)(seq.size());
			if(errors.at(RandIndex)=="err"){
				j = rand() % L;
				if(seq.at(j)=='A'){
					err = notA[rand() % 3];
				}else if(seq.at(j)=='C'){
					err = notC[rand() % 3];
				}else if(seq.at(j)=='G'){
					err = notG[rand() % 3];
				}else if(seq.at(j)=='T'){
					err = notT[rand() % 3];
				}
				seq.replace(seq.begin()+j,seq.begin()+j+1,err);
				if((pt+n) -> error_pos.length()>1){
					(pt+n) -> error_pos.append(",");
					(pt+n) -> error_list.append(",");
				}
				(pt+n) -> error_pos.append(to_string(j));
				(pt+n) -> error_list.append(1,seq.at(j));
				(pt+n) -> error_list.append("->"+err);
			}else if(errors.at(RandIndex)=="del"){
				del_len = 1 + del_distr(generator);
				do{
					j = rand() % L;
				}while(j+del_len>L);
				err = "";
				seq.replace(seq.begin()+j,seq.begin()+j+del_len,err);
				if((pt+n) -> deletion_pos.length()>1){
					(pt+n) -> deletion_pos.append(",");
					(pt+n) -> deletion_list.append(",");
				}
				(pt+n) -> deletion_pos.append(to_string(j));
				(pt+n) -> deletion_list.append(to_string(del_len));
			}else if(errors.at(RandIndex)=="ins"){
				ins_len = 1 + ins_distr(generator);
				j = rand() % L;
				err = "";
				for(k=0;k<ins_len;k++){
					err += nuclAlphabet[rand() % 4];
				}
				seq.replace(seq.begin()+j,seq.begin()+j,err);
				if((pt+n) -> insertion_pos.length()>1){
					(pt+n) -> insertion_pos.append(",");
					(pt+n) -> insertion_list.append(",");
				}
				(pt+n) -> insertion_pos.append(to_string(j));
				(pt+n) -> insertion_list.append(to_string(ins_len)+"->");
				(pt+n) -> insertion_list.append(err);
			}
			errors.erase(errors.begin()+RandIndex);
		}
		
		(pt+n) -> seq = seq;
		
		(pt+n) -> error_pos.append("]");
		(pt+n) -> error_list.append("]");
		(pt+n) -> deletion_pos.append("]");
		(pt+n) -> deletion_list.append("]");
		(pt+n) -> insertion_pos.append("]");
		(pt+n) -> insertion_list.append("]");
		
	}
	
	return;
	
}

void write_generate_files(int N_seqs_synth, struct synth_seq *pt, string batch_name){
	
	int n;
	ofstream synth_seqs_fasta, synth_seqs_anchored, synth_scenarios_file;
	
	// Open files for generated sequences and their scenarios
	synth_seqs_fasta.open(batch_name + "generate/synthetic_seqs.fasta");
	synth_seqs_anchored.open(batch_name + "generate/synthetic_seqs_anchored.csv");
	synth_scenarios_file.open(batch_name + "generate/synthetic_scenarios.txt");
	
	synth_seqs_anchored << "seq_ID;aligned_seq_nt;V_best;V_best_start;V_best_end" << endl;
	synth_scenarios_file << "seq_ID;mu_err;V_choice;error_positions;error_list;deletion_pos;deletion_list;insertion_pos;insertion_list" << endl;
	
	for(n=0;n<N_seqs_synth;n++){
		
		synth_scenarios_file << n << ";";
		synth_scenarios_file << (pt+n) -> mu_err << ";";
		synth_scenarios_file << (pt+n) -> V_best_idx << ";";
		
		synth_seqs_fasta << ">" << n << "\n" << (pt+n) -> seq << "\n";
		synth_seqs_anchored << n << ";" << (pt+n) -> seq << ";" << (pt+n) -> V_best << ";" << (pt+n) -> V_start << ";" << (pt+n) -> V_end << endl;
		synth_scenarios_file << (pt+n) -> error_pos << ";" << (pt+n) -> error_list << ";";
		synth_scenarios_file << (pt+n) -> deletion_pos << ";" << (pt+n) -> deletion_list << ";";
		synth_scenarios_file << (pt+n) -> insertion_pos << ";" << (pt+n) -> insertion_list;
		synth_scenarios_file << endl;
	}
	
	synth_seqs_fasta.close();
	synth_seqs_anchored.close();
	synth_scenarios_file.close();
	
	return;
}
