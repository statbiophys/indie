/*
 *
 *  funcs.cpp
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

std::ostream& bold_on(std::ostream& os){
	return os << "\e[1m";
}

std::ostream& bold_off(std::ostream& os){
	return os << "\e[0m";
}

bool sort_using_smaller_than(double u, double v){
	return u < v;
}

bool sort_using_greater_than(double u, double v){
	return u > v;
}

vector <string> split_string(string line_str, char separator){
	stringstream str_str("");
	string sub_str = "";
	vector <string> vec_str;
	vec_str.clear();
	
	str_str << line_str;
	while(str_str.good()){
		getline(str_str,sub_str, separator);
		vec_str.push_back(sub_str);
	}
	return vec_str;
}

vector <pair<const int, const string>> read_fasta(string filename){
	//TODO: Check for \r,\n\s stuff
	ifstream infile(filename);
	if(!infile){
		throw runtime_error("File not found: "+filename);
	}
	string seq_str;
	string temp_str;
	int seq_count = -1;
	vector<pair<const int, const std::string>> sequence_vect;
	
	while (getline(infile,temp_str)){
		if(temp_str[temp_str.size()-1] == '\r'){
			temp_str.erase(temp_str.size()-1);
		}
		if(temp_str[0] == '>'){
			
			if(seq_count>(-1)){
				//Read sequences in upper case
				transform(seq_str.begin(),seq_str.end(),seq_str.begin(),::toupper);
				sequence_vect.push_back(pair<const int , const string> (seq_count , seq_str));
			}
			
			seq_str = string();
			seq_count++;
		}
		else{
			seq_str+=temp_str;
		}
		
	}
	if(seq_count>(-1)){
		//Read sequences in upper case
		transform(seq_str.begin(),seq_str.end(),seq_str.begin(),::toupper);
		sequence_vect.push_back(pair<const int , const string> (seq_count , seq_str));
	}
	
	return sequence_vect;
}

vector <pair<const int , const string>> read_indexed_csv(string filename){
	ifstream infile(filename);
	if(!infile){
		throw runtime_error("File not found: "+filename);
	}
	string line_str;
	vector <pair<const int, const std::string>> sequence_vect;
	getline(infile,line_str);
	while(getline(infile,line_str)){
		size_t semi_col_index = line_str.find(";");
		int index = stoi(line_str.substr(0,semi_col_index));
		string seq_str = line_str.substr(semi_col_index+1 , string::npos);
		transform(seq_str.begin() , seq_str.end() , seq_str.begin() , ::toupper);
		sequence_vect.push_back(pair<const int , const string >(index , seq_str));
	}
	return sequence_vect;
}

vector <double> educated_init_params(vector <double> p){
	int i;
	double sum;
	
	p.clear();
	
	p.push_back(0.05);    // repertoire-averaged error rate
	p.push_back(0.01);    // ratio between del_rate and error_rate
	p.push_back(0.01);    // ratio between ins_rate and error_rate
	
	sum = 0.;
	for(i=1;i<=gap_bound;i++){
		p.push_back(exp(-0.1*i));
		sum += exp(-0.1*i);
	}
	for(i=0;i<gap_bound;i++){
		p[del_params_begin+i] /= sum;
	}
	
	sum = 0.;
	for(i=1;i<=gap_bound;i++){
		p.push_back(exp(-0.1*i));
		sum += exp(-0.1*i);
	}
	for(i=0;i<gap_bound;i++){
		p[ins_params_begin+i] /= sum;
	}
	
	return p;
}

void update_gammas(){
	int n;
	
	for(n=0;n<gap_bound;n++){
		gamma_del_array[n] = params[1] * params[del_params_begin+n];
		log_gamma_del_array[n] = log(gamma_del_array[n]);
		gamma_ins_array[n] = params[2] * params[ins_params_begin+n];
		gamma_ins_array[n] *= pow(0.25,n+1);
		log_gamma_ins_array[n] = log(gamma_ins_array[n]);
		log_gamma_ins_array[n] += log(pow(0.25,n+1));
	}
	
	return;
}

vector <double> project_by_hand(vector <double> p){
	int i;
	double sum;
	
	// deletion profile
	sum = 0.;
	for(i=0;i<gap_bound;i++){
		sum += p[del_params_begin+i];
	}
	for(i=0;i<gap_bound;i++){
		p[del_params_begin+i] /= sum;
	}
	
	// insertion profile
	sum = 0.;
	for(i=0;i<gap_bound;i++){
		sum += p[ins_params_begin+i];
	}
	for(i=0;i<gap_bound;i++){
		p[ins_params_begin+i] /= sum;
	}
	
	return p;
}

vector <double> simplex_projection(vector <double> p){
	int i, rho;
	double theta, sum;
	vector <double> profile;
	
	// deletion profile
	profile.clear();
	for(i=0;i<gap_bound;i++){
		profile.push_back(p[del_params_begin+i]);
	}
	stable_sort(profile.begin(), profile.end(), sort_using_greater_than);
	rho = 0;
	sum = 0.;
	for(i=1;i<=gap_bound;i++){
		sum += profile[i-1];
		if( profile[i-1] > ((sum-1)/i) ){
			rho = i;
		}else{
			break;
		}
	}
	theta = (accumulate(profile.begin(), profile.begin()+rho, 0.0) - 1.0)/rho;
	for(i=0;i<gap_bound;i++){
		p[del_params_begin+i] = max(p[del_params_begin+i]-theta,0.);
	}
	
	// insertion profile
	profile.clear();
	for(i=0;i<gap_bound;i++){
		profile.push_back(p[ins_params_begin+i]);
	}
	stable_sort(profile.begin(), profile.end(), sort_using_greater_than);
	rho = 0;
	sum = 0.;
	for(i=1;i<=gap_bound;i++){
		sum += profile[i-1];
		if( profile[i-1] > ((sum-1)/i) ){
			rho = i;
		}else{
			break;
		}
	}
	theta = (accumulate(profile.begin(), profile.begin()+rho, 0.0) - 1.0)/rho;
	for(i=0;i<gap_bound;i++){
		p[ins_params_begin+i] = max(p[ins_params_begin+i]-theta,0.);
	}
	
	return p;
}

void check_model(vector <double> p){
	int params_size = (int)params.size();
	if(p[0] < 0. or p[0] > 1.){
		cout << endl << bold_on << "ERROR: " << bold_off;
		cout << "Average error rate does not belong to the [0,1] interval!" << endl << endl;
		exit(EXIT_FAILURE);
	}else if(p[1] < 0. or p[1] > 1.){
		cout << endl << bold_on << "ERROR: " << bold_off;
		cout << "Deletion rate does not belong to the [0,1] interval!" << endl << endl;
		exit(EXIT_FAILURE);
	}else if(p[2] < 0. or p[2] > 1.){
		cout << endl << bold_on << "ERROR: " << bold_off;
		cout << "Insertion rate does not belong to the [0,1] interval!" << endl << endl;
		exit(EXIT_FAILURE);
	}else if(p[0]*(1+p[1]+p[2]) > 1.){
		cout << endl << bold_on << "ERROR: " << bold_off;
		cout << "Sum of error, deletion and insertion rates is larger than 1!" << endl << endl;
		exit(EXIT_FAILURE);
	}else if(gap_bound != (params_size-del_params_begin)/2){
		cout << endl << bold_on << "ERROR: " << bold_off;
		cout << "Deletion and/or insertion profiles have a length not coherent with gap_bound value!" << endl << endl;
		exit(EXIT_FAILURE);
	}else if(abs(accumulate(p.begin()+del_params_begin, p.begin()+del_params_begin+gap_bound, 0.0) - 1.) > 1e-6){
		cout << endl << bold_on << "ERROR: " << bold_off;
		cout << "Deletion profile is not properly normalized!" << endl << endl;
		exit(EXIT_FAILURE);
	}else if(abs(accumulate(p.begin()+ins_params_begin, p.begin()+ins_params_begin+gap_bound, 0.0) - 1.) > 1e-6){
		cout << endl << bold_on << "ERROR: " << bold_off;
		cout << "Insertion profile is not properly normalized!" << endl << endl;
		exit(EXIT_FAILURE);
	}
	
}

void write_model_file(const string filename, double frac_mutated, vector <double> p){
	
	int i;
	ofstream out_file;
	
	out_file.open(filename);
	
	out_file << "@Rates" << endl;
	
	out_file << "#frac_mutated" << endl;
	out_file << frac_mutated << endl;
	
	out_file << "#aver_err_rate;del_ratio;ins_ratio" << endl;
	out_file << p[0] << ";" << p[1] << ";" << p[2] << endl;
	
	out_file << "@Profiles" << endl;
	
	out_file << "#max_gap_length" << endl;
	out_file << gap_bound << endl;
	
	out_file << "#del_profile" << endl;
	for(i=0;i<gap_bound;i++){
		out_file << p[del_params_begin+i];
		if(i<gap_bound-1){
			out_file << ";";
		}
	}
	out_file << endl;
	
	out_file << "#ins_profile" << endl;
	for(i=0;i<gap_bound;i++){
		out_file << p[ins_params_begin+i];
		if(i<gap_bound-1){
			out_file << ";";
		}
	}
	out_file << endl;
	
	out_file.close();
	
}

double gammapdf(double x, double alpha, double beta) {
	double value;
	//value = pow(beta, alpha) * pow(x,alpha-1) * pow(M_E,-beta*x) / tgamma(alpha);
	value = alpha*log(beta) + (alpha-1)*log(x) - beta*x - lgamma(alpha);
	if(value<-100){
		return 0.;
	}else{
		return exp(value);
	}
}

double gausspdf(double x, double mean, double sigma) {
	static const double inv_sqrt_2pi = 0.3989422804014327;
	double a = (x-mean)/sigma;
	return inv_sqrt_2pi/sigma * exp(-0.5*a*a);
}

double gausscdf(double x, double mean, double sigma){
	return 0.5 * erfc(-(x-mean)/(sigma*sqrt(2)));
}

double gausscdf_complement(double x, double mean, double sigma){
	return 0.5 * erfc(+(x-mean)/(sigma*sqrt(2)));
}
