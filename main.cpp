/*
 *
 *  main.cpp
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
#include "generate.h"
#include "greedy.h"

string last_version = "1.2.3";
string last_date = "14/12/2021";

// Recalling global variables from funcs.h
int gap_bound = max_gap_bound;
int del_params_begin = 3;
int ins_params_begin = del_params_begin + max_gap_bound;
int N_seqs;
vector <pair<const int, const string>> genomic_V;
vector <pair<const int, const string>> indexed_seqList;
vector <double> params;
double gamma_del_array[max_gap_bound], gamma_ins_array[max_gap_bound];
double log_gamma_del_array[max_gap_bound], log_gamma_ins_array[max_gap_bound];
double frac_mutated;
vector <double> prior, prior_mean_history, prior_std_history;
double prior_min = 0.;
double prior_max = 1.0;
double d_prior = 0.0005;
int N_prior = (int)((prior_max-prior_min+0.5*d_prior)/d_prior) + 1;
double *new_prior;
string smooth_method = "Gaussian_KDE";  // to be chosen in ['Gaussian_KDE', 'Gaussian_KDE_log', 'Gaussian_KDE_log_rescaled', 'Gaussian_KDE_soft_log']
double Gauss_KDE_sigma = 0.001;
string batch_name;
double frac_post_update = 0.25;  // fraction of posteriors to be updated at each inference step
int N_sampled_mu_errs = (int)(1./frac_post_update*5);  // number of sampled values for mu_err to be stored for each sequence
int W = 15;  // number of nodes to be used for the spline approximation of the posteriors (on the left for posteriors peaked in mu=0, on both sides for posteriors peaked in mu>0)

// Recalling global variables from full_lik.h
int max_T = 1000;
double eps = 1e-10;
double alpha = 1e-5;  // initial learning rate for error/deletion/insertion rates
double M = 100;  // constant to multiply the initial learning rate alpha for deletion/insertion profiles
double prune_thr = 1e-5;
double p_mutated;

// Random generator init
unsigned seed = chrono::system_clock::now().time_since_epoch().count();
minstd_rand0 generator (seed);

/* **************** */
/* ***** Main ***** */
/* **************** */

int main (int argc, char *argv[]) {
	
	int count, genomic_V_size, i, n, N_seqs_synth, N_seqs_unmutated, params_size, t, T;
	int old_precision, barWidth;
	double f0, f1, sum, decreaseRate, rm2var, prior_mean, prior_mean2, prior_std;
	double progress, omega;
	double del_sum, ins_sum;
	double frac_mutated_t0;
	double p;
	string filename, input_file, init_model_file, germline, str;
	string line_str;
	ifstream infile;
	ofstream greedy_file, aux_infer_file, aux_prior_file, aux_posterior_file;
	vector <double> sub_grad;
	vector <double> greedy_params, old_old_params, old_params, params_infer_t0, grad, pseudo_grad, alphas;
	vector <double> old_prior, prior_t0, prior_temp;
	vector <vector <double>> grad_history;
	vector <string> vec_str;
	struct synth_seq *synth_seqs;
	struct seq_align *all_seqs;
	
	#if AG
		double diff, tau;
	#endif
	#if GAMMA
		double mu;
	#endif
	
	// Specific variables when using greedy.cpp
	int Lx;
	double mu_err_greedy, mu_del_greedy, mu_ins_greedy, beta_del_greedy, beta_ins_greedy, greedy_frac_mutated;
	double P_del_greedy[max_gap_bound], P_ins_greedy[max_gap_bound];
	
	// Beginning of getopt section
	
	int c;
	
	// Default values for args
	batch_name = "";
	static int f_alpha = 0;
	static int f_generate = 0;
	static int f_inference = 0;
	static int f_init_model = 0;
	static int f_log = 0;
	static int f_T = 0;
	static int f_read = 0;
	
	while(1){
		int option_index = 0;
		
		static struct option long_options[] = {
			{"alpha",       required_argument,  NULL,  'a'},
			{"batch_name",  required_argument,  NULL,  'b'},
			{"generate",    required_argument,  NULL,  'g'},
			{"help",        no_argument,        NULL,  'h'},
			{"inference",   no_argument,        NULL,  'i'},
			{"log",         no_argument,        NULL,  'l'},
			{"parameters",  required_argument,  NULL,  'p'},
			{"read",        required_argument,  NULL,  'r'},
			{"max_iter",    required_argument,  NULL,  'T'},
			{"version",     no_argument,        NULL,  'v'},
			{NULL,          0,                  NULL,    0}
		};
		
		c = getopt_long(argc, argv, "-:a:b:g:hilp:r:T:v", long_options, &option_index);
		if(c==-1){
			break;
		}
		switch (c){
			case 'a':
				f_alpha = 1;
				// Hack to recognize argument starting with '-' as the next option
				// and not as the required argument for the current option
				if(optarg[0]!='-'){
					alpha = atof(optarg);
				}else{
					cout << endl << bold_on << "ERROR: " << bold_off;
					if(argv[optind-2][1]=='-'){
						cout << "Missing required argument for the option \'" << argv[optind-2] << "\'." << endl;
					}else{
						cout << "Missing required argument for the option \'" << argv[optind-2] << "\'." << endl;
					}
					cout << endl;
					exit(EXIT_FAILURE);
				}
				break;
			case 'b':
				// Hack to recognize argument starting with '-' as the next option
				// and not as the required argument for the current option
				if(optarg[0]!='-'){
					batch_name = optarg;
				}else{
					cout << endl << bold_on << "ERROR: " << bold_off;
					if(argv[optind-2][1]=='-'){
						cout << "Missing required argument for the option \'" << argv[optind-2] << "\'." << endl;
					}else{
						cout << "Missing required argument for the option \'" << argv[optind-2] << "\'." << endl;
					}
					cout << endl;
					exit(EXIT_FAILURE);
				}
				break;
			case 'g':
				f_generate = 1;
				// Hack to recognize argument starting with '-' as the next option
				// and not as the required argument for the current option
				if(optarg[0]!='-'){
					N_seqs_synth = atoi(optarg);
				}else{
					cout << endl << bold_on << "ERROR: " << bold_off;
					if(argv[optind-2][1]=='-'){
						cout << "Missing required argument for the option \'" << argv[optind-2] << "\'." << endl;
					}else{
						cout << "Missing required argument for the option \'" << argv[optind-2] << "\'." << endl;
					}
					cout << endl;
					exit(EXIT_FAILURE);
				}
				break;
			case 'h':
				//printf("List of command-line options\n");
				cout << endl << bold_on << "List of command-line options" << bold_off << endl << endl;
				printf("  -a, --alpha         [req. arg.] sets the global learning rate parameter to the specified value\n");
				printf("  -b, --batch_name    [req. arg.] creates a directory with the specified name\n");
				printf("                      and puts all the output files in it\n");
				printf("  -g, --generate      [req. arg.] generate N synthetic sequences (where N has to be inserted\n");
				printf("                      as argument\n");
				printf("  -h, --help          [no arg.] prints all the possible command-line options\n");
				printf("                      and if they require any (optional) argument\n");
				printf("  -i, --inference     [no arg.] performs the inference on model parameters\n");
				printf("                      by maximizing the total alignment likelihood\n");
				printf("  -l, --log           [no arg.] prints on auxiliary files optional info about the convergence\n");
				printf("                      of the algorithm during the inference step\n");
				printf("  -p, --parameters    [req. arg.] reads from the specified file the parameters to be used\n");
				printf("                      for the generation step or as the initial condition for the inference step\n");
				printf("  -r, --read          [req. arg.] reads from the specified file the dataset\n");
				printf("                      of sequences to analyze\n");
				printf("  -T, --max_iter      [req. arg.] sets the maximum number of gradient descent iterations\n");
				printf("\n");
				exit(EXIT_SUCCESS);
				break;
			case 'i':
				f_inference = 1;
				break;
			case 'l':
				f_log = 1;
				break;
			case 'p':
				f_init_model = 1;
				// Hack to recognize argument starting with '-' as the next option
				// and not as the required argument for the current option
				if(optarg[0]!='-'){
					init_model_file = optarg;
				}else{
					cout << endl << bold_on << "ERROR: " << bold_off;
					if(argv[optind-2][1]=='-'){
						cout << "Missing required argument for the option \'" << argv[optind-2] << "\'." << endl;
					}else{
						cout << "Missing required argument for the option \'" << argv[optind-2] << "\'." << endl;
					}
					cout << endl;
					exit(EXIT_FAILURE);
				}
				break;
			case 'r':
				f_read = 1;
				// Hack to recognize argument starting with '-' as the next option
				// and not as the required argument for the current option
				if(optarg[0]!='-'){
					input_file = optarg;
				}else{
					cout << endl << bold_on << "ERROR: " << bold_off;
					if(argv[optind-2][1]=='-'){
						cout << "Missing required argument for the option \'" << argv[optind-2] << "\'." << endl;
					}else{
						cout << "Missing required argument for the option \'" << argv[optind-2] << "\'." << endl;
					}
					cout << endl;
					exit(EXIT_FAILURE);
				}
				break;
			case 'T':
				f_T = 1;
				// Hack to recognize argument starting with '-' as the next option
				// and not as the required argument for the current option
				if(optarg[0]!='-'){
					T = atoi(optarg);
				}else{
					cout << endl << bold_on << "ERROR: " << bold_off;
					if(argv[optind-2][1]=='-'){
						cout << "Missing required argument for the option \'" << argv[optind-2] << "\'." << endl;
					}else{
						cout << "Missing required argument for the option \'" << argv[optind-2] << "\'." << endl;
					}
					cout << endl;
					exit(EXIT_FAILURE);
				}
				break;
			case 'v':
				cout << endl << "Current version " << last_version << " (last updated on " << last_date << ")." << endl << endl;
				exit(EXIT_SUCCESS);
				break;
			case '?':
				cout << endl << bold_on << "ERROR: " << bold_off;
				if(argv[optind-1][1]=='-'){
					cout << "The option \'" << argv[optind-1] << "\' is unknown." << endl;
				}else{
					cout << "The option \'-" << (char)optopt << "\' is unknown." << endl;
				}
				cout << endl;
				exit(EXIT_FAILURE);
				break;
			case ':':
				printf("Missing option for %c\n", optopt);
				cout << endl << bold_on << "ERROR: " << bold_off;
				if(argv[optind-1][1]=='-'){
					cout << "Missing required argument for the option \'" << argv[optind-1] << "\'." << endl;
				}else{
					cout << "Missing required argument for the option \'-" << (char)optopt << "\'." << endl;
				}
				cout << endl;
				exit(EXIT_FAILURE);
				break;
			default:
				cout << endl << bold_on << "ERROR: " << bold_off;
				cout << "The option corresponding to character code " << oct << c << " is unknown!" << endl;
				cout << endl;
				exit(EXIT_FAILURE);
				break;
		}
	}
	
	// Raise errors and exit
	if(f_generate==0 and f_read==0){
		cout << endl << bold_on << "ERROR: " << bold_off;
		cout << "Choose one option among \'-r\' [\'--read\'] or \'-g\' [\'--generate\']." << endl;
		cout << endl;
		exit(EXIT_FAILURE);
	}else if(f_generate==1 and f_read==1){
		cout << endl << bold_on << "ERROR: " << bold_off;
		cout << "Choose only one option among \'-r\' [\'--read\'] or \'-g\' [\'--generate\']." << endl;
		cout << endl;
		exit(EXIT_FAILURE);
	}
	if(f_read==1 and input_file.substr(input_file.find_last_of(".") + 1) != "csv") {
		cout << endl << bold_on << "ERROR: " << bold_off;
		cout << "Choose a valid csv file as input." << endl;
		cout << endl;
		exit(EXIT_FAILURE);
	}
	if(f_inference==0){
		if(f_alpha==1){
			cout << endl << bold_on << "ERROR: " << bold_off;
			cout << "\'-a\' [\'--alpha\'] option cannot be used without \'-f\' [\'--full_lik\'] option." << endl;
			cout << endl;
			exit(EXIT_FAILURE);
		}
		if(f_T==1){
			cout << endl << bold_on << "ERROR: " << bold_off;
			cout << "\'-T\' [\'--max_iter\'] option cannot be used without \'-f\' [\'--full_lik\'] option." << endl;
			cout << endl;
			exit(EXIT_FAILURE);
		}
	}
	
	// Print chosen settings
	cout << endl << bold_on << "Options" << bold_off << endl << endl;
	printf(" @ Input mode: ");
	if(f_read==0){
		printf("generate\n");
	}else{
		printf("read (from \"%s\")\n", input_file.c_str());
	}
	printf(" @ Batch name: \"%s\"\n", (batch_name!="") ? batch_name.c_str() : "None");
	if(f_inference==1){
		printf(" @ Alpha: ");
		if(f_alpha==0){
			printf("%g (default value)\n", alpha);
		}else{
			printf("%g (custom value)\n", alpha);
		}
		if(f_T==1){
			printf(" @ T: %d (custom value)\n", T);
		}else{
			printf(" @ T: %d (default value)\n", max_T);
			T = max_T;
		}
	}
	printf(" @ gap bound: %d\n", gap_bound);
	cout << endl;
	
	// End of getopt section
	
	if(batch_name!=""){
		system(&("mkdir -p " + batch_name)[0]);
		if(batch_name.back()!='/'){
			batch_name = batch_name + "/";
		}
		system(&("rm -rf " + batch_name + "*")[0]);
	}
	
	bool use_greedy = true;
	bool full_greedy = false;
	bool reduced_greedy = true;
	bool use_greedy_init = true;
	bool project_simplex = true;
	bool reduceMaxGap = false;
	bool use_momentum = true;
	
	clock_t start1, start2, end1, end2;
	
	/* ********** Loading the model ********** */
	
	cout << bold_on << "Model loading" << bold_off << endl << endl;
	
	// Import genomic templates
	filename = "./models/imgt_templates_IGHV_Homo_Sapiens_F.csv";
	
	// New method
	unordered_map <string, int> genomic_V_map = {};
	unordered_map <int, string> genomic_V_map_inverse = {};
	infile.open(filename);
	if(!infile){
		throw runtime_error("File not found: " + filename);
	}
	getline(infile,line_str);  // if we want to get rid of headers
	n = 0;
	while(getline(infile,line_str)){
		vec_str = split_string(line_str, ';');
		std::transform(vec_str[1].begin(), vec_str[1].end(), vec_str[1].begin(), ::toupper);
		genomic_V.push_back(pair <const int, const string> (n, vec_str[1]));
		genomic_V_map.insert(pair <string, int> (vec_str[0], n));
		genomic_V_map_inverse.insert(pair <int, string> (n, vec_str[0]));
		n++;
	}
	infile.close();
	
	genomic_V_size = (int)genomic_V.size();
	cout << "Genomic template(s) of " << genomic_V_size << " V gene(s) correctly imported!\n" << endl;
	
	if(f_init_model==0){
		filename = "./models/parms.txt";
	}else{
		filename = init_model_file;
	}
	
	infile.open(filename);
	if(!infile){
		throw runtime_error("File not found: " + filename);
	}
	
	params.clear();
	while(getline(infile,line_str)){
		if(line_str=="@Rates"){
			getline(infile,line_str);
			getline(infile,line_str);
			frac_mutated = stod(line_str);
			getline(infile,line_str);
			getline(infile,line_str);
			vec_str = split_string(line_str, ';');
			params.push_back(stod(vec_str[0]));
			params.push_back(stod(vec_str[1]));
			params.push_back(stod(vec_str[2]));
		}else if(line_str=="@Profiles"){
			getline(infile,line_str);
			getline(infile,line_str);
			gap_bound = stoi(line_str);
			getline(infile,line_str);
			getline(infile,line_str);
			vec_str = split_string(line_str, ';');
			for(i=0;i<(int)vec_str.size();i++){
				params.push_back(stod(vec_str[i]));
			}
			getline(infile,line_str);
			getline(infile,line_str);
			vec_str = split_string(line_str, ';');
			for(i=0;i<(int)vec_str.size();i++){
				params.push_back(stod(vec_str[i]));
			}
		}
	}
	infile.close();
	
	if(gap_bound > max_gap_bound){
		cout << endl << bold_on << "ERROR: " << bold_off;
		cout << "Input model has a max gap size larger than the maximum one allowed!" << endl;
		cout << endl;
		exit(EXIT_FAILURE);
	}
	
	ins_params_begin = del_params_begin + gap_bound;  // we need to re-init this value
	params_size = (int)params.size();
	check_model(params);
	
	params = simplex_projection(params);
	update_gammas();
	
	cout << "Model parameters:" << endl;
	cout << " - frac_mutated: " << frac_mutated << endl;
	cout << " - aver_err_rate: " << params[0] << endl;
	cout << " - del_ratio: " << params[1] << endl;
	cout << " - ins_ratio: " << params[2] << endl;
	cout << " - gap_bound: " << gap_bound << endl;
	cout << " - del_profile: [";
	for(i=0;i<gap_bound;i++){
		cout << params[del_params_begin+i];
		if(i<gap_bound-1){
			cout << ";";
		}else{
			cout << "]" << endl;
		}
	}
	cout << " - ins_profile: [";
	for(i=0;i<gap_bound;i++){
		cout << params[ins_params_begin+i];
		if(i<gap_bound-1){
			cout << ";";
		}else{
			cout << "]" << endl;
		}
	}
	cout << endl;
	
	// Write on a file the loaded model parameters
	write_model_file(batch_name + "loaded_params.txt", frac_mutated, params);
	
	/* ********** Generating a synthetic repertoire with hyper-indels ********** */
	
	if(f_generate){
		cout << bold_on << "Generation" << bold_off << endl << endl;
		
		system(&("mkdir -p " + batch_name + "generate")[0]);
		
		// Write on a file the model parameters used for the generation
		write_model_file(batch_name + "generate/generative_params.txt", frac_mutated, params);
		
		synth_seqs = new synth_seq[N_seqs_synth];
		
		init_synth_seqs(N_seqs_synth, synth_seqs, genomic_V_map_inverse, genomic_V_map);
		generate_synth_seqs_type_1(N_seqs_synth, synth_seqs);
		write_generate_files(N_seqs_synth, synth_seqs, batch_name);
		
		delete [] synth_seqs;
		
		infile.open(batch_name + "generate/synthetic_seqs_anchored.csv");
		if(!infile){
			throw runtime_error("File not found: " + filename);
		}
		getline(infile,line_str);  // if we want to get rid of headers
		n = 0;
		while(getline(infile,line_str)){
			n++;
		}
		infile.close();
		if(n != N_seqs_synth){
			cout << bold_on << "ERROR: " << bold_off;
			cout << "Number of generated sequences is not the one set as input!" << endl;
			cout << endl;
			exit(EXIT_FAILURE);
		}
		input_file = batch_name + "generate/synthetic_seqs_anchored.csv";
	}
	
	/* ********** Reading the file with sequences (either given as input or synthetically generated) ********** */
	
	if(f_read){
		cout << bold_on << "Input" << bold_off << endl << endl;
		system(&("mkdir -p " + batch_name + "read")[0]);
		system(&("cp " + input_file + " " + batch_name + "read/")[0]);
	}
	
	// New method
	infile.open(input_file);
	if(!infile){
		throw runtime_error("File not found: " + input_file);
	}
	getline(infile,line_str);  // if we want to get rid of headers
	N_seqs = 0;
	while(getline(infile,line_str)){
		N_seqs++;
	}
	infile.close();
	
	// Allocate the memory where to store the alignment info
	all_seqs = new seq_align[N_seqs];
	infile.open(input_file);
	getline(infile,line_str);  // if we want to get rid of headers
	n = 0;
	while(getline(infile,line_str)){
		vec_str = split_string(line_str, ';');
		std::transform(vec_str[1].begin(), vec_str[1].end(), vec_str[1].begin(), ::toupper);
		std::transform(vec_str[2].begin(), vec_str[2].end(), vec_str[2].begin(), ::toupper);
		all_seqs[n].seq_ID = n;
		all_seqs[n].seq = vec_str[1];
		all_seqs[n].V_best = vec_str[2];
		all_seqs[n].V_best_idx = genomic_V_map[all_seqs[n].V_best];
		all_seqs[n].V_start = stoi(vec_str[3]);
		all_seqs[n].V_end = stoi(vec_str[4]);
		n++;
	}
	infile.close();
	
	if(f_generate){
		cout << "Synthetic repertoire of " << N_seqs << " sequence(s) correctly generated!\n" << endl;
	}else if(f_read){
		cout << "Repertoire of " << N_seqs << " sequence(s) correctly imported!\n" << endl;
	}
	
	/* ********** Greedy ********** */
	
	if(use_greedy==true){
		// Greedy algorithm to evaluate the best initial condition for the full-likelihood maximization
		
		system(&("mkdir -p " + batch_name + "greedy")[0]);
		
		frac_mutated = 0.8;
		params = educated_init_params(params);
		update_gammas();
		write_model_file(batch_name + "greedy/initial_params.txt", frac_mutated, params);
		
		// Parallel computing
		start1 = clock();
		start2 = omp_get_wtime();
		if(full_greedy==true){
			// Here we use the full greedy algorithm (runs over all the V germlines and finds the most likely one)
			cout << bold_on << "Full greedy alignment" << bold_off << endl << endl;
			#pragma omp parallel for
			for(n=0;n<N_seqs;n++){
				all_seqs[n] = greedy_likelihood(indexed_seqList[n].second);
			}
		}else if(reduced_greedy==true){
			// Here we use the reduced greedy algorithm (runs only over the V germline given in the input file)
			cout << bold_on << "Reduced greedy alignment" << bold_off << endl << endl;
			#pragma omp parallel for
			for(n=0;n<N_seqs;n++){
				all_seqs[n] = reduced_greedy_likelihood(all_seqs+n);
			}
		}
		end1 = clock();
		end2 = omp_get_wtime();
		cout << "Execution time: " << double(end2-start2) << " s" << endl;
		cout << "Avg. time per seq: " << double(end1-start1) / CLOCKS_PER_SEC / N_seqs << " s" << endl;
		cout << endl;
		
		// Compute greedy parameters
		greedy_frac_mutated = 0.;
		mu_err_greedy = 0.;
		mu_del_greedy = 0.;
		mu_ins_greedy = 0.;
		for(i=0;i<gap_bound;i++){
			P_del_greedy[i] = 0.;
			P_ins_greedy[i] = 0.;
		}
		for(n=0;n<N_seqs;n++){
			Lx = all_seqs[n].V_end - all_seqs[n].V_start;
			if(all_seqs[n].N_err>0 or all_seqs[n].N_del>0 or all_seqs[n].N_ins>0){
				greedy_frac_mutated += 1.0;
				all_seqs[n].mutated = true;
			}else{
				all_seqs[n].mutated = false;
			}
			mu_err_greedy += (double)(all_seqs[n].N_err)/Lx;
			mu_del_greedy += (double)(all_seqs[n].N_del)/Lx;
			mu_ins_greedy += (double)(all_seqs[n].N_ins)/Lx;
			if(all_seqs[n].vec_del.size()>0){
				for(i=0;i<(int)all_seqs[n].vec_del.size();i++){
					P_del_greedy[all_seqs[n].vec_del[i]-1] += 1.;
				}
			}
			if(all_seqs[n].vec_ins.size()>0){
				for(i=0;i<(int)all_seqs[n].vec_ins.size();i++){
					P_ins_greedy[all_seqs[n].vec_ins[i]-1] += 1.;
				}
			}
		}
		greedy_frac_mutated /= N_seqs;
		mu_err_greedy /= N_seqs;
		mu_del_greedy /= N_seqs;
		mu_ins_greedy /= N_seqs;
		beta_del_greedy = mu_del_greedy/mu_err_greedy;
		beta_ins_greedy = mu_ins_greedy/mu_err_greedy;
		
		for(n=0;n<N_seqs;n++){
			Lx = all_seqs[n].V_end - all_seqs[n].V_start;
			all_seqs[n].mu_err = (double)(max(1,all_seqs[n].N_err+all_seqs[n].N_del+all_seqs[n].N_ins))/Lx/(1.+beta_del_greedy+beta_ins_greedy);
		}
		
		// Write on file the alignment for each sequence
		filename = batch_name + "greedy/greedy_alignment_summary.csv";
		greedy_file.open(filename);
		greedy_file << "seq_ID;seq_len;N_err;N_del;N_ins;del_len_list;ins_len_list" << endl;
		for(n=0;n<N_seqs;n++){
			Lx = all_seqs[n].V_end - all_seqs[n].V_start;
			greedy_file << n << ";" << Lx << ";" << all_seqs[n].N_err << ";" << all_seqs[n].N_del << ";" << all_seqs[n].N_ins << ";" << "[";
			for(i=0;i<(int)all_seqs[n].vec_del.size();i++){
				greedy_file << all_seqs[n].vec_del[i];
				if(i<(int)all_seqs[n].vec_del.size()-1){
					greedy_file << ",";
				}
			}
			greedy_file << "];[";
			for(i=0;i<(int)all_seqs[n].vec_ins.size();i++){
				greedy_file << all_seqs[n].vec_ins[i];
				if(i<(int)all_seqs[n].vec_ins.size()-1){
					greedy_file << ",";
				}
			}
			greedy_file << "]" << endl;
		}
		greedy_file.close();
		
		cout << "Model parameters from greedy algorithm:" << endl;
		cout << " - frac_mutated: " << greedy_frac_mutated << endl;
		cout << " - aver_err_rate: " << mu_err_greedy << endl;
		cout << " - del_ratio: " << beta_del_greedy << endl;
		cout << " - ins_ratio: " << beta_ins_greedy << endl;
		cout << " - gap_bound: " << gap_bound << endl;
		cout << " - del_counts: [";
		del_sum = 0.;
		for(i=0;i<gap_bound;i++){
			cout << P_del_greedy[i];
			if(i<gap_bound-1){
				cout << ";";
			}else{
				cout << "]" << endl;
			}
			del_sum += P_del_greedy[i];
		}
		cout << " - ins_counts: [";
		ins_sum = 0.;
		for(i=0;i<gap_bound;i++){
			cout << P_ins_greedy[i];
			if(i<gap_bound-1){
				cout << ";";
			}else{
				cout << "]" << endl;
			}
			ins_sum += P_ins_greedy[i];
		}
		cout << " - del_profile: [";
		for(i=0;i<gap_bound;i++){
			P_del_greedy[i] /= del_sum;
			cout << P_del_greedy[i];
			if(i<gap_bound-1){
				cout << ";";
			}else{
				cout << "]" << endl;
			}
		}
		cout << " - ins_profile: [";
		for(i=0;i<gap_bound;i++){
			P_ins_greedy[i] /= ins_sum;
			cout << P_ins_greedy[i];
			if(i<gap_bound-1){
				cout << ";";
			}else{
				cout << "]" << endl;
			}
		}
		
		cout << endl;
		
		greedy_params.clear();
		greedy_params.push_back(mu_err_greedy);
		greedy_params.push_back(mu_del_greedy/mu_err_greedy);
		greedy_params.push_back(mu_ins_greedy/mu_err_greedy);
		for(i=0;i<gap_bound;i++){
			greedy_params.push_back(P_del_greedy[i]);
		}
		for(i=0;i<gap_bound;i++){
			greedy_params.push_back(P_ins_greedy[i]);
		}
		
		// Write on a file the final model parameters got through the greedy algorithm
		write_model_file(batch_name + "greedy/" + "greedy_params.txt", greedy_frac_mutated, greedy_params);
		
		// We reduce the largest allowed size for gaps according to the greedy alignment
		if(reduceMaxGap==true){
			i = gap_bound - 1;
			while(P_del_greedy[i]==0. and P_ins_greedy[i]==0.){
				i -= 1;
			}
			gap_bound = min(gap_bound, i + 3);  // let's have some more space, just to be sure
			ins_params_begin = del_params_begin + gap_bound;  // we need to reset this value
			printf("New value for gap_bound: %d\n", gap_bound);
		}
		
	}
	
	/* ********** Full-likelihood inference ********** */
	
	// Inference of hyper-indel parameters through full-likelihood maximization
	if(f_inference==1){
		
		cout << bold_on << "Full likelihood maximization" << bold_off << endl << endl;
		
		system(&("mkdir -p " + batch_name + "inference/")[0]);
		
		if(use_greedy_init==true){
			frac_mutated = greedy_frac_mutated;
			if(params[0]>0. and params[1]>0. and params[2]>0.){
				// Set as initial condition the values got from greedy alignment
				cout << "Initial condition for the inference taken from greedy alignment." << endl << endl;
				params = greedy_params;
				// Add pseudo-counts to regularize length profiles for the inference
				for(i=0;i<gap_bound;i++){
					params[i+del_params_begin] += 1./N_seqs;
					params[i+ins_params_begin] += 1./N_seqs;
				}
				params = simplex_projection(params);
				update_gammas();
			}else{
				// Greedy alignment has detected no errors and/or no indels
				cout << "Greedy alignment assigned zero-valued rates." << endl;
				cout << "Rates are safely set to a minimal value; length profiles are taken from initially loaded model." << endl << endl;
				params[0] = max(0.001,params[0]);
				params[1] = max(0.001,params[1]);
				params[2] = max(0.001,params[2]);
			}
			if(project_simplex==false){
				params = project_by_hand(params);
			}else{
				params = simplex_projection(params);
			}
			update_gammas();
			// Set the initial prior as uniform
			prior.clear();
			new_prior = new double[N_prior];
			for(i=0;i<N_prior;i++){
				prior.push_back(1./N_prior);
				new_prior[i] = prior.at(i);
			}
			old_prior.clear();
			old_prior = prior;
		}
		
		// Prune the alignment matrix before starting the iterative procedure
		start1 = clock();
		start2 = omp_get_wtime();
		cout << "Pruning alignment matrices...";
		cout.flush();
		#pragma omp parallel for
		for(n=0;n<N_seqs;n++){
			prune_by_new_method(all_seqs+n);
		}
		end1 = clock();
		end2 = omp_get_wtime();
		cout << " done." << endl;
		cout << "Execution time: " << double(end2-start2) << " s" << endl;
		cout << "Avg. time per seq: " << double(end1-start1) / CLOCKS_PER_SEC / N_seqs << " s" << endl << endl;
		
		N_seqs_unmutated = 0;
		frac_mutated = 0.;
		for(n=0;n<N_seqs;n++){
			if(all_seqs[n].mutated == false){
				all_seqs[n].p_mutated = 0.25;
				N_seqs_unmutated++;
			}else{
				all_seqs[n].p_mutated = 0.9;
			}
			frac_mutated += all_seqs[n].p_mutated;
		}
		frac_mutated /= N_seqs;
		cout << "Fraction of (apparently) unmutated sequences: " << (double)N_seqs_unmutated/N_seqs << endl << endl;
		
		// Write on a file the initial model parameters for the inference step
		write_model_file(batch_name + "inference/" + "initial_params.txt", frac_mutated, params);
		
		// If we want to list all the likely scenarios
		if(false){
			explore_scenarios(all_seqs); // or all_seqs+n for the (n+1)-th sequence
		}
		
		// Store the initial condition
		frac_mutated_t0 = frac_mutated;
		params_infer_t0 = params;
		prior_t0 = prior;
		
		bool whileFlag = false;    // remains false until the inference is not properly accomplished (i.e. until inference_status remains false)
		int whileCount = 0;
		bool nanFlag;
		
		while(!whileFlag){
			
			if(whileCount>0){
				// reduce the global learning rate
				alpha *= 0.75;
				cout << endl << "Learning rate alpha too large. Now reduced to: " << alpha << endl;
			}else{
				cout << "Initial value for the learning rate alpha: " << alpha << endl;
			}
			
			whileCount++;
			nanFlag = false;    // remains false unless the likelihood becomes nan, exiting with error from the inference
			
			// Set/reset the initial condition
			frac_mutated = frac_mutated_t0;
			p_mutated = frac_mutated;
			params = params_infer_t0;
			params = simplex_projection(params);
			update_gammas();
			old_old_params = params;
			old_params = params;
			alphas.clear();
			alphas.resize(del_params_begin,alpha);
			alphas.resize(del_params_begin+2*gap_bound,alpha*M);
			prior = prior_t0;
			old_prior = prior;
			prior_mean_history.clear();
			prior_std_history.clear();
			for(n=0;n<N_seqs;n++){
				Lx = all_seqs[n].V_end - all_seqs[n].V_start;
				all_seqs[n].mu_err = (double)(max(1,all_seqs[n].N_err+all_seqs[n].N_del+all_seqs[n].N_ins))/Lx/(1.+beta_del_greedy+beta_ins_greedy);
			}
			#pragma omp parallel for
			for(n=0;n<N_seqs;n++){
				#if GAMMA
					init_gamma_posterior(all_seqs+n);
				#endif
				#if SPLINE
					init_spline_posterior(all_seqs+n);
				#endif
			}
			// Allocate memory for auxiliary vectors
			grad.clear();
			grad.resize(params_size,0.);
			pseudo_grad.clear();
			pseudo_grad.resize(params_size,0.);
			grad_history.clear();
			for(i=0;i<params_size;i++){
				sub_grad.clear();
				sub_grad.resize(5, 0.);
				grad_history.push_back(sub_grad);
			}
			
			if(f_log==1){
				// Open output file for parameters
				aux_infer_file.open(batch_name + "inference/aux_infer_file_attempt_" + to_string(whileCount) + ".txt");
				aux_prior_file.open(batch_name + "inference/aux_prior_file_attempt_" + to_string(whileCount) + ".txt");
				aux_posterior_file.open(batch_name + "inference/aux_posterior_file_attempt_" + to_string(whileCount) + ".txt");
				// Write the headers
				aux_infer_file << "t;negLogLik;rm2var;frac_mutated;aver_err_rate;del_ratio;ins_ratio;";
				for(i=1;i<=gap_bound;i++){
					aux_infer_file << "P_del_" << i << ";";
				}
				for(i=1;i<=gap_bound;i++){
					aux_infer_file << "P_ins_" << i << ";";
				}
				aux_infer_file << "grad_aver_err_rate;grad_del_ratio;grad_ins_ratio;";
				for(i=1;i<=gap_bound;i++){
					aux_infer_file << "grad_P_del_" << i << ";";
				}
				for(i=1;i<=gap_bound;i++){
					aux_infer_file << "grad_P_ins_" << i;
					if(i<gap_bound)
					aux_infer_file << ";";
					else
					aux_infer_file << endl;
				}
				aux_prior_file << "t;";
				for(i=0;i<N_prior;i++){
					aux_prior_file << "Prior_" << i;
					if(i<N_prior-1){
						aux_prior_file << ";";
					}
				}
				aux_prior_file << endl;
				aux_prior_file << "mu;";
				for(i=0;i<N_prior;i++){
					aux_prior_file << prior_min + i*d_prior;
					if(i<N_prior-1){
						aux_prior_file << ";";
					}
				}
				aux_prior_file << endl;
				#if GAMMA
					aux_posterior_file << "t;p_mutated;alphas;betas;drawn_mu" << endl;
				#endif
				#if SPLINE
					aux_posterior_file << "t;p_mutated;mode;mean;sigma;drawn_mu" << endl;
				#endif
			}
			
			progress = 0.;
			barWidth = 50;
			rm2var = 0.;
			
			for(t=1;(t<=T and !nanFlag);t++){
				
				cout << "Inference progress:" << "\r";
				cout.flush();
				
				start1 = clock();
				start2 = omp_get_wtime();
				
				// print all posteriors at a certain time step
				/*
				if(t==50){
					double old_mu_err_0 = all_seqs[0].mu_err;
					double old_mu_err_1 = all_seqs[1].mu_err;
					ofstream tmp_file;
					for(n=0;n<N_seqs;n++){
						if(n==0 or n==1){
							tmp_file.open(batch_name + "example_of_posterior_n" + to_string(n) + "_t" + to_string(t) + ".csv");
							for(i=0;i<N_prior;i+=5){
								tmp_file << prior_min + i*d_prior;
								if(i<N_prior-5){
									tmp_file << ";";
								}
							}
							tmp_file << endl;
							for(i=0;i<N_prior;i+=5){
								all_seqs[n].mu_err = prior_min + i*d_prior;
								tmp_file << log(partial_likelihood_new_method(all_seqs+n));
								if(i<N_prior-5){
									tmp_file << ";";
								}
							}
							tmp_file << endl;
							for(i=0;i<N_prior;i+=5){
								tmp_file << log(prior.at(i));
								if(i<N_prior-5){
									tmp_file << ";";
								}
							}
							tmp_file << endl;
							tmp_file.close();
						}
					}
					all_seqs[0].mu_err = old_mu_err_0;
					all_seqs[1].mu_err = old_mu_err_1;
				}
				*/
				
				// Schedule for the decreasing of the learning rate
				decreaseRate = exp(-2.*t/T);
				
				// Update parameters to maximize likelihood
				if(use_momentum==true){
					omega = 0.9 * (double)t/(t+3);
					for(i=1;i<params_size;i++){
						params[i] = old_params[i] + omega * (old_params[i] - old_old_params[i]);
					}
				}
				
				params = simplex_projection(params);
				update_gammas();
				
				// I can sample a point-mutation rate for each sequence from its posterior
				for(n=0;n<N_seqs;n++){
					if(frand > all_seqs[n].p_mutated){
						all_seqs[n].mu_err = 0.;
					}else{
						#if GAMMA
							gamma_distribution<double> distribution(all_seqs[n].alpha,1./all_seqs[n].beta);
							do{
								all_seqs[n].mu_err = distribution(generator);
							}while(all_seqs[n].mu_err<=prior_min or all_seqs[n].mu_err>=prior_max);
						#endif
						#if SPLINE
							all_seqs[n].mu_err = all_seqs[n].sampled_mu_errs[rand() % (int)all_seqs[n].sampled_mu_errs.size()];
						#endif
					}
				}
				
				// Here I refresh the pruning for the alignment matrices
				/*
				if(t>0 and t%50==0){
					#pragma omp parallel for
					for(n=0;n<N_seqs;n++){
						prune_by_new_method(all_seqs+n);
					}
				}
				*/
				
				f0 = negLogLikelihood(all_seqs,"pruned");
				if(!isfinite(f0)){
					nanFlag = true;
				}
				
				grad[0] = 0.; // average mutation rate is changed with a different method
				for(i=1;i<params_size;i++){
					params[i] += eps;
					update_gammas();
					f1 = negLogLikelihood(all_seqs,"pruned");
					params[i] -= eps;
					update_gammas();
					grad[i] = (f1-f0)/eps;
					grad[i] /= (N_seqs*frac_mutated);
				}
				
				// Project by hand the gradient orthogonally
				// to the normalization-preserving manifold
				if(true){
					del_sum = 0.;
					count = 0;
					for(i=0;i<gap_bound;i++){
						del_sum += grad[del_params_begin+i];
						count ++;
					}
					del_sum /= count;
					for(i=0;i<gap_bound;i++){
						grad[del_params_begin+i] -= del_sum;
					}
					ins_sum = 0.;
					count = 0;
					for(i=0;i<gap_bound;i++){
						ins_sum += grad[ins_params_begin+i];
						count ++;
					}
					ins_sum /= count;
					for(i=0;i<gap_bound;i++){
						grad[ins_params_begin+i] -= ins_sum;
					}
				}
				
				// Write the parameters on the file
				//   The time-index convention is such that at iteration t, index time is t-1,
				//   parameters are the one obtained at the end of iteration t-1 (with the corresponding total likelihood)
				//   and gradients are the ones computed at iteration t, starting from the parameters above
				if(f_log==1){
					aux_infer_file << t-1 << ";";
					aux_infer_file << f0 << ";";
					aux_infer_file << scientific << rm2var << ";" << defaultfloat;
					aux_infer_file << frac_mutated << ";";
					for(i=0;i<params_size;i++){
						aux_infer_file << old_params[i] << ";";
					}
					for(i=0;i<params_size;i++){
						aux_infer_file << grad[i];
						if(i<params_size-1){
							aux_infer_file << ";";
						}else{
							aux_infer_file << endl;
						}
					}
					
					aux_prior_file << t-1 << ";";
					for(i=0;i<N_prior;i++){
						aux_prior_file << old_prior[i];
						if(i<N_prior-1){
							aux_prior_file << ";";
						}else{
							aux_prior_file << endl;
						}
					}
					
					aux_posterior_file << t-1 << ";" << "[";
					for(n=0;n<N_seqs;n++){
						aux_posterior_file << all_seqs[n].p_mutated;
						if(n<N_seqs-1){
							aux_posterior_file << ",";
						}else{
							aux_posterior_file << "];[";
						}
					}
					#if GAMMA
						for(n=0;n<N_seqs;n++){
							aux_posterior_file << all_seqs[n].alpha;
							if(n<N_seqs-1){
								aux_posterior_file << ",";
							}else{
								aux_posterior_file << "];[";
							}
						}
						for(n=0;n<N_seqs;n++){
							aux_posterior_file << all_seqs[n].beta;
							if(n<N_seqs-1){
								aux_posterior_file << ",";
							}else{
								aux_posterior_file << "];[";
							}
						}
					#endif
					#if SPLINE
						for(n=0;n<N_seqs;n++){
							aux_posterior_file << all_seqs[n].mode;
							if(n<N_seqs-1){
								aux_posterior_file << ",";
							}else{
								aux_posterior_file << "];[";
							}
						}
						for(n=0;n<N_seqs;n++){
							aux_posterior_file << all_seqs[n].mean;
							if(n<N_seqs-1){
								aux_posterior_file << ",";
							}else{
								aux_posterior_file << "];[";
							}
						}
						for(n=0;n<N_seqs;n++){
							aux_posterior_file << all_seqs[n].sigma;
							if(n<N_seqs-1){
								aux_posterior_file << ",";
							}else{
								aux_posterior_file << "];[";
							}
						}
					#endif
					for(n=0;n<N_seqs;n++){
						aux_posterior_file << all_seqs[n].mu_err;
						if(n<N_seqs-1){
							aux_posterior_file << ",";
						}else{
							aux_posterior_file << "]" << endl;
						}
					}
					
				}
				
				// If we want to avoid too big jumps
				for(i=1;i<params_size;i++){
					if(i>=0 and i<del_params_begin){
						pseudo_grad[i] = alphas[i] * grad[i] * decreaseRate;
						while(abs(pseudo_grad[i])>0.2*params[i]){
							pseudo_grad[i] *= 0.5;
						}
					}else{
						pseudo_grad[i] = alphas[i] * grad[i] * decreaseRate;
					}
				}
				
				for(i=1;i<params_size;i++){
					params[i] = params[i] - pseudo_grad[i];
				}
				
				if(project_simplex==true){
					params = simplex_projection(params);
				}else{
					params = project_by_hand(params);
				}
				update_gammas();
				// Check boundaries and normalizations
				del_sum = 0.;
				ins_sum = 0.;
				for(i=0;i<params_size;i++){
					if(params[i] < 0){
						cout << "params[" << i << "]: " << params[i] << endl;
						params[i] = eps;
					}else if(params[i] > 1){
						cout << "params[" << i << "]: " << params[i] << endl;
						params[i] = 1 - eps;
					}
					if(i>=del_params_begin and i<ins_params_begin){
						del_sum += params[i];
					}else if(i>=ins_params_begin and i<ins_params_begin+gap_bound){
						ins_sum += params[i];
					}
				}
				if(abs(del_sum-1.)>eps or abs(ins_sum-1.)>eps){
					cout << "Deletion and insertion profiles are not being properly normalized!" << endl << endl;
					exit(EXIT_FAILURE);
				}
				
				rm2var = 0.;
				for(i=0;i<params_size;i++){
					rm2var += (old_params[i]-params[i]) * (old_params[i]-params[i]);
					old_old_params[i] = old_params[i];
					old_params[i] = params[i];
				}
				rm2var = sqrt(rm2var/params_size);
				update_gammas();
				
				// Before updating posteriors and prior, check the new parameters
				f0 = negLogLikelihood(all_seqs,"pruned");
				if(!isfinite(f0)){
					nanFlag = true;
				}else{
					
					for(i=0;i<N_prior;i++){
						new_prior[i] = 0.;
					}
					
					// Here I update the posterior for each sequence
					sum = 0.;
					#pragma omp parallel for reduction (+:sum) shared(new_prior) private(p)
					for(n=0;n<N_seqs;n++){
						#if GAMMA
							if(all_seqs[n].mutated == true or all_seqs[n].mutated == false){
								sum += update_gamma_posterior(all_seqs+n, t);
							}
						#endif
						#if SPLINE
							p = frand;
							if(p < frac_post_update){
								update_spline_posterior(all_seqs+n, t);
							}
						#endif
					}
					if(!isfinite(sum)){
						nanFlag = true;
					}else{
						// Here I update the prior at the repertoire-level
						p_mutated = 0.;
						for(n=0;n<N_seqs;n++){
							p_mutated += all_seqs[n].p_mutated;
						}
						p_mutated /= N_seqs;
						frac_mutated = p_mutated;
						#if GAMMA
							for(i=0;i<N_prior;i++){
								mu = prior_min + i*d_prior;
								sum = 0.;
								#pragma omp parallel for reduction (+:sum)
								for(n=0;n<N_seqs;n++){
									if(all_seqs[n].alpha > 1.){
										sum += all_seqs[n].p_mutated * gammapdf(mu, all_seqs[n].alpha, all_seqs[n].beta);
									}
									else{
										sum += all_seqs[n].p_mutated * all_seqs[n].beta * exp(-all_seqs[n].beta*mu);
									}
								}
								prior[i] = sum / (N_seqs*frac_mutated);
							}
						#endif
						#if SPLINE
							sum = 0.;
							for(i=0;i<N_prior;i++){
								sum += new_prior[i];
							}
							sum -= 0.5*new_prior[0];
							sum *= d_prior;
							for(i=0;i<N_prior;i++){
								new_prior[i] /= sum;
								prior.at(i) = 0.5*(1.0-frac_post_update)*prior.at(i) + 0.5*frac_post_update*new_prior[i];
							}
						#endif
						// Here I smooth the prior obtained above, and then check its mean and std
						if(t==1){
							prior_mean = 0.;
							prior_mean2 = 0.;
							prior_std = 0.;
							for(i=0;i<N_prior;i++){
								prior_mean += (prior_min + i*d_prior) * prior.at(i);
								prior_mean2 += (prior_min + i*d_prior) * (prior_min + i*d_prior) * prior.at(i);
							}
							prior_mean *= d_prior;
							prior_mean2 *= d_prior;
							prior_std = sqrt(prior_mean2-prior_mean*prior_mean);
							if(smooth_method=="Gaussian_KDE"){
								Gauss_KDE_sigma = Gauss_KDE_sigma;
							}else if(smooth_method=="Gaussian_KDE_log"){
								Gauss_KDE_sigma = prior_std;
							}else if(smooth_method=="Gaussian_KDE_log_rescaled"){
								Gauss_KDE_sigma = prior_std;
							}else if(smooth_method=="Gaussian_KDE_soft_log"){
								Gauss_KDE_sigma = prior_std;
							}
							cout << "Gauss_KDE_sigma initially set to: " << Gauss_KDE_sigma  << " vs prior_std=" << prior_std << endl;
						}
						prior = smooth_prior(prior, prior_mean, smooth_method);
						old_prior = prior;
						prior_mean = 0.;
						prior_mean2 = 0.;
						prior_std = 0.;
						for(i=0;i<N_prior;i++){
							prior_mean += (prior_min + i*d_prior) * prior.at(i);
							prior_mean2 += (prior_min + i*d_prior) * (prior_min + i*d_prior) * prior.at(i);
						}
						prior_mean *= d_prior;
						prior_mean2 *= d_prior;
						prior_std = sqrt(prior_mean2-prior_mean*prior_mean);
						params[0] = prior_mean;
						prior_mean_history.push_back(prior_mean);
						prior_std_history.push_back(prior_std);
					}
					
				}
				
				end1 = clock();
				end2 = omp_get_wtime();
				
				// Print the progress bar
				cout << "Inference progress: [";
				int pos = barWidth * progress;
				for (i=0;i<barWidth;i++) {
					if(i <= pos){
						cout << "#";
					}else{
						cout << "_";
					}
				}
				progress += 1./T;
				old_precision = cout.precision();
				cout.precision(old_precision);
				cout << "] " << fixed << setprecision(2) << (progress * 100.0) << " %" << fixed << setprecision(3) << "    (Total iteration time: " << int(end2-start2) << " s; Avg. time per seq: " << double(end1-start1) / CLOCKS_PER_SEC / N_seqs << " s; Root mean square variation: " << scientific << rm2var << defaultfloat << ")          " << "\r";
				cout.precision(old_precision);
				cout.flush();
				
			} // closes the for loop with index t
			
			// If for loop is properly finished
			if(nanFlag==false){
				// Update the flag to exit the while loop
				whileFlag = true;
				// Last update of the progress bar
				cout << "Inference progress: [";
				for (i=0;i<barWidth;i++) {
					cout << "#";
				}
				cout << "] " << "Done.    ";
				cout.flush();
			}else{
				cout << endl;
			}
			
			if(f_log==1){
				aux_infer_file.close();
				aux_prior_file.close();
				aux_posterior_file.close();
			}
			
		} // closes the while loop
		
		// Write on a file the inferred model parameters
		write_model_file(batch_name + "inference/" + "inferred_params.txt", frac_mutated, params);
		
		
		cout << endl << endl;
		
		cout << "Model parameters from inference:" << endl;
		cout << " - frac_mutated: " << frac_mutated << endl;
		cout << " - aver_err_rate: " << params.at(0) << endl;
		cout << " - del_ratio: " << params.at(1) << endl;
		cout << " - ins_ratio: " << params.at(2) << endl;
		cout << " - gap_bound: " << gap_bound << endl;
		cout << " - del_profile: [";
		for(i=0;i<gap_bound;i++){
			cout << params.at(del_params_begin+i);
			if(i<gap_bound-1){
				cout << ";";
			}else{
				cout << "]" << endl;
			}
		}
		cout << " - ins_profile: [";
		for(i=0;i<gap_bound;i++){
			cout << params.at(ins_params_begin+i);
			if(i<gap_bound-1){
				cout << ";";
			}else{
				cout << "]" << endl;
			}
		}
		cout << endl << endl;
	}
	
	return 0;
	
}
