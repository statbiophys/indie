# *indie* (INference on Deletions and InsErtions) documentation
----

## Cosimo Lupo &#169; 2019-2021

Probabilistic alignment software to deal with indels in Ig receptor sequences.

## Version

Latest released version: 1.2.3 (dated 14/12/2021).

## License

Free use of *indie* is granted under the terms of the GNU General Public License version 3 (GPLv3).

## Contact

For any issue, question or bug, please write [us](mailto:cosimo.lupo89@gmail.com) an email.

## Before the compilation

Before proceeding to the compilation, make sure your platform is equipped with a C compiler supporting OpenMP. The GNU C compiler (ie gcc) is recommended, and usually already present on Linux systems. If on MacOs, please be careful and read the section below before going on with the compilation. *indie* has not yet been tested on Windows.

No libraries are required to be installed separately; all the necessary files are already shipped with *indie* distribution.

#### ⚠️ Mac OS

Mac OS platforms are now equipped with clang compiler (installed together with Xcode application), aliased as gcc/g++. This means that when using the Makefile on Mac Os, the g++ invocation will actually refer to clang. Unfortunately, the latter does not support OpenMP, so the compilation of *indie* will fail.

Recommended procedure on Mac OS platforms is hence to install the GNU C compiler through [Macports](https://www.macports.org) or [Homebrew](https://brew.sh) (to be installed in turn) and then use it for compiling *indie*.

It's a good practice to avoid assigning the same alias `gcc` to the freshly installed GNU compiler as Apple clang; otherwise, [many potential failures](https://docs.brew.sh/Custom-GCC-and-cross-compilers.html) may occur during compilation of other packages.

## Compiling *indie*

Download the latest version of the *indie* source code from [Github](https://github.com/statbiophys/indie) and extract the files from the archive, if not already done automatically by the OS.

Copy the directory containing the source files in the desired position and then proceed to the compilation:

```
make
```

Makefile invokes by default `g++` compiler; if such command/alias actually refers to a compiler that does not support OpenMP, explicit reference has to be passed when using the make command through `CXX` flag. As an example, if wanting to use a compiler named `g++-9`:

```
make CXX=g++-9
```

*indie* executable then appears in the current directory.

No full installation is currently provided, so that indie executable and supplied models won't be copied in the system default directory. Every invocation to the *indie* software should be done in the folder where the executable actually is, together with `models` directory.

## Using *indie*

*indie* has several command-line options, that can be listed by using the `-h` or `--help` option:

```
./indie -h
```

Here a summary:

```
List of command-line options

  -a, --alpha         [req. arg.] sets the global learning rate parameter to the specified value
  -b, --batch_name    [req. arg.] creates a directory with the specified name
                      and puts all the output files in it
  -g, --generate      [req. arg.] generate N synthetic sequences (where N has to be inserted
                      as argument
  -h, --help          [no arg.] prints all the possible command-line options
                      and if they require any (optional) argument
  -i, --inference     [no arg.] performs the inference on model parameters
                      by maximizing the total alignment likelihood
  -l, --log           [no arg.] prints on auxiliary files optional info about the convergence
                      of the algorithm during the inference step
  -p, --parameters    [req. arg.] reads from the specified file the parameters to be used
                      for the generation step or as the initial condition for the inference step
  -r, --read          [req. arg.] reads from the specified file the dataset
                      of sequences to analyze
  -T, --max_iter      [req. arg.] sets the maximum number of gradient descent iterations
```

All the options need to be specified at the same time, with a single invocation to `./indie`.

A typical command-line string to read from a file, perform a deterministic annotation and then infer the model parameters would be something like:

```
./indie -b wk_dir -r input_file.csv -i
```

where the flag `-b` specifies the working directory where to store output files. If such directory is not yet existing, it will be created, together with the necessary intermediate folders along the path.

Instead of reading an existing file, a custom synthetic dataset of size N can be generated, and then the inference can be launched on it as well by combining the `-g` flag with the `-i` one:

```
./indie -b wk_dir -g N -i
```

Further details about the various options can be found below.

### Reading from file

*indie* can be executed on an existing dataset, with input file to be specified by means of the `-r` flag

```
./indie -r input_file.csv
```

This file has to be formatted according to a specific comma-separated syntax:

| seq\_ID | aligned\_seq\_nt                  | V\_best       | V\_best\_start | V\_best\_end |
| :------ | :-------------------------------- | :------------ | :------------- | :----------- |
| 0       | CAGGTGCAGCTGGTGGAGTCTGGGGGGGGC... | IGHV3\-30\*12 | 0              | 296          |
| 1       | GAGGTGCAGCTGGTGGAGTCTGGGGGAGGC... | IGHV3\-48\*04 | 0              | 296          |
| 2       | TAGGTGAAGCTCGCCGAGGTGAAGAAGCGT... | IGHV1\-2\*01  | 0              | 296          |
| 3       | CAGGTGCAGTTGCAGAGTTCGCTCCCAGGT... | IGHV4\-61\*03 | 0              | 299          |
| ...     | ...                               | ...           | ...            | ...          |

Last columns refer to the name of the best-scoring V template and the limits of the aligned region along the template itself. Any Ig alignment tool could be used to obtain these info; we recommend [igBlast](https://www.ncbi.nlm.nih.gov/igblast) for its reliability and speed, but other tools can be exploited as well.

Once read the input file, *indie* creates a folder `read` and stores in it a copy of the input file.

### Generating a custom dataset

Instead of reading sequences from an existing dataset, the `generate` feature of *indie* allows to create a synthetic dataset of human Ig sequences, according to a model fully customable. Default generative model is stored in the `models` folder, and can be invoked to create a repertoire of size N through the command:

```
./indie -g N
```

A specific model parameters can be used by means of the `-p` flag:

```
./indie -g N -p custom_parms.txt
```

Parameters file syntax has to be the following:

```
@Rates
#frac_mutated
0.90
#aver_err_rate;del_ratio;ins_ratio
0.10;0.02;0.02
@Profiles
#max_gap_length
50
#del_profile
0.095808132;0.086690783;0.078441064;0.07097641;0.064222112;0.05811057;0.052580618;0.04757691;0.043049369;0.03895268;0.035245842;0.031891757;0.028856855;0.026110762;0.023625994;0.021377684;0.019343328;0.017502567;0.015836978;0.01432989;0.012966221;0.011732322;0.010615844;0.0096056125;0.0086915177;0.0078644104;0.0071160128;0.0064388346;0.0058260985;0.0052716719;0.004770006;0.0043160799;0.0039053506;0.0035337074;0.0031974307;0.0028931549;0.0026178348;0.0023687149;0.0021433019;0.0019393397;0.0017547871;0.0015877971;0.0014366982;0.0012999783;0.001176269;0.0010643322;0.00096304761;0.00087140151;0.00078847669;0.00071344321
#ins_profile
0.095808132;0.086690783;0.078441064;0.07097641;0.064222112;0.05811057;0.052580618;0.04757691;0.043049369;0.03895268;0.035245842;0.031891757;0.028856855;0.026110762;0.023625994;0.021377684;0.019343328;0.017502567;0.015836978;0.01432989;0.012966221;0.011732322;0.010615844;0.0096056125;0.0086915177;0.0078644104;0.0071160128;0.0064388346;0.0058260985;0.0052716719;0.004770006;0.0043160799;0.0039053506;0.0035337074;0.0031974307;0.0028931549;0.0026178348;0.0023687149;0.0021433019;0.0019393397;0.0017547871;0.0015877971;0.0014366982;0.0012999783;0.001176269;0.0010643322;0.00096304761;0.00087140151;0.00078847669;0.00071344321
```

`@Rates` section refers to the frequency of hyper-mutations events:

- `frac_mutated`: the probability for a generated Ig sequence to be naive, ie not mutated at all;
- `aver_err_rate`: the average probability per base pair of having a point substitution (in fact, such uniform rate is only used in the deterministic annotation of indels, but not to generate synthetic sequences, see below);
- `del_ratio`: the base pair probability of hosting a deletion at a given position, rescaled by the point substitution rate of that sequence;
- `ins_ratio`: the base pair probability of hosting an insertion at a given position, rescaled by the point substitution rate of that sequence.

`@Profiles` section instead codes for the length distribution of indels:

- `max_gap_length`: it's the maximum length allowed for a single deletion/insertion event during the *indie* execution; it has to be smaller than (or equal to) the maximum value of 60 base pairs hard-coded in the software;
- `del_profile `: it's the normalized vector of probabilities ruling the length of deletions; first component is the probability of deleting one base pair, second component is the probability of deleting two base pairs, and so on; the length of such vector has be `max_gap_length`;
- `ins_profile `: same as `del_profile `, but for insertions.

⚠️ Other features of the generative model are hard-coded in the software and can't be simply changed through the command line. Among them, the set of V templates to pick from (default is `models/imgt_templates_IGHV_Homo_Sapiens_F.csv`), the corresponding frequency for each V family (default is `models/V_usage.csv`), and the repertoire distribution of point mutation rates (default is `models/prior.csv`). Please be careful when trying to modify such files, or directly the source code relative to the generation (and in this latter case, remember to recompile at the end of changes).

Generated sequences are stored in a standard fasta file (`synthetic_seqs.fasta`) under the `generate` folder, together with a formatted file (`synthetic_seqs_anchored.csv`) with the same syntax required for input files, a copy of the model parameters file (`generative_params.txt`), a copy of the file containing the point mutation distribution (`generative_mu_err_distr.txt`), and finally a file with the set of events behind the generation of each synthetic sequence, named `synthetic_scenarios.txt`:

| seq\_ID | mu\_err  | V\_choice | error\_positions  | error\_list               | deletion\_pos | deletion\_list | insertion\_pos | insertion\_list |
| :------ | :------- | --------- | :---------------- | :------------------------ | :------------ | :------------- | :------------- | :-------------- |
| 0       | 0.054128 | 112       | [34,58,70,74,...] | [T->C,T->C,C->T,T->A,...] | []            | []             | [289]          | [6->TGGGCG]     |
| 1       | 0.000359 | 138       | [195]             | [G->C]                    | []            | []             | []             | []              |
| 2       | 0.133590 | 4         | [0,6,23,26,...]   | [C->T,C->A,G->C,T->C,...] | [11]          | [12]           | []             | []              |
| 3       | 0.139510 | 244       | [9,15,16,17,...]  | [C->T,G->A,A->G,G->T,...] | [211,261]     | [19,1]         | []             | []              |
| ...     | ...      | ...       | ...               | ...                       | ...           | ...            | ...            | ...             |

### Deterministic alignment

After having read the input file (or having generated a synthetic dataset) *indie* software moves to performing a deterministic alignment against V templates suggested in the suitably formatted input file.

Alignment is done through the standard Needleman-Wunsch (NW) algorithm, and results are stored in a new folder called `greedy`. Model parameters used for the initialization cabn be found in `greedy/initial_params.txt`, while the ones resulting from the deterministic alignment are in `greedy/greedy_params.txt`. More info about the best-scoring alignment for each sequence can be found in `greedy/greedy_alignment_summary.csv`.

### Probabilistic alignment \& model inference

Once performed the deterministic annotation of the Ig dataset, corresponding results should be compared with those inferred within the probabilistic framework provided by *indie*. To this aim, it's necessary to provide the flag `-i` when invoking `./indie`. By command line can be also changed the maximum number of inference steps, `-T`, and the learning rate alpha, `-a`.

Inference results will be stored in a dedicated folder, `inference`. Again, initial and final values for model parameters can be found in `inference/initial_params.txt` and `inference/inferred_params.txt`, respectively.

Using the flag for a detailed report of the inference step by step, `-l`, will add further files, `inference/aux_infer_file_attempt_1.txt`, `inference/aux_posterior_file_attempt_1.txt` and `inference/aux_prior_file_attempt_1.txt`, with extra info about the values of model parameters, repertoire and sequence-specific point mutation distributions at each inference time step. The integer number at the end of the file names stands for the number of times the inference procedure has been launched. According to the examined dataset, indeed, it may be necessary to decrease the initial learning rate and hence restart the inference a few times.
