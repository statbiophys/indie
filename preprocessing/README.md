# Auxiliary scripts for the *indie* paper
----

## Cosimo Lupo, Natanael Spisak &#169; 2019-2022

Preprocessing and annotation scripts for data used in the *indie* paper: Cosimo Lupo, Natanael Spisak, Aleksandra M. Walczak, Thierry Mora, *"Learning the statistics and landscape of somatic mutation-induced insertions and deletions in antibodies"*, [arXiv (2021) 2112.07953](https://arxiv.org/abs/2112.07953).

## Requirements

Below a list of third-party softwares used inside our scripts for raw reads preprocessing, data formatting and annotation. Please be sure to have them properly installed before running our scripts.

- [pRESTO](https://presto.readthedocs.io/en/stable/), for preprocessing raw reads from high-throughput sequencing of B cell and T cell repertoires; it is part of the [Immcantation](https://immcantation.readthedocs.io/en/stable/) framework. We relied on software version **0.5.13**.
- [Seqtk](https://github.com/lh3/seqtk), for fastq to fasta formatting.
- [igBlast](https://ncbi.github.io/igblast/), for the annotation of B cell and T cell repertoire sequences. We relied on software version **1.13.0**.

All the following scripts are meant to be launched from command line in a bash shell. While some of them are directly written in bash, others are written in [Python3](https://www.python.org); be sure to have it installed as well.

## License

Free use of present scripts is granted under the terms of the GNU General Public License version 3 (GPLv3). Credits for required third-party softwares to their respective owners.

## Contact

For any issue, question or bug, please write [us](mailto:cosimo.lupo89@gmail.com) an email.

## Preprocessing from raw reads

The *indie* software is applied on a publicly available IgM and IgG repertoire from a recent ultra-deep repertoire study of Immunoglobulin heavy chains (IGH), including approximately 3 billions sequences in total, coming from 9 different healthy donors ([Briney *et al*, 2019](https://www.nature.com/articles/s41586-019-0879-y)).

We obtained IgM and IgG receptor sequences by reprocessing raw reads, available on the [NCBI Sequencing Read Archive](https://www.ncbi.nlm.nih.gov/sra) with access number PRJNA406949. The primers for the V and C regions are provided in the *Supplementary Information* section of the original article introducing the dataset.

The preprocessing script `preprocessing.sh` can be modified by specifying directories for raw read files (`R1.fastq` and `R2.fastq`) and for primers in the first lines of the script. Then, can be executed through the command:

```
./preprocessing.sh
```

As output, `.fastq` and `.fasta` files segregated by read counts per UMI will be produced. We set 3 as threshold, but this parameter can be changed in the script, too. Tab-separated files will also be produced, with extra information inside (e.g. number of reads per UMI).

## Annotation and quality filtering

Once `.fasta` files with BCR sequences are ready, they can be annotated, quality-filtered and eventually sorted into in-frame and out-of-frame subsets.

Main script in charge of these tasks is `Briney_annotation.py `, that can be invoked from the command line:

```
python3 Briney_annotation.py
```

The script is modular and quite flexible, with several customizable parameters in the first lines of the script: data and output directories, type of Ig chain, quality filters, and so on.

There are 4 main execution blocks, that should be launched sequentially (or even independently by acting on the related boolean variables, provided the necessary input files are already present):

1. A preliminar section with some further manipulation on sequence files (split large fasta files into chunks, create new sequence IDs, and so on).
2. The proper annotation stage, where igBlast software (version 1.13.0) is exploited and files with igBlast raw output are produced, starting from plain `.fasta` files. Reference templates have been extracted by [IMGT](https://www.imgt.org), considering only functional (F) genes, and including different alleles for each gene.
3. A parsing of the `.igBlast_raw_output` files obtained in the previous stage, where new features for each sequence are computed (e.g. the list of mutations, the unmutated sequence ancestor, ...) and sequences are also quality-filtered (we required V annotation for at least 200 nt counted along the germline, J annotation and correct 5' to 3' direction of read). Output `.igBlast_statistics` files are formatted as `.csv` files, separated by `;`.
4. Starting from annotation results, sequences are segregated into in-frame and out-of-frame sequences, possibly applying further filters (e.g. keeping only sequences with a properly annotated CDR3, or above a minimum number of reads per UMI, and so on). Final sequences are written into both `.fasta` and `.csv` files.

It's in this last block that sequences are also formatted into the peculiar comma-separated files (with semicolons as field delimiters) required by the *indie* software:

| seq\_ID | aligned\_seq\_nt                  | V\_best       | V\_best\_start | V\_best\_end |
| :------ | :-------------------------------- | :------------ | :------------- | :----------- |
| 0       | CAGGTGCAGCTGGTGGAGTCTGGGGGGGGC... | IGHV3\-30\*12 | 0              | 296          |
| 1       | GAGGTGCAGCTGGTGGAGTCTGGGGGAGGC... | IGHV3\-48\*04 | 0              | 296          |
| 2       | TAGGTGAAGCTCGCCGAGGTGAAGAAGCGT... | IGHV1\-2\*01  | 0              | 296          |
| 3       | CAGGTGCAGTTGCAGAGTTCGCTCCCAGGT... | IGHV4\-61\*03 | 0              | 299          |
| ...     | ...                               | ...           | ...            | ...          |

where last columns refer to the name of the best-scoring V template and the limits of the aligned region along the template itself.
