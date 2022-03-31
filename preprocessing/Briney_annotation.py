#
#  coding: utf-8
#
#  Briney_annotation.py
#
#  ---------------------------------------------------------------------------
#
#  Copyright (C) 2019-2022 Cosimo Lupo
#
#  This source code is distributed as part of the 'indie' software.
#  'indie' (INference on Deletion and InsErtions) is a versatile software
#  for evaluation and inference of indel and point substitution hypermutations
#  in high-throughput Ig antibody sequencing data, as well as
#  for the generation of synthetic repertoires with custom models.
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  You should have received a copy of the GNU General Public License
#  along with this program. If not, see <https://www.gnu.org/licenses/>.
#
#  ---------------------------------------------------------------------------
#
#  For any issue or question, please send an email to <cosimo.lupo89@gmail.com>.
#

import numpy as np
import pandas as pd
import glob
import os
import sys
import shutil
import multiprocessing as mp
from pandarallel import pandarallel
pd.set_option('display.max_columns',100)
from scipy.stats import poisson
from scipy.stats import binom
from scipy.optimize import curve_fit
from scipy.optimize import minimize
import Bio
from Bio import SeqIO
import datetime
import re

sys.path.insert(1, '/home/lupo/')
from annotation_funcs import natural_keys, make_csv_from_fasta, make_fasta_from_csv
from annotation_funcs import run_igBlast, parse_igBlast, revert_seq, reconstruct_gapped_seqs, write_anchored_seqs

########## settings ##########

preProcessFiles = False
runIgBlast = False
parseIgBlast = False
sortAndCohort = False

# For Briney data, we have:
# - Samples: ['316188','326650','326651','326713','326737','326780','326797','326907','327059']
# - Cell types: ['IgM', 'IgG']
# - Chain types: ['heavy']

splitSeqsFiles = True
subsetString = 'subset??'

species = "Homo_Sapiens"
cellType = "IgG"
chainType = "heavy"
UMIfilter = False
UMI_thr = 3
onlyAnnotatedCDR3 = True
V_len_min = 200

blast_database = "/home/lupo/igBlast/blast_database/Homo_Sapiens/"
home_dir = "/home/lupo/BCR_Briney_data/"
sample_dirs = glob.glob(home_dir + '3?????')

chainType_dict = {"heavy": "HC", "kappa": "KC", "lambda": "LC"}

########## defs ##########

def make_IDs(sample,cellType,chainType,row):
  if(row['SEQUENCE_INPUT']==row['SEQUENCE_INPUT']):
    seq_ID = sample + ':' + cellType + ':' + chainType + ':' + row['SEQUENCE_ID'] + ':N' + str(row['CONSCOUNT']) + ':D' + str(row['DUPCOUNT']) + ':R' + str(row.name)
    return seq_ID
  else:
    return np.nan

def remove_initial_Ns(row):
  if(row['SEQUENCE_INPUT']==row['SEQUENCE_INPUT']):
    seq = row['SEQUENCE_INPUT']
    while(seq[0]=='N'):
      seq = seq[1:]
    return seq
  else:
    return np.nan

########## main ##########

print('\n' + '\033[1m' + " ***** Pre-processing, annotation through igBlast and post-processing ***** " + '\033[0m' + '\n')

# Step 1: check which format seqs have to be read and (if needed) produce the other format

if(preProcessFiles):
  
  print("\nStep 1: check which format seqs have to be read and (if needed) produce the other format\n")
  
  sample_dirs.sort(key=natural_keys)
  for i,sample_dir in enumerate(sample_dirs):

    sample = sample_dir.split('/')[-1]
    print("    sample: " + sample)
    filenames = []

    # add patient name to file names
    files = glob.glob(sample_dir + '/' + '*.*')
    files.sort(key=natural_keys)
    for j,fullfilename in enumerate(files):
      filename = fullfilename.split('/')[-1]
      new_fullfilename = "/".join(fullfilename.split('/')[:-1]) + '/' + sample + '_' + cellType + '_' + chainType_dict[chainType] + '_' + fullfilename.split('/')[-1]
      if(filename.find(sample)==-1):
        os.rename(fullfilename,new_fullfilename)

    # extract ID and seqs from old files (and split them if necessary)
    out_file = sample_dir + '/' + sample + '_' + cellType + '_' + chainType_dict[chainType] + '.csv'
    if(os.path.isfile(out_file)==False):
      #in_file = sample_dir + '/' + sample + '_' + cellType + '_' + chainType_dict[chainType] + '_all_igblast_db-pass.tab'
      in_file = sample_dir + '/' + sample + '_' + cellType + '_' + chainType_dict[chainType] + '_all.tab'
      df = pd.read_csv(in_file, index_col=False, sep='\t')
      df = df[['SEQUENCE_INPUT','SEQUENCE_ID','CONSCOUNT','DUPCOUNT']]
      df['seq_ID'] = -1
      df['raw_seq_nt'] = np.nan
      pandarallel.initialize()
      df['seq_ID'] = df.parallel_apply(lambda row: make_IDs(sample,cellType,chainType_dict[chainType],row), axis=1)
      df['raw_seq_nt'] = df.parallel_apply(lambda row: remove_initial_Ns(row), axis=1)
      if(len(df)!=len(df['seq_ID'].unique())):
        print('ATTENTION: generated IDs are not unique! ' + str(len(df['seq_ID'].unique())) + ' vs ' + str(len(df)))
      else:
        print('OK: generated IDs are unique...')
      df = df[['seq_ID','raw_seq_nt']]
      if(splitSeqsFiles==False):
        df.to_csv(out_file, index=False, sep=';')
      else:
        s = 500000
        L = len(df)
        N = int(L/s) + 1
        #print('L: ',L)
        #print('N: ',N)
        for n in range(N):
          out_file = sample_dir + '/' + sample + '_' + cellType + '_' + chainType_dict[chainType] + '_subset' + '{:02d}'.format(n+1) + '.csv'
          df[s*n:s*(n+1)].to_csv(out_file, index=False, sep=';')
    
    # check csv
    if(splitSeqsFiles==False):
      csv_files = glob.glob(sample_dir + '/' + '*' + cellType + '_' + chainType_dict[chainType] + '.csv')
    else:
      csv_files = glob.glob(sample_dir + '/' + '*' + cellType + '_' + chainType_dict[chainType] + '_subset??' + '.csv')
    if(len(csv_files)>0):
      for j,fullfilename in enumerate(csv_files):
        filename = fullfilename.split('/')[-1].split('.')[0]
        filenames.append(filename)
  
    # check fasta
    if(splitSeqsFiles==False):
      fasta_files = glob.glob(sample_dir + '/' + '*' + cellType + '_' + chainType_dict[chainType] + '.fasta')
    else:
      fasta_files = glob.glob(sample_dir + '/' + '*' + cellType + '_' + chainType_dict[chainType] + '_subset??' + '.fasta')
    if(len(fasta_files)>0):
      for j,fullfilename in enumerate(fasta_files):
        filename = fullfilename.split('/')[-1].split('.')[0]
        filenames.append(filename)
    
    # make unique
    filenames = list(set(filenames))
    filenames.sort(key=natural_keys)
    
    # produce the missing format (if any)
    for i,filename in enumerate(filenames):
      fullfilename = sample_dir + '/' + filename
      if(os.path.isfile(fullfilename + '.csv')==True and os.path.isfile(fullfilename + '.fasta')==False):
        try:
          make_fasta_from_csv(fullfilename + '.fasta', headers=True, sep=';')
        except BaseException as err:
          print(err)
      elif(os.path.isfile(fullfilename + '.csv')==False and os.path.isfile(fullfilename + '.fasta')==True):
        try:
          make_csv_from_fasta(fullfilename + '.csv', headers=['seq_ID','raw_seq_nt'], sep=';')
        except BaseException as err:
          print(err)

# Step 2: run igBlast

if(runIgBlast):
  
  print("\nStep 2: run igBlast")

  sample_dirs.sort(key=natural_keys)
  for i,sample_dir in enumerate(sample_dirs):
    sample = sample_dir.split('/')[-1]
    print("    sample: " + sample)
    if(splitSeqsFiles==False):
      fasta_files = glob.glob(sample_dir + '/' + '*' + cellType + '_' + chainType_dict[chainType] + '.fasta')
    else:
      fasta_files = glob.glob(sample_dir + '/' + '*' + cellType + '_' + chainType_dict[chainType] + '_' + subsetString + '.fasta')
    if(len(fasta_files)>0):
      fasta_files.sort(key=natural_keys)
      for j,fullfilename in enumerate(fasta_files):
        in_file = fullfilename.split('.')[0] + '.fasta'
        out_file = fullfilename.split('.')[0] + '.igBlast_raw_output'
        if(splitSeqsFiles==True):
          subset = fullfilename.split('.')[0][-2:]
          print("        subset: " + subset)
          
        t1 = datetime.datetime.now()
        try:
          run_igBlast(in_file, species, chainType)
        except BaseException as err:
          print(err)
        t2 = datetime.datetime.now()
        print('            igBlast running time:', t2-t1)
      
# Step 3: parse the raw output from igBlast

if(parseIgBlast):
  
  print("\nStep 3: parse the raw output from igBlast")

  sample_dirs.sort(key=natural_keys)
  for i,sample_dir in enumerate(sample_dirs):
    sample = sample_dir.split('/')[-1]
    print("    sample: " + sample)
    if(splitSeqsFiles==False):
      igBlast_raw_output_files = glob.glob(sample_dir + '/' + '*' + cellType + '_' + chainType_dict[chainType] + '.igBlast_raw_output')
    else:
      igBlast_raw_output_files = glob.glob(sample_dir + '/' + '*' + cellType + '_' + chainType_dict[chainType] + '_' + subsetString + '.igBlast_raw_output')
    if(len(igBlast_raw_output_files)>0):
      igBlast_raw_output_files.sort(key=natural_keys)
      for j,fullfilename in enumerate(igBlast_raw_output_files):
        in_file = fullfilename.split('.')[0] + '.igBlast_raw_output'
        out_file = fullfilename.split('.')[0] + '.igBlast_statistics'
        if(splitSeqsFiles==True):
          subset = fullfilename.split('.')[0][-2:]
          print("        subset: " + subset)
        
        t1 = datetime.datetime.now()
        
        try:
          df = parse_igBlast(in_file, chainType, requireJ=True)
        except BaseException as err:
          print(err)
        
        # UMI counts
        df['UMI_counts'] = df['seq_ID'].apply(lambda x: int(x.split(':')[4][1:]))
        
        # Quality filtering
        df = df[df['V_best_align_length']>=V_len_min]    # V gene should align at least for V_len_min nt
        df = df[df['strand']=="+"]    # Only sequences read in the correct direction
        #df = df.query('UMI_counts'>@UMI_thr)    # Only UMIs counted at least UMI_thr times
        
        # mapping of nt sequences through IDs (included primers and C segment, if any)
        seqs_file = fullfilename.split('.')[0] + '.csv'
        raw_seqs = pd.read_csv(seqs_file, sep=';', low_memory=False, keep_default_na=True, index_col=False)
        raw_seqs['seq_ID'] = raw_seqs['seq_ID'].astype('str')
        df['raw_seq_nt'] = df['seq_ID'].map(raw_seqs.set_index('seq_ID')['raw_seq_nt'])
        df['seq_nt'] = df.apply(lambda row: row['raw_seq_nt'][row['V_best_align_start_seq']-1:row['J_best_align_end_seq']], axis=1)
        df['seq_nt_len'] = df['seq_nt'].apply(lambda x: len(x))
        
        # Clonal abundance
        counts = df.groupby('seq_nt').count().reset_index()
        counts = counts[['seq_nt','seq_ID']]
        counts.rename(columns={'seq_ID': 'seq_nt_counts'}, inplace=True)
        df['seq_nt_counts'] = df['seq_nt'].map(counts.set_index('seq_nt')['seq_nt_counts'])
        df = df.drop_duplicates(subset=['seq_nt'],keep='first')    # Drop seq_nt duplicates
        
        # Reconstruct gapped query and germlines
        df['temp'] = df.apply(lambda row: reconstruct_gapped_seqs(row), axis=1)
        df['gapped_query'] = df['temp'].apply(lambda x: x.split(';')[0])
        df['gapped_germline'] = df['temp'].apply(lambda x: x.split(';')[1])
        
        # Revert indels in the sequence
        df['reverted_seq_nt'] = df.apply(lambda row: revert_seq(row), axis=1)
        
        # Drop unnecessary data
        df = df.drop(['raw_seq_nt', 'V_best_aligned_query', 'V_best_aligned_germline', \
                      'J_best_aligned_query', 'J_best_aligned_germline', 'temp'], axis=1)
        
        # Export on file
        df.to_csv(out_file, index=False, sep=';')

        t2 = datetime.datetime.now()
        if(splitSeqsFiles==False):
          print('        parsing & formatting time:', t2-t1)
        else:
          print('            parsing & formatting time:', t2-t1)

# Step 4: sort NP-P into new files and gather cohortwide data

if(sortAndCohort):
  
  print("\nStep 4: sort NP-P into new files and gather data if split")
  
  sample_dirs.sort(key=natural_keys)
  for i,sample_dir in enumerate(sample_dirs):
    
    sample = sample_dir.split('/')[-1]
    print("    sample: " + sample + ' ' + cellType, end="")
    
    if(splitSeqsFiles==False):
      igBlast_statistics_files = glob.glob(sample_dir + '/' + '*' + cellType + '_' + chainType_dict[chainType] + '.igBlast_statistics')
    else:
      igBlast_statistics_files = glob.glob(sample_dir + '/' + '*' + cellType + '_' + chainType_dict[chainType] + '_subset*' + '.igBlast_statistics')
    igBlast_statistics_files.sort(key=natural_keys)
    
    if(len(igBlast_statistics_files)>0):

      # Open .igBlast_statistics file(s)
	  keys = ['seq_ID','seq_nt','V_J_frame','stop_codon_CDR3','CDR3_nt','N_indels','V_best_identity',\
			  'V_best_align_start_seq','V_best_align_end_seq','V_best_align_start_gene','V_best_align_end_gene',\
			  'UMI_counts']
      if(splitSeqsFiles==False):
        for j,fullfilename in enumerate(igBlast_statistics_files):
          in_file = fullfilename.split('.')[0] + '.igBlast_statistics'
          df = pd.read_csv(in_file, sep=';', low_memory=False, keep_default_na=True, index_col=False)
          df = df[keys]
      else:
        df = pd.DataFrame();
        for j,fullfilename in enumerate(igBlast_statistics_files):
          in_file = fullfilename.split('.')[0] + '.igBlast_statistics'
          df2 = pd.read_csv(in_file, sep=';', low_memory=False, keep_default_na=True, index_col=False)
          df2 = df2[keys]
          df = pd.concat([df,df2]).reset_index(drop=True)
      
      # Drop duplicate sequences
      df = df.drop_duplicates(subset=['seq_nt'],keep='first').reset_index(drop=True)    # Drop seq_nt duplicates

      print("  [data_size: " + str(len(df)) + "]")
      df['cellType'] = cellType
      df['sample'] = sample
      
      # Trim sequences at the two edges
      n_l = 3    # N. of nucleotides trimmed on the left (Nmer size is N = 2*n+1, I'm also taking into account the multialignment between the V germlines of +- 3)
      n_r = 3    # N. of nucleotides trimmed on the right (Nmer size is N = 2*n+1, the J germlines end at the same position)
      df['trimmed_seq_nt'] = df['seq_nt'].apply(lambda x: x[n_l:-n_r])
      
      # Keep only sequences above UMI count threshold
      if(UMIfilter):
        df = df[df['UMI_counts']>=UMI_thr]
      
      # Keep only sequences with annotated CDR3
      if(onlyAnnotatedCDR3):
        df = df.query('CDR3_nt==CDR3_nt')
      
      for filterOutHyperIndels in [False,True]:
        for trimEdges in [False,True]:
          
          # Filter P/NP sequences
          # Here I'm using the following definition:
          # - P:  Both V,J are in-frame and no stop codon in the CDR3
          # - NP: Either V,J are out-frame or there is a stop codon in the CDR3
          #df_NP = df.loc[(df['V_J_frame']=="Out-of-frame") | (df['stop_codon_CDR3']=="Yes")]
          df_NP = df.loc[df['V_J_frame']=="Out-of-frame"]
          df_P = df[(df['V_J_frame']=="In-frame") & (df['stop_codon_CDR3']=="No")]
          
          # Keep only sequences with no hyper-indels
          if(filterOutHyperIndels):
            df_NP = df_NP.query('N_indels==0')
            df_P = df_P.query('N_indels==0')
          
          # Export sequences on file
          if(splitSeqsFiles==False):
            out_file = fullfilename.split('.')[0]
            out_file_NP = out_file + '_uniqueSeqs_NP'
            out_file_P = out_file + '_uniqueSeqs_P'
          else:
            out_file = fullfilename.split('.')[0]
            out_file = "_".join(out_file.split('_')[:-1])
            out_file_NP = out_file + '_uniqueSeqs_NP'
            out_file_P = out_file + '_uniqueSeqs_P'
          if(UMIfilter):
            out_file_NP += '_UMIthr_' + str(UMI_thr)
            out_file_P += '_UMIthr_' + str(UMI_thr)
          if(trimEdges):
            out_file_NP += '_trimmed_' + str(n_l) + '_' + str(n_r)
            out_file_P += '_trimmed_' + str(n_l) + '_' + str(n_r)
			keys2 = ['seq_ID','trimmed_seq_nt','V_best_identity','V_best_align_start_seq','V_best_align_end_seq','V_best_align_start_gene','V_best_align_end_gene']
            df_NP = df_NP[keys2]
            df_P = df_P[keys2]
			df_NP = df_NP.rename(columns={'trimmed_seq_nt': 'seq_nt'})
            df_P = df_P.rename(columns={'trimmed_seq_nt': 'seq_nt'})
          else:
			keys2 = ['seq_ID','seq_nt','V_best_identity','V_best_align_start_seq','V_best_align_end_seq','V_best_align_start_gene','V_best_align_end_gene']
            df_NP = df_NP[keys2]
            df_P = df_P[keys2]
          if(filterOutHyperIndels):
            out_file_NP += '_noHyperIndels'
            out_file_P += '_noHyperIndels'
          
          # Write csv
          df_NP[['seq_ID','seq_nt']].to_csv(out_file_NP + '.csv', index=False, sep=';')
          df_P[['seq_ID','seq_nt']].to_csv(out_file_P + '.csv', index=False, sep=';')
          
          # Write fasta
          ofile = open(out_file_NP + '.fasta', "w")
          for ii in range(len(df_NP)):
            ofile.write(">" + df_NP.iloc[ii,0] + "\n" + df_NP.iloc[ii,1] + "\n")
          ofile.close()
          ofile = open(out_file_P + '.fasta', "w")
          for ii in range(len(df_P)):
            ofile.write(">" + df_P.iloc[ii,0] + "\n" + df_P.iloc[ii,1] + "\n")
          ofile.close()
          
          if(filterOutHyperIndels==False and trimEdges==False):
            
            # Formatted csv for indels software
            # When numbering starting from 0:
            #  - initial position is given by 'start' index - 1
            #  - final position (not included) is given by 'end' index, so that last included is given by 'end' index - 1
            
            anchored_subset_size = 100000
			
			for produc in ['NP','P']:
            
			  out_file = (produc == 'NP' ? out_file_NP : out_file_P)
			  out_file += '_anchored.csv'
			  out_f = open(out_file, "w")
              out_f.write("seq_ID;aligned_seq_nt;V_best;V_best_start;V_best_end\n")
			  df2 = (produc == 'NP' ? df_NP : df_P)
              for r,row in df2.iterrows():
                out_f.write(write_anchored_seqs(row) + "\n")
              out_f.close()
            
              if(len(df2)>anchored_subset_size):
                out_file = (produc == 'NP' ? out_file_NP : out_file_P)
				out_file += '_anchored_subset100K.csv'
				out_f = open(out_file, "w")
                out_f.write("seq_ID;aligned_seq_nt;V_best;V_best_start;V_best_end\n")
                df2 = df2.sample(n=anchored_subset_size).sort_index().reset_index(drop=True)
                for r,row in df2.iterrows():
                  out_f.write(write_anchored_seqs(row) + "\n")
                out_f.close()

quit()
