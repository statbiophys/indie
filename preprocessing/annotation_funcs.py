#
#  coding: utf-8
#
#  annotation_funcs.py
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
import os
import pandas as pd
import Bio
from Bio import SeqIO
import re

########## settings ##########

blast_database = "/home/lupo/igBlast/blast_database/"
  
########## defs ##########

def make_csv_from_fasta(in_file, headers=[], sep=';'):
  
  """
  Arguments:
    - in_file:    [string] is the full address of the csv file to be converted into fasta
    - headers:    [array of strings, optional] is the full address of the csv file to be
                  converted into fasta, default value is an empty array
    - sep:        [string, optional] is the field separator for the csv file, default value is ';'
  Returns:
    - void
  """
  
  func_name = 'make_csv_from_fasta'
  
  if(not os.path.isfile(in_file)):
    raise ValueError(' *** Error in the \'' + func_name + '\' function. Input file does not exist. *** ')
  elif(in_file.split('.')[-1]!='fasta'):
    raise ValueError(' *** Error in the \'' + func_name + '\' function. Input file does not have the required \'fasta\' extension. *** ')
  else:
    out_file = ".".join(in_file.split('.')[:-1]) + '.csv'
  
  fp = open(out_file, "w")
  if headers != []:
    fp.write(sep.join(headers) + "\n")
  for record in SeqIO.parse(in_file, 'fasta'):
    fp.write(record.id + sep + str(record.seq.upper()) + "\n")
  fp.close()
  
  return

def make_fasta_from_csv(in_file, headers=True, sep=';'):
  
  """
  Arguments:
    - in_file:    [string] is the full address of the fasta file to be converted into csv
    - headers:    [boolean, optional] tells if the csv has explicit headers or not,
                  default value is True (i.e. first row is skipped)
    - sep:        [string, optional] is the field separator for the csv file, default value is ';'
  Returns:
    - void
  """
  
  func_name = 'make_fasta_from_csv'
  
  if(not os.path.isfile(in_file)):
    raise ValueError(' *** Error in the \'' + func_name + '\' function. Input file does not exist. *** ')
  elif(in_file.split('.')[-1]!='csv'):
    raise ValueError(' *** Error in the \'' + func_name + '\' function. Input file does not have the required \'csv\' extension. *** ')
  elif(headers not in [True, False]):
    raise ValueError(' *** Error in the \'' + func_name + '\' function. \'headers\' argument has to be a boolean. *** ')
  else:
    out_file = ".".join(in_file.split('.')[:-1]) + '.fasta'
    if(headers):
      head = 0
    else:
      head = None
  
  fp = open(out_file, "w")
  df = pd.read_csv(in_file, index_col=False, header=head, sep=sep)
  for i in df.index:
    fp.write(">" + str(df.iat[i,0]) + "\n" + str(df.iat[i,1]).upper() + "\n")
  fp.close()
  
  return

def pairwise_comparison(query,germline,context):
	
    """
    Arguments:
        - query [string]: is the gapped query, as aligned by igBlast
        - germline [string]: is the gapped germline, as aligned by igBlast
        - context [integer]: is the size of the context for hyper-indels
    Returns:
        - a string containing a list of errors (point mutations, deletions and insertions)
    """
	
    L = len(query)
    Lside = context
    Rside = context
    fixed_seq = ""
    indel_count = 0
    ins_count = 0
    ins_flag = 0
    del_count = 0
    del_flag = 0
    df = pd.DataFrame({'type': [], 'composition': [], 'Lcontext': [], 'Rcontext': []})
    error_list = []
    error_str = ""
    for i in range(L):
        if(query[i]=='-'):
            # deletion
            if(del_flag==0):
                dct = {'type': "del", 'composition': "", 'Lcontext': str(query[max(0,i-Lside):i]), 'Rcontext': ""}
                df = df.append(dct, ignore_index=True)
                error_str = "[d|" + str(i) + "|" + str(query[max(0,i-Lside):i]) + "|"
            del_flag += 1
            df.loc[indel_count]['composition'] = df.iloc[indel_count]['composition'] + germline[i]
            error_str += germline[i]
        elif(germline[i]=='-'):
            # insertion
            if(ins_flag==0):
                dct = {'type': "ins", 'composition': "", 'Lcontext': str(query[max(0,i-Lside):i]), 'Rcontext': ""}
                df = df.append(dct, ignore_index=True)
                error_str = "[i|" + str(i) + "|" + str(query[max(0,i-Lside):i]) + "|"
            ins_flag += 1
            df.iloc[indel_count]['composition'] = df.iloc[indel_count]['composition'] + query[i]
            error_str += query[i]
        else:
            # match/mismatch
            if(ins_flag>0):
                df.iloc[indel_count]['Rcontext'] = query[i:min(L,i+Rside)]
                indel_count += 1
                ins_flag = 0
                error_str += "|" + str(query[i:min(L,i+Rside)]) + "]"
                error_list.append(error_str)
                error_str = ""
            elif(del_flag>0):
                df.iloc[indel_count]['Rcontext'] = query[i:min(L,i+Rside)]
                indel_count += 1
                del_flag = 0
                error_str += "|" + str(query[i:min(L,i+Rside)]) + "]"
                error_list.append(error_str)
                error_str = ""
            if(query[i]!=germline[i]):
                error_str = "[e|" + str(i) + "]"
                error_list.append(error_str)
                error_str = ""
    List = []
    for i in range(len(df.index)):
        List.append([df.iloc[i]['type'],df.iloc[i]['composition'],df.iloc[i]['Lcontext'],df.iloc[i]['Rcontext']])
    if(len(error_list)>0):
      err_str = '['
      for el in error_list:
          err_str += str(el) + ','
      err_str = err_str[:-1]
      err_str += ']'
      return err_str
    else:
      return "[]"

def nt2aa(ntseq):
	
    """
	Translate a nucleotide sequence into an amino acid sequence.

    Parameters
    ----------
    ntseq : str
        Nucleotide sequence composed of A, C, G, or T (uppercase or lowercase)

    Returns
    -------
    aaseq : str
        Amino acid sequence
    
    Example
    --------
    >>> nt2aa('TGTGCCTGGAGTGTAGCTCCGGACAGGGGTGGCTACACCTTC')
    'CAWSVAPDRGGYTF'
        
    """
    nt2num = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'a': 0, 'c': 1, 'g': 2, 't': 3}
    aa_dict ='KQE*TPASRRG*ILVLNHDYTPASSRGCILVFKQE*TPASRRGWMLVLNHDYTPASSRGCILVF'
    return ''.join([aa_dict[nt2num[ntseq[i]] + 4*nt2num[ntseq[i+1]] + 16*nt2num[ntseq[i+2]]] for i in range(0, len(ntseq), 3) if i+2 < len(ntseq)])

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

def atoi(text):
  return int(text) if text.isdigit() else text

def natural_keys(text):
	
  """
  alist.sort(key=natural_keys) sorts in human order
  http://nedbatchelder.com/blog/200712/human_sorting.html
  (See Toothy's implementation in the comments)
  """
  return [ atoi(c) for c in re.split(r'(\d+)', text) ]

def find_templates(species, chainType):
	
  """
  Arguments:
    - species:    [string] is the species of interest;
                  allowed values are: ["Homo_Sapiens"]
    - chainType:  [string] is the type of chain to be analyzed by igBlast;
                  allowed values are: ["heavy","HC","kappa","KC","lambda","LC"]
  Returns:
    - a triple with V, D and J databases
  """
  
  func_name = 'find_templates'
  allowed_species = ["Homo_Sapiens"]
  allowed_chains = ["heavy","HC","kappa","KC","lambda","LC"]
  
  if(species not in allowed_species):
    raise ValueError(' *** Error in the \'' + func_name + '\' function. Argument \'species\' does not match any of the allowed values. *** ')
  elif(chainType not in allowed_chains):
    raise ValueError(' *** Error in the \'' + func_name + '\' function. Argument \'chainType\' does not match any of the allowed values. *** ')
  
  if(chainType in ["heavy","HC"]):
    Vdatabase = blast_database + species + "/BCR_Heavy/forIgBlast_IGHV_" + species + "_F.fasta"
    Ddatabase = blast_database + species + "/BCR_Heavy/forIgBlast_IGHD_" + species + "_F.fasta"
    Jdatabase = blast_database + species + "/BCR_Heavy/forIgBlast_IGHJ_" + species + "_F.fasta"
  elif(chainType in ["kappa","KC"]):
    Vdatabase = blast_database + species + "/BCR_Kappa/forIgBlast_IGKV_" + species + "_F.fasta"
    Ddatabase = blast_database + species + "/BCR_Heavy/forIgBlast_IGHD_" + species + "_F.fasta"
    Jdatabase = blast_database + species + "/BCR_Kappa/forIgBlast_IGKJ_" + species + "_F.fasta"
  elif(chainType in ["lambda","LC"]):
    Vdatabase = blast_database + species + "/BCR_Lambda/forIgBlast_IGLV_" + species + "_F.fasta"
    Ddatabase = blast_database + species + "/BCR_Heavy/forIgBlast_IGHD_" + species + "_F.fasta"
    Jdatabase = blast_database + species + "/BCR_Lambda/forIgBlast_IGLJ_" + species + "_F.fasta"
    
  return Vdatabase,Ddatabase,Jdatabase

def run_igBlast(in_file, species, chainType):
  
  """
  Arguments:
    - in_file:    [string] is the full address of the fasta file to pass through igBlast
    - chainType:  [string] is the type of chain to be analyzed by igBlast;
                  allowed values are: ["heavy","HC","kappa","KC","lambda","LC"]
  Returns:
    - void
  """
  
  func_name = 'run_igBlast'
  
  if(not os.path.isfile(in_file)):
    raise ValueError(' *** Error in the \'' + func_name + '\' function. Input file does not exist. *** ')
  elif(in_file.split('.')[-1]!='fasta'):
    raise ValueError(' *** Error in the \'' + func_name + '\' function. Input file does not have the required \'fasta\' extension. *** ')
  
  out_file = ".".join(in_file.split('.')[:-1]) + '.igBlast_raw_output'
  
  try:
    Vdatabase,Ddatabase,Jdatabase = find_templates(species, chainType) # 
  except BaseException as err:
    raise ValueError(err)
  
  org = {"Homo_Sapiens": "human"}
  
  igBlast_command = "igblastn" + \
                    " -germline_db_V " + Vdatabase + \
                    " -germline_db_D " + Ddatabase + \
                    " -germline_db_J " + Jdatabase + \
                    " -organism " + org[species] + \
                    " -query " + in_file + \
                    " -auxiliary_data optional_file/" + org[species] + "_gl.aux" + \
                    " -outfmt '7 std qseq sseq'" + \
                    " -show_translation" + \
                    " > " + out_file
  os.system(igBlast_command)
  
  return

def parse_igBlast(in_file, chainType, requireJ=True):
  
  """
  Arguments:
    - in_file:    [string] is the full address of the igBlast raw output file to be parsed
    - chainType:  [string] is the type of chain to be analyzed by igBlast;
                  allowed values are: ["heavy","HC","kappa","KC","lambda","LC"]
  Returns:
    - a pandas dataframe, containing all the info for each sequence
  """
  
  func_name = 'parse_igBlast'
  allowed_chains = ["heavy","HC","kappa","KC","lambda","LC"]
  
  if(not os.path.isfile(in_file)):
    raise ValueError(' *** Error in the \'' + func_name + '\' function. Input file does not exist. *** ')
  elif(chainType not in allowed_chains):
    raise ValueError(' *** Error in the \'' + func_name + '\' function. Argument \'chainType\' does not match any of the allowed values. *** ')
  elif(requireJ not in [True, False]):
    raise ValueError(' *** Error in the \'' + func_name + '\' function. \'requireJ\' argument has to be a boolean. *** ')
  
  dict_keys = ['seq_ID', 'V_best_identity', 'D_best_identity', 'J_best_identity', 'stop_codon', \
               'V_J_frame', 'productive', 'strand', 'CDR3_nt', 'V_CDR3len', 'V_best_frac_of_matches', \
               'V_best_align_length', 'V_best_align_length_beforeCDR3', 'V_best_N_of_mm', \
               'V_best_gap_opens', 'V_best_total_gap_length', \
               'V_best_align_start_seq', 'V_best_align_end_seq', 'V_best_align_start_gene', 'V_best_align_end_gene', \
               'V_best_igBlast_score', 'V_best_aligned_query', 'V_best_aligned_germline', \
               'D_best_frac_of_matches', 'D_best_align_length', 'D_best_N_of_mm', 'D_best_gap_opens', \
               'D_best_total_gap_length', 'D_best_align_start_seq', 'D_best_align_end_seq', \
               'D_best_align_start_gene', 'D_best_align_end_gene', 'D_best_igBlast_score', 'D_best_aligned_query', \
               'D_best_aligned_germline', 'J_best_frac_of_matches', 'J_best_align_length', \
               'J_best_N_of_mm', 'J_best_gap_opens', 'J_best_total_gap_length', \
               'J_best_align_start_seq', 'J_best_align_end_seq', 'J_best_align_start_gene', 'J_best_align_end_gene', \
               'J_best_igBlast_score', 'J_best_aligned_query', 'J_best_aligned_germline']
  
  row_list = []
  df = pd.DataFrame()
  seq_count = 0
  
  with open(in_file, 'r') as temp_f:
    ParseType = 0
    for line in temp_f:
      line = line.strip()
      
      # Step 1
      if(line.find("# Query:")!=-1 and ParseType==0):
        seq_count += 1
        if(seq_count>1):
          # add the row for the sequence previously read
          for key in seq_dict.keys():
            if(seq_dict[key]=="N/A"):
              seq_dict[key] = np.nan
          row_list.append([seq_dict[key] for key in seq_dict.keys()])
        # reset the dictionary for the new sequence
        seq_dict = dict((key,np.nan) for key in dict_keys)
        seq_dict['seq_ID'] = line.split()[2]
      
      # Step 2
      elif(line.find("rearrangement summary")!=-1 and ParseType==0):
        ParseType = 2
      elif(ParseType==2):
        ParseType = 0
        v = line.split()
        if(chainType in ["heavy","HC"] and v[3]=="VH"):
          seq_dict['V_best_identity'] = v[0]
          seq_dict['D_best_identity'] = v[1]
          seq_dict['J_best_identity'] = v[2]
          seq_dict['stop_codon'] = v[4]
          seq_dict['V_J_frame'] = v[5]
          seq_dict['productive'] = v[6]
          seq_dict['strand'] = v[7]
        elif(chainType in ["kappa","KC"] and v[2]=="VK"):
          seq_dict['V_best_identity'] = v[0]
          seq_dict['D_best_identity'] = np.nan
          seq_dict['J_best_identity'] = v[1]
          seq_dict['stop_codon'] = v[3]
          seq_dict['V_J_frame'] = v[4]
          seq_dict['productive'] = v[5]
          seq_dict['strand'] = v[6]
        elif(chainType in ["lambda","LC"] and v[2]=="VL"):
          seq_dict['V_best_identity'] = v[0]
          seq_dict['D_best_identity'] = np.nan
          seq_dict['J_best_identity'] = v[1]
          seq_dict['stop_codon'] = v[3]
          seq_dict['V_J_frame'] = v[4]
          seq_dict['productive'] = v[5]
          seq_dict['strand'] = v[6]
      
      # Step 3
      elif(line.find("region sequence details")!=-1 and ParseType==0):
        ParseType = 3
      elif(ParseType==3):
        ParseType = 0
        seq_dict['CDR3_nt'] = line.split()[1]
      
      # Step 4
      elif(line.find("CDR3-IMGT")!=-1 and ParseType==0):
        seq_dict['V_CDR3len'] = float(line.split()[4])
      
      # Step 5
      elif(line.find("hits found")!=-1 and ParseType==0):
        ParseCount = (int)(line.split()[1])
        if(ParseCount>0):
          ParseType = 5
        else:
          ParseType = 0
        flagV = False
        flagD = False
        flagJ = False
      elif(ParseType==5 and ParseCount>0):
        ParseCount = ParseCount-1
        if(ParseCount==0):
          ParseType = 0
        v = line.split()
        if(v[0]=="V" and flagV==False):
          flagV = True
          seq_dict['V_best_frac_of_matches'] = float(v[3])
          seq_dict['V_best_align_length'] = float(v[4])
          if(seq_dict['V_CDR3len']==seq_dict['V_CDR3len']):
            seq_dict['V_best_align_length_beforeCDR3'] = seq_dict['V_best_align_length'] - seq_dict['V_CDR3len']
          else:
            seq_dict['V_best_align_length_beforeCDR3'] = np.nan
          seq_dict['V_best_N_of_mm'] = float(v[5])
          seq_dict['V_best_gap_opens'] = float(v[6])
          seq_dict['V_best_total_gap_length'] = float(v[7])
          seq_dict['V_best_align_start_seq'] = float(v[8])
          seq_dict['V_best_align_end_seq'] = float(v[9])
          seq_dict['V_best_align_start_gene'] = float(v[10])
          seq_dict['V_best_align_end_gene'] = float(v[11])
          seq_dict['V_best_igBlast_score'] = float(v[13])
          seq_dict['V_best_aligned_query'] = v[14]
          seq_dict['V_best_aligned_germline'] = v[15]
        if(v[0]=="D" and flagD==False):
          flagD = True
          seq_dict['D_best_frac_of_matches'] = float(v[3])
          seq_dict['D_best_align_length'] = float(v[4])
          seq_dict['D_best_N_of_mm'] = float(v[5])
          seq_dict['D_best_gap_opens'] = float(v[6])
          seq_dict['D_best_total_gap_length'] = float(v[7])
          seq_dict['D_best_align_start_seq'] = float(v[8])
          seq_dict['D_best_align_end_seq'] = float(v[9])
          seq_dict['D_best_align_start_gene'] = float(v[10])
          seq_dict['D_best_align_end_gene'] = float(v[11])
          seq_dict['D_best_igBlast_score'] = float(v[13])
          seq_dict['D_best_aligned_query'] = v[14]
          seq_dict['D_best_aligned_germline'] = v[15]
        if(v[0]=="J" and flagJ==False):
          flagJ = True
          seq_dict['J_best_frac_of_matches'] = float(v[3])
          seq_dict['J_best_align_length'] = float(v[4])
          seq_dict['J_best_N_of_mm'] = float(v[5])
          seq_dict['J_best_gap_opens'] = float(v[6])
          seq_dict['J_best_total_gap_length'] = float(v[7])
          seq_dict['J_best_align_start_seq'] = float(v[8])
          seq_dict['J_best_align_end_seq'] = float(v[9])
          seq_dict['J_best_align_start_gene'] = float(v[10])
          seq_dict['J_best_align_end_gene'] = float(v[11])
          seq_dict['J_best_igBlast_score'] = float(v[13])
          seq_dict['J_best_aligned_query'] = v[14]
          seq_dict['J_best_aligned_germline'] = v[15]
  # here I close the 'with'
  
  # add the row for the last sequence read
  for key in seq_dict.keys():
    if(seq_dict[key]=="N/A"):
      seq_dict[key] = np.nan
  #df = df.append(seq_dict, ignore_index=True)
  row_list.append([seq_dict[key] for key in seq_dict.keys()])
  
  for c,col in enumerate(seq_dict.keys()):
    df[col] = [row_list[i][c] for i in range(len(row_list))]
  
  row_list = []

  # Preliminar quality filtering (minimal requirements)
  df = df[df['V_best_identity']==df['V_best_identity']]
  if requireJ:
    df = df[df['J_best_identity']==df['J_best_identity']]
  
  # Cast to int
  # The 'Int64' type by pandas is able to store nan values
  for col in ['V_CDR3len', 'V_best_align_length', 'V_best_align_length_beforeCDR3', 'V_best_N_of_mm', \
              'V_best_gap_opens', 'V_best_total_gap_length', 'V_best_align_start_seq', 'V_best_align_end_seq', \
              'V_best_align_start_gene', 'V_best_align_end_gene', 'D_best_align_length', 'D_best_N_of_mm', \
              'D_best_gap_opens', 'D_best_total_gap_length', 'D_best_align_start_seq', 'D_best_align_end_seq', \
              'D_best_align_start_gene', 'D_best_align_end_gene', 'J_best_align_length', 'J_best_N_of_mm', \
              'J_best_gap_opens', 'J_best_total_gap_length', 'J_best_align_start_seq', 'J_best_align_end_seq', \
              'J_best_align_start_gene', 'J_best_align_end_gene']:
    df[col] = df[col].astype('Int64')
  
  # reindex to change the column ordering
  df = df.reindex(columns=dict_keys)
  
  # Extraction of other info
  # High-quality sequences may still have no annotated CDR3, likely because of missing anchors on V or J (or both)
  df['V_best_family'] = df['V_best_identity'].apply(lambda x: x.split(',')[0].split('*')[0] if type(x) is str else np.nan)
  df['D_best_family'] = df['D_best_identity'].apply(lambda x: x.split(',')[0].split('*')[0] if type(x) is str else np.nan)
  df['J_best_family'] = df['J_best_identity'].apply(lambda x: x.split(',')[0].split('*')[0] if type(x) is str else np.nan)
  df['CDR3_nt_len'] = df['CDR3_nt'].apply(lambda x: len(x) if type(x) is str else np.nan)
  df['CDR3_aa'] = df.apply(lambda row: nt2aa(row['CDR3_nt']) if (row['V_J_frame']=="In-frame" and type(row['CDR3_nt']) is str) else np.nan, axis=1)
  df['CDR3_aa_len'] = df['CDR3_aa'].apply(lambda x: len(x) if type(x) is str else np.nan)
  df['stop_codon_CDR3'] = df.apply(lambda row: row['CDR3_aa'].find('*') if type(row['CDR3_aa']) is str else np.nan, axis=1)
  df['stop_codon_CDR3'] = df['stop_codon_CDR3'].apply(lambda x: "No" if x==-1 else ("Yes" if (x==x and x>=0) else np.nan))

  # Hyper-indels analysis - for list of lists, then becoming a string 
  df['errors_list'] = df.apply(lambda row: pairwise_comparison(row['V_best_aligned_query'],row['V_best_aligned_germline'],3), axis=1)
  df['errors_list_split'] = df['errors_list'].apply(lambda x: x[1:-1].split(',') if len(x)>2 else [])
  df['N_del'] = df['errors_list_split'].apply(lambda x: len([y for y in x if y[1]=="d"]))
  df['N_ins'] = df['errors_list_split'].apply(lambda x: len([y for y in x if y[1]=="i"]))
  df['N_indels'] = df.apply(lambda row: row['N_del']+row['N_ins'], axis=1)
  df['Total_L_del'] = df['errors_list_split'].apply(lambda x: sum([len(y[1:-1].split('|')[3]) for y in x if y[1]=="d"]))
  df['Total_L_ins'] = df['errors_list_split'].apply(lambda x: sum([len(y[1:-1].split('|')[3]) for y in x if y[1]=="i"]))
  
  # Drop unnecessary data
  df = df.drop(['V_CDR3len', 'D_best_aligned_query', 'D_best_aligned_germline', 'errors_list_split'], axis=1)
  
  return df

def reconstruct_gapped_seqs(row):
  
  raw_seq = row['raw_seq_nt']
  #germline = V_dict[row['V_best_identity'].split(',')[0]]
  
  # V segment
  reconstructed_query = row['V_best_aligned_query']
  reconstructed_germline = row['V_best_aligned_germline']
  
  if(row['J_best_identity']==row['J_best_identity']):
    
    # CDR3 portion
    start = row['V_best_align_end_seq'] - 1 + 1                # igBlast indexes are 1-based;
                                                               # then, we need to move 1 position further
                                                               # w.r.t. the end of the aligned V region
    end = row['J_best_align_start_seq'] - 1 - 1 + 1            # igBlast indexes are 1-based;
                                                               # then, we need to move 1 position backward
                                                               # w.r.t. the beginning of the aligned J region;
                                                               # finally, we need to include the last position, too
    
    subseq = raw_seq[start:end]
    reconstructed_query += subseq
    reconstructed_germline += subseq
    
    # J segment
    reconstructed_query += row['J_best_aligned_query']
    reconstructed_germline += row['J_best_aligned_germline']
  
  return reconstructed_query + ';' + reconstructed_germline

def revert_seq(row):
  
  query = row['gapped_query']
  germline = row['gapped_germline']
  
  reverted_seq = ""
  
  for i in range(len(query)):
    if(query[i]!='-' and germline[i]!='-'):
      reverted_seq += query[i]
    elif(query[i]=='-' and germline[i]!='-'):
      reverted_seq += germline[i]
    #elif(query[i]!='-' and germline[i]=='-'):
      # do nothing
    elif(query[i]=='-' and germline[i]=='-'):
      print('Dashes on both gapped query and germline!')
    
  return reverted_seq

def write_anchored_seqs(row):

  V_best = row['V_best_identity'].split(',')[0]
  # if sequence has not been trimmed already
  # seq = row['seq_nt'][row['V_best_align_start_seq']-1:row['V_best_align_end_seq']]
  # if sequence has been already trimmed once
  seq = row['seq_nt'][0:row['V_best_align_end_seq']-row['V_best_align_start_seq']+1]
  # germ = V_dict[V_best][row['V_best_align_start_gene']-1:row['V_best_align_end_gene']]

  return str(row['seq_ID']) + ";" + seq.upper() + ";" + V_best + ";" + str(row['V_best_align_start_gene']-1) + ";" + str(row['V_best_align_end_gene'])
