# Module metadata variables
__author__ = ["Victor M. Guerrero-Sanchez", "Rafael Barrero"]
__credits__ = ["Rafael Barrero", "Jesus Vazquez"]
__license__ = "Apache License 2.0 https://www.apache.org/licenses/LICENSE-2.0"
__version__ = "1.0.0"
__maintainer__ = "Victor M. Guerrero-Sanchez"
__email__ = "vmguerrero@cnic.es"
__status__ = "Development"

import pandas as pd
import numpy as np
import argparse
import os
import yaml
import re
from Bio import SeqIO
from datetime import datetime

def get_args():
    parser = argparse.ArgumentParser(description="PTM analysis in hypermodified regions.")
    
    parser.add_argument("-f", "--fasta", type=str, required=True, help="Path to the FASTA file")
    parser.add_argument("-i", "--infile", type=str, required=True, help="Path to the TSV report file")
    parser.add_argument("-c", "--contrast", type=str, nargs='+', required=True, help="Contrast name(s) for analysis (can be a single value or a list)")
    parser.add_argument("-y", "--config", type=str, required=True, help="Path to the YAML configuration file")

    args = parser.parse_args()

    if not os.path.exists(args.fasta):
        raise FileNotFoundError(f"The FASTA file '{args.fasta}' does not exist.")
    if not os.path.exists(args.infile):
        raise FileNotFoundError(f"The report file '{args.infile}' does not exist.")
    if not os.path.exists(args.config):
        raise FileNotFoundError(f"The YAML configuration file '{args.config}' does not exist.")
    
    return args

def load_config(config_path):
    with open(config_path, 'r') as file:
        config = yaml.safe_load(file)
    return config

def create_log(output_dir, fasta_path, infile, SigMetric, MetricThr, contrast):
    log_file = os.path.join(output_dir, f"analysis_log_{contrast}.txt")
    with open(log_file, "w") as log:
        log.write(f"Analysis performed: {datetime.now()}\n")
        log.write(f"FASTA file used: {fasta_path}\n")
        log.write(f"Report file used: {infile}\n")
        log.write(f"Contrast used: {contrast}\n")
        log.write(f"Metric used: {SigMetric}\n")
        log.write(f"Cutoff value used: {MetricThr}\n")
    print(f"Log saved at: {log_file}")

def clean_column_names(df, contrast):
    df.columns = [re.sub(f"_{contrast}$", "", col) if isinstance(col, str) else col for col in df.columns]
    return df
print("######################################", flush=True)
print("Searching Hypermodified zones...", flush=True)
args = get_args()
fasta_path = args.fasta
infile = args.infile
contrasts = args.contrast
config_path = args.config

config = load_config(config_path)

SigMetric = config['SigMetric']
MetricThr = config['MetricThr']

print(f"Reading FASTA file: {fasta_path}", flush=True)
myfasta = {i.description: i.seq for i in SeqIO.parse(fasta_path, 'fasta')}

print(f"Reading report file: {infile}", flush=True)
rep = pd.read_csv(infile, sep='\t', header=[0, 1], low_memory=False)

infile_dir = os.path.abspath(os.path.dirname(args.infile))
output_dir = os.path.join(infile_dir, "HyperLocator_Analysis")
os.makedirs(output_dir, exist_ok=True)
print(f"Output directory created at: {output_dir}", flush=True)

all_results = []

for contrast in contrasts:
    print(f"Processing contrast: {contrast}", flush=True)

    print(f"Filtering data for contrast: {contrast}", flush=True)
    qv_NM = (f'Z_pgm2p_limma_NM_ONLY_{contrast}', SigMetric)
    pgm2p_dX = (f'Z_pgm2p_dX_{contrast}', 'dX')

    myPeptides = rep[np.logical_and.reduce([ 
        rep[qv_NM] < MetricThr, 
        rep[pgm2p_dX] < 0 
    ])][('p', 'LEVEL')].drop_duplicates().tolist()

    rep2 = rep[np.isin(rep[('p', 'LEVEL')], myPeptides)]
    rep3 = rep2[[ 
        ('pgmFreq', 'REL'),
        ('pgm', 'LEVEL'),
        ('p', 'LEVEL'),
        ('g', 'REL'),
        ('a', 'REL'),
        ('q', 'LEVEL'),
        ('description', 'REL'),
        ('m', 'REL'),
        ('b', 'REL'),
        ('e', 'REL'),
        ('n', 'REL'),
        (f'Z_pgm2p_dNM_dX_{contrast}', 'dX'),
        (f'Z_pgm2p_dNM_limma_{contrast}', SigMetric),
        (f'Z_pgm2p_dNM_limma_{contrast}', 'qvalue')
    ]].rename(columns={
        'm': 'pos_in_p',
        'b': 'q_begin',
        'e': 'q_end',
        'n': 'pos_in_q',
        (f'Z_pgm2p_dNM_limma_{contrast}', SigMetric): (f'pgm2p_dNM-{SigMetric}', '')
    }).drop_duplicates().copy()

    rep3.columns = pd.MultiIndex.from_tuples([(f'{j}(pgm2p_dNM)', '') if j in ['pvalue', 'qvalue'] else (i, j) for i, j in rep3.columns])

    print("Getting neighboring amino acids...", flush=True)
    aa = [
        [
            '' if ii < 0 or ii >= len(myfasta[k]) else myfasta[k][ii] 
            for ii in range(int(str(i).split(';')[0]) - 1 - 5, int(str(i).split(';')[0]) - 1 + 6)
        ] if not j == 'U' else 11 * ['']
        for i, j, k in zip(rep3[('pos_in_q', 'REL')].tolist(), rep3[('a', 'REL')].tolist(), rep3[('description', 'REL')].tolist())
    ]

    for n, i in enumerate(zip(*aa)):
        rep3[f'aa_{-5+n}'] = i

    rep4 = rep3[rep3[('g', 'REL')] != 'NM'].droplevel(1, axis=1)

    rep4['Comparison'] = contrast

    rep4 = clean_column_names(rep4, contrast)

    all_results.append(rep4)

    output_file = os.path.join(output_dir, f'PTMs_in_hypermodified_regions_{contrast}.xlsx')
    rep4.to_excel(output_file, index=False)

    print(f"Analysis for {contrast} completed. Results saved at: {output_file}", flush=True)

    create_log(output_dir, fasta_path, infile, SigMetric, MetricThr, contrast)

concated_results = pd.concat(all_results, ignore_index=True)

concated_output_file = os.path.join(output_dir, "Concated_HyperLocator_Results.xlsx")
concated_results.to_excel(concated_output_file, index=False)
print("Hyperlocator analysis finished.", flush=True)

print(f"Concatenated results saved at: {concated_output_file}", flush=True)
print("All contrasts have been processed.", flush=True)
