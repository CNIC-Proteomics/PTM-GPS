import argparse
import pandas as pd
import numpy as np
import os
import subprocess
import sys
import shutil
import tempfile

idx = pd.IndexSlice

def read_group_table(group_table_path):
    
    group_table = pd.read_csv(group_table_path, sep='\t')

    groups = {col: group_table[col].dropna().tolist() for col in group_table.columns}

    return groups

def preprocess_input(input_file, groups):
    
    integrations = [('Z_pgm2p', 'pgm'), ('Z_pgm2p_dNM', 'pgm'),
                    ('Z_p2qf', 'p'), ('Z_qf2q', 'qf'), ('Z_q2all', 'q')]

    input_dir = os.path.dirname(input_file)
    temp_dir = tempfile.mkdtemp(prefix="tmp-", dir=input_dir)
    for g, samples in groups.items():
        for intgr, lowL in integrations:
            g2=os.path.join(temp_dir, g)
            group_dir = os.path.join(g2, intgr)
            os.makedirs(group_dir, exist_ok=True)
            
            selected_data = df.loc[:, [(lowL, 'LEVEL')]].join(
                df.loc[:, pd.IndexSlice[intgr, :]].dropna(axis=0, how='all'),
                how='right'
            ).drop_duplicates().droplevel(0, axis=1).loc[:, ['LEVEL'] + samples]

            output_file_path = os.path.join(group_dir, 'input.tsv')
            selected_data.to_csv(output_file_path, sep='\t', index=False)
    return temp_dir

def main(input_file, group_table, output_file_final):
    output_dir = os.path.dirname(input_file)
    script_dir = os.path.dirname(os.path.abspath(__file__))
    r_script_file = os.path.join(script_dir, "LIMMA_script.R")
    
    groups = read_group_table(group_table)

    tmpfolder=preprocess_input(input_file, groups)
    
    tmpfolder2 = tmpfolder.replace("\\", "/")

    create_r_script(r_script_file, groups, tmpfolder2)

    call_r_script(r_script_file, script_dir)
    
    process_limma_outputs(groups, tmpfolder, output_dir, output_file_final)
    
    shutil.rmtree(tmpfolder)
    

    

def create_r_script(r_script_file, groups, output_dir):
    with open(r_script_file, 'w') as f:
        f.write("library(limma)\n")
        f.write(f"setwd('{output_dir}')\n\n")
        f.write("integrations <- c('Z_pgm2p', 'Z_pgm2p_dNM', 'Z_p2qf', 'Z_qf2q', 'Z_q2all')\n\n")
        f.write("intgr <- integrations[1]\n\n")
        f.write("groups <- c('" + "', '".join(groups.keys()) + "')\n\n")
        
        f.write("for (g in groups) {\n")
        f.write("  for (intgr in integrations) {\n")
        f.write(f"    df <- read.csv(paste0('{output_dir}', '/', g, '/', intgr, '/', 'input.tsv'), header=T, sep='\t')\n")
        f.write("    rownames(df) <- df$LEVEL\n")
        f.write("    df$LEVEL <- NULL\n\n")
        f.write("    samples <- colnames(df)\n")
        f.write("    design <- matrix(rep(1, length(samples)))\n")
        f.write("    rownames(design) <- samples\n")
        f.write("    colnames(design) <- g\n\n")
        f.write("    fit <- lmFit(df,design)\n")
        f.write("    fit <- eBayes(fit)\n\n")
        f.write("    pvalue <- as.data.frame(fit$p.value)\n")
        f.write("    meanRow <- as.data.frame(fit$Amean)\n\n")
        f.write("    output <- merge(meanRow, pvalue, 'row.names', all=T)\n")
        f.write("    colnames(output) <- c('LEVEL', 'Mean', 'pvalue')\n\n")
        f.write(f"    write.table(output, paste0('{output_dir}', '/', g, '/', intgr, '/', 'output.tsv'), sep='\t', row.names=F, col.names=T)\n")
        f.write("  }\n")
        f.write("}\n")

def call_r_script(r_script_file, script_dir):

    rscript_path = os.path.join(script_dir, "R-Portable", "App", "R-Portable", "bin", "x64", "Rscript.exe")
    
    subprocess.call([rscript_path, r_script_file])

def process_limma_outputs(groups, tmp_folder, output_dir, final_file):
    integrations = [('Z_pgm2p', 'pgm'), ('Z_pgm2p_dNM', 'pgm'), ('Z_p2qf', 'p'), ('Z_qf2q', 'qf'), ('Z_q2all', 'q')]
    intgr, lowL = integrations[0]

    df_out = df.copy()
    
    for g in groups:

        for intgr, lowL in integrations:
            g2=os.path.join(tmp_folder, g)


            limma_output = pd.read_csv(os.path.join(g2, intgr, 'output.tsv'), sep='\t', low_memory=False)
            limma_output.columns = pd.MultiIndex.from_tuples([(lowL, 'LEVEL'), (f'{intgr}_dX_{g}', 'dX'), (f'{intgr}_limma_{g}', 'pvalue')])
            limma_output[f'{intgr}_logLimma_{g}', 'LPS'] = -np.log10(limma_output[(f'{intgr}_limma_{g}', 'pvalue')])*np.sign(limma_output[(f'{intgr}_dX_{g}', 'dX')])

            df_out = pd.merge(
                df_out, limma_output,
                how='left',
                on=[(lowL, 'LEVEL')]
                )

    output_file = os.path.join(output_dir, final_file)
    df_out.to_csv(output_file, sep='\t', index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="One Sample LIMMA test. VMGS")
    parser.add_argument("input_file", help="Path to the input CSV file.")
    parser.add_argument("group_table", help="Path to the group table CSV file.")
    parser.add_argument("output_table", default="LIMMA1Sample_NM_Tabla_final.tsv", help="Output file")
    args = parser.parse_args()
    
    df = pd.read_csv(args.input_file, sep='\t', header=[0,1], low_memory=False)
    list(set([i for i,j in df.columns if j=='LEVEL']))
    
    main(args.input_file, args.group_table, args.output_table)
