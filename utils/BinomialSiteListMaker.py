# -*- coding: utf-8 -*-
"""
Created on Wed Nov  2 10:46:08 2022
"""


#
# Import Libraries
#

import argparse
import configparser
import logging
import os
import pandas as pd
import re
from scipy.stats import binom
import sys
from statsmodels.stats.multitest import multipletests


#
# Constants
#

PEAK = 'PEAK'


#
# Local Functions
#


def getBinom(wdf, col):
    '''
    '''
    p, d, a, m, x = col
    

    # Get aa freq
    afreq = pd.Series(list(zip(*[
        (k1, k2-j)
        for i, j in zip(wdf[p].tolist(), wdf[m].tolist())
        for k1, k2 in zip(list(i), range(len(i))) if abs(k2-j) <= abs(x)
    ]))[0]).value_counts()

    afreq = {
        i: j/afreq.sum()
        for i, j in afreq.to_dict().items()
    }

    # Get modification freq
    dfreq = wdf[d].value_counts()

    dfreq = {
        i: j/dfreq.sum()
        for i, j in dfreq.to_dict().items()
    }

    # Get binom table
    fdf = wdf.loc[:, [d, a]].groupby([d, a]).size(
    ).to_frame().reset_index().rename(columns={0: 'x'})

    # P1 --> p = p(a)*p(d) & n = Total pdm in ref table
    fdf['p1'] = [dfreq[i]*afreq[j] for i, j in zip(fdf[d], fdf[a])]
    fdf['n1'] = wdf.shape[0]

    # P2 --> p = p(a) & ni = Total pdm with i mod
    fdf['p2'] = [afreq[i] for i in fdf[a]]

    d_size = wdf[d].value_counts().to_frame().reset_index()#.rename(columns={'index': 'd', 'd': 'n2'})
    d_size.columns = [d, 'n2']

    fdf = pd.merge(
        fdf,
        d_size,
        how='left',
        on=d
    )

    # binom = P( Bi(n,p) >= x )
    fdf['binom1'] = [1-binom.cdf(i-1, j, k)
                     for i, j, k in zip(fdf['x'], fdf['n1'], fdf['p1'])]
    fdf['binom2'] = [1-binom.cdf(i-1, j, k)
                     for i, j, k in zip(fdf['x'], fdf['n2'], fdf['p2'])]

    return fdf


#
# Main
#

def main(args):
    '''
    Parameters
    ----------
    args : Dictionary
            - infile: String to input table or Pandas dataframe
            - outfile: String to output folder
            - peptidoform_column: String indicating name of the column with pdm
            - x: Integer indicating amino acid window size (left and right)
            - peakorph_column: String indicating column peak assignation (PEAK constants contains the value).
                If None, all pdm will be considered.
            - scanfreq_column: String indicating column with scanfreq

    Returns
    -------
    biS: Pandas dataframe with Binomial pvalues obtained using PSM level
    biP: Pandas dataframe with Binomial pvalues obtained using pdm level
    * If outfile is indicated, biS and biP will be saved
    '''

    # Set column names
    pdm, p, d, a, m, x = args['peptidoform_column'], args['peptide_column'], \
        args['modifcation_column'], args['modified_residue_column'], \
            args['modified_position_column'], args['x']

    # Read infile
    if type(args['infile']) == pd.DataFrame:
        df = args['infile']
    else:
        logging.info(f"Reading infile: {args['infile']}")
        df = pd.read_csv(args['infile'], sep='\t', low_memory=False)

    if args['scanfreq_column']:
        logging.info(f"Duplicating pdm based on {args['scanfreq_column']}")
        df = df.loc[df.index.repeat(
            df[args['scanfreq_column']])].reset_index(drop=True)

    # Build working df
    if d=='' or p=='' or a=='' or m=='':
        p, d, a, m = 'p', 'd', 'a', 'm'
        
        pdmList = df[pdm].tolist()
    
        if args['peakorph_column']:
            logging.info(f"Filtering NM based on {args['peakorph_column']}")
            pdmList = df.loc[df[args['peakorph_column']] == PEAK, pdm].tolist()
    
        # if not args['include_nm']:
            # logging.info("Excluding NM (pdm without [Mod])")
        
        pdmListNM = [i for i in pdmList if '[' not in i]
        unassigned = pd.Series([i.split('_')[1] for i in pdmListNM]).value_counts().to_frame()
        unassigned.columns = ['Unnasigned']
        
        pdmList = [i for i in pdmList if '[' in i]
    
        logging.info("Obtaining working dataframe")
        wdf = [
            (i, re.search(r'(.)\[([^]]+)\]', i))
            for i in pdmList
        ]
    
        wdf = [
            # (i, *i.split('_'), 'U', int(len(i)/2)) if j == None else
            (i, re.sub(r'\[[^]]+\]', '', i), j.groups()[1],
             j.groups()[0], i.index('[')-1)  # m index is 0-based
            for i, j in wdf
        ]
        
    
        wdf = pd.DataFrame(wdf, columns=[pdm, p, d, a, m])
        
    else:
        wdf = pd.DataFrame(df, columns=[pdm, p, d, a, m])
        unassigned = wdf[wdf.m.isna()][[d]].value_counts().to_frame()
        unassigned.columns = ['Unnasigned']
        wdf = wdf[~wdf.m.isna()]

    if wdf.shape[0] == 0:
        logging.error('No modified peptidoform was detected. Exiting program...')
        return None, None


    logging.info("Calculating binomial pvalues at PSM level")
    biS = getBinom(wdf, [p, d, a, m, x])

    logging.info("Calculating binomial pvalues at PDM level")
    biP = getBinom(wdf.drop_duplicates(), [p, d, a, m, x])

    logging.info("Merging binomial tables")
    biS.columns = [i if i in [d, a] else f'{i}-PSM' for i in biS.columns]
    biP.columns = [i if i in [d, a] else f'{i}-PDM' for i in biP.columns]

    bi = pd.merge(
        biP,
        biS,
        on=[d, a],
        how='outer'
    )

    # Add FDR
    for i in ['binom1-PDM', 'binom2-PDM', 'binom1-PSM', 'binom2-PSM']:
        bi[f'{i}-qvalue'] = multipletests(bi[i], method='fdr_bh')[1]

    # Pivot table
    binom = args['binom'] #'binom1-PSM'
    q_thr = float(args['q_thr']) # 0.01
    values_pivot = args['values_pivot'] #'x-PSM'

    biPivot = pd.pivot_table(bi[bi[binom]<q_thr], index=d, columns=a, values=values_pivot)
    
    if args['show_unassigned']:
        biPivot = pd.concat([biPivot, unassigned])
    
    biPivot['total'] = biPivot.sum(axis=1)
    biPivot = biPivot.sort_values('total', ascending=False)

    if args['outfile']:
        logging.info(f"Writing outfile: {args['outfile']}")
        #bi.to_csv(args['outfile'], sep='\t', index=False)
        # Write output
        with pd.ExcelWriter(args['outfile']) as writer:
            bi.to_excel(writer, sheet_name='Raw', index=False)
            biPivot.to_excel(writer, sheet_name=f'PIVOT-{binom}-FDR-{q_thr}-{values_pivot}')

    else:
        return bi, biPivot


if __name__ == '__main__':

    # Argument parsing
    parser = argparse.ArgumentParser(
        description='BinomialResMod',
        epilog='''
        Example:
            python BinomialResMod.py
        ''')

    parser.add_argument('-c', '--config', dest='config',
                        required=False, default=None, help='Path to config file')

    parser.add_argument('-i', '--infile', dest='infile',
                        required=False, help='Path to input file')  # required
    parser.add_argument('-o', '--outfile', dest='outfile',
                        required=False, default=None, help='Path to output file')
    parser.add_argument('-p', '--pepcol', dest='peptidoform_column', required=False,
                        help='Column name with peptidoform: AAAA[mod]AAAA')  # required
    parser.add_argument('-k', '--peakcol', dest='peakorph_column', required=False, default=None,
                        help='Column name indicating PEAK/ORPHAN. If none, all are used.')  # convenient
    parser.add_argument('-s', '--scancol', dest='scanfreq_column', required=False, default=None,
                        help='Column name with scan frequency. If none, infile is considered to be PSM.')
    parser.add_argument('-x', '--x', dest='x', required=False,
                        default=5, help='Window size used to calculate frequency aa.')
    # parser.add_argument('-nm', '--incnm', dest='include_nm', action='store_true', default=False, help='Include NM')

    args = parser.parse_args()

    if args.config:
        config = configparser.ConfigParser()
        config.read(args.config)
        params = dict(config.items('Params'))
        params['x'] = int(params['x'])
        params['show_unassigned'] = params['show_unassigned'].lower() == 'true'

    else:
        params = args.__dict__

    # Logging
    logging.basicConfig(level=logging.INFO,
                        format="BinomialResMod"+' - ' +
                        str(os.getpid()) +
                        ' - %(asctime)s - %(levelname)s - %(message)s',
                        datefmt='%m/%d/%Y %I:%M:%S %p',
                        handlers=[
                            logging.FileHandler(os.path.splitext(
                                params['infile'])[0]+'.log'),
                            logging.StreamHandler()
                        ]
                        )

    # Main
    logging.info('start script: ' +
                 "{0}".format(" ".join([x for x in sys.argv])))
    main(params)
    logging.info('End script')
