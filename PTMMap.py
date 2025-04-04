#!/usr/bin/env python
# coding: utf-8

# Module metadata variables
__author__ = "Cristina Amparo Devesa Arbiol"
__credits__ = ["Cristina Amparo Devesa Arbiol", "Jose Rodriguez", "Jesus Vazquez"]
__license__ = "Creative Commons Attribution-NonCommercial-NoDerivs 4.0 Unported License https://creativecommons.org/licenses/by-nc-nd/4.0/"
__version__ = "0.0.1"
__maintainer__ = "Jose Rodriguez"
__email__ = "cristinaamparo.devesa@cnic.es;jmrodriguezc@cnic.es"
__status__ = "Development"


# Imports
import pandas as pd
import numpy as np
import plotly
import plotly.express as px
import plotly.graph_objects as go
import argparse
import logging
import os
import sys
import yaml
import warnings
warnings.filterwarnings("ignore")



###################
# Local functions #
###################
def readInfile(infile, pgm):
    '''    
    Read input file to dataframe.
    '''
    df = pd.read_csv(infile,sep="\t", encoding='latin',header=[0,1], low_memory=False)
    df.columns = df.columns.map('_'.join)

    df = df.replace('#NUM!',np.nan)
    df = df.drop_duplicates(subset=[pgm])

    return df

def applyStructure(row,New_FDR,g,NM,a,first_b,New_LPS,LPS_pgm2p,LPS_pgm2p_NM,n,FDR_pgm,p,e,d,a_g_d,pgmFreq,dicc_FDRNM,pgmFreqThreshold):

    row[e] = float(str(row[e]).split(";")[0])
    row[n] = float(str(row[n]).split(";")[0])
    if row[p] in dicc_FDRNM[New_FDR].keys(): 
        row[New_FDR] = dicc_FDRNM[New_FDR][row[p]] 
    if np.isnan(row[New_FDR]) == True: 
        row[New_FDR] = 1
    if np.isnan(row[FDR_pgm]) ==True: 
        row[FDR_pgm] = 1

    if np.isnan(row[New_FDR]) ==True and np.isnan(row[FDR_pgm]) ==True: 
        row[New_FDR] = np.nan
    elif row[New_FDR]< row[FDR_pgm]: 
        row[New_FDR] = row[New_FDR]
    elif row[New_FDR]> row[FDR_pgm]: 
        row[New_FDR] = row[FDR_pgm]    
    if row[g] ==NM: 
        row[a] = "U"
        if str(row[g])!= "nan": 
            dm = "|" +str(row[g])
        else:
            dm = str(row[d])
            dm ="|" + dm[:dm.find(".")+3]
        row[n] = (float(row[first_b])+float(row[e]))/2 # NM will be represented at the average position between first_b and first_e
        row[New_LPS] = row[LPS_pgm2p] 
            
        if row[p] in dicc_FDRNM[New_FDR].keys(): 
            row[New_FDR] = dicc_FDRNM[New_FDR][row[p]]  
    else:  
        row[New_LPS] = row[LPS_pgm2p_NM]
        if str(row[g])!= "nan": 
            dm = "|" +str(row[g])
        else:
            dm = str(row[d])
            dm ="|" + dm[:dm.find(".")+3]
        if row[a] == "U":
            row[n] = (float(row[first_b])+float(row[e]))/2 # NM will be represented at the average position between first_b and first_e
            

    if row[pgmFreq]>= pgmFreqThreshold: 
        row[a_g_d]= row[a]+dm
    else: 
        row[a_g_d] = row[a]
    return row

def obtaindf (df,New_FDR, g,a,n,first_b,LPS_pgm2p,LPS_pgm2p_NM,FDR_NM,FDR_pgm,FDR_p2qf,FDR_qf2q,Missing_Cleavages,LPS_p2qf,LPS_qf2q,e,description, p,q,qf,pFreq,pgmFreq, qfFreq,d,NM,threshold_pgm2p,pgmFreqThreshold):

    df = df.loc[:,[g, a,n,first_b,LPS_pgm2p,LPS_pgm2p_NM, FDR_NM, FDR_pgm,FDR_p2qf, FDR_qf2q,Missing_Cleavages,LPS_p2qf, LPS_qf2q, e , description,p, q, qf, pFreq, qfFreq, pgmFreq,d]]

    df = df.rename(columns={New_FDR:"New_FDR",FDR_pgm:"FDR_pgm",g:"g", a: "a",n:"n",first_b : "first_b", LPS_pgm2p: "LPS_pgm2p",
                        LPS_pgm2p_NM: "LPS_pgm2p_NM", FDR_NM: "FDR_NM",FDR_p2qf : "FDR_p2qf", FDR_qf2q:"FDR_qf2q",
                      Missing_Cleavages: "Missing_Cleavages",LPS_p2qf:"LPS_p2qf", LPS_qf2q:"LPS_qf2q", e : "e", description: "description",
                       p :"p", q:"q", qf:"qf", pFreq:"pFreq", qfFreq: "qfFreq", pgmFreq: "pgmFreq",d :"d"})
    df = df.astype({'FDR_NM': 'float64', 'FDR_pgm': 'float64','LPS_pgm2p_NM': 'float64','LPS_pgm2p': 'float64','FDR_p2qf': 'float64',
                'FDR_qf2q': 'float64','LPS_p2qf': 'float64','LPS_qf2q': 'float64','pgmFreq': 'float64','pFreq': 'float64','qfFreq': 'float64',
                'Missing_Cleavages': 'float64', "first_b": 'float64'})
    
    df =df[df['first_b'].notnull()] # if an incomplete input it will discard using this parameter

    dfpeptide= df[df.g.eq(NM)]
    dfpeptide_FDR= dfpeptide[dfpeptide.FDR_NM.le(threshold_pgm2p)]
    dfpeptide_FDR = dfpeptide_FDR[['p','FDR_NM']].rename(columns={"FDR_NM":"New_FDR"})
    dicc_FDRNM = dfpeptide_FDR.set_index('p').to_dict()
    dicc_FDRNM["New_FDR"].keys()
 
    df["New_LPS"] = np.nan
    df["New_FDR"] = np.nan
    df["a_g_d"] = ""

    
    df_final = df.apply(lambda y: applyStructure(y,"New_FDR","g",NM,"a","first_b","New_LPS","LPS_pgm2p","LPS_pgm2p_NM","n","FDR_pgm","p","e","d","a_g_d","pgmFreq",dicc_FDRNM,pgmFreqThreshold), axis = 1)
    
    return df_final

def p2qfMaker (dfp,listproteins,threshold_p2qf):
    df_p2qf = pd.DataFrame(columns=["p","qf","q","pFreq",'Missing_Cleavages',"LPS_p2qf","position","description","FDR_p2qf"],dtype=float) # Dataframe 2 is created with the aim of 
    dfp_proteins=dfp [dfp ['q'].isin (listproteins)]
    dfp_proteins = dfp_proteins.drop_duplicates(subset=['p'])
    cont = 0
    for index, row in dfp_proteins.iterrows():
        int_b = int(row["first_b"])
        int_e = int(row["e"])

        for i in range(int_b,int_e+1):
            cont = cont+1
            df_p2qf.loc[cont,"p"] = row["p"]
            df_p2qf.loc[cont,"pFreq"] = row["pFreq"]
            df_p2qf.loc[cont,"q"] = row["q"]
            df_p2qf.loc[cont,"Missing_Cleavages"] = row["Missing_Cleavages"]
            df_p2qf.loc[cont,"LPS_p2qf"] = row["LPS_p2qf"]   
            df_p2qf.loc[cont,"description"] = row["description"]
            df_p2qf.loc[cont,"FDR_p2qf"] = row["FDR_p2qf"]
            df_p2qf.loc[cont,"position"] = i 
    df_p2qf_filtered = df_p2qf[df_p2qf.FDR_p2qf.le(threshold_p2qf)].reset_index()
    return df_p2qf,df_p2qf_filtered

def qf2qMaker (dfqf, listproteins,threshold_qf2q):
    dfqf_2 = pd.DataFrame(columns=["qf","q","qfFreq","LPS_qf2q","position","description","FDR_qf2q"],dtype=float)
    dfqf_proteins=dfqf [dfqf ['q'].isin (listproteins)]
    dfqf_proteins = dfqf_proteins.drop_duplicates(subset=['qf'])
    cont = 0
    for index, row in dfqf_proteins.iterrows():
        W = int(row["qf"].replace("-","_").split(";")[0].split(":")[1].split("_")[0])
        V = int(row["qf"].replace("-","_").split(";")[0].split(":")[1].split("_")[1])
        for i in range(W,V+1):
            cont = cont+1
            dfqf_2.loc[cont,"qfFreq"] = row["qfFreq"]
            dfqf_2.loc[cont,"q"] = row["q"]
            dfqf_2.loc[cont,"LPS_qf2q"] = row["LPS_qf2q"]  
            dfqf_2.loc[cont,"FDR_qf2q"] = row["FDR_qf2q"] 
            dfqf_2.loc[cont,"description"] = row["description"]
            dfqf_2.loc[cont,"position"] = i 
            dfqf_2.loc[cont,"qf"] = row["qf"]
        
    dfqf_2_filtered= dfqf_2[dfqf_2.FDR_qf2q.le(threshold_qf2q)].reset_index()
    return dfqf_2, dfqf_2_filtered

def TablesMaker (df_final,threshold_p2qf,NM,New_LPS,New_FDR,threshold_pgm2p,FDR_qf2q,threshold_qf2q):
    # p2qf Table
    df_final["Missing_Cleavages"] = df_final["Missing_Cleavages"].replace(0,"DT").replace(1,"DP").replace(2,"DP").replace(3,"DP").replace(4,"DP")
    df_final["first_b"] = df_final["first_b"].astype(int)
    dfp =df_final[df_final["LPS_p2qf"].notnull()] 
    dfp_filtered= dfp[dfp.FDR_p2qf.le(threshold_p2qf)].reset_index()
     
    # pgm2p Table
    dfpgm = df_final[df_final[New_LPS].notnull()]
    dfpgm.loc[dfpgm["g"] !=NM, "g"] = 'Mod'
    dfpgm_filtered= dfpgm[dfpgm.New_FDR.le(threshold_pgm2p)].reset_index()
    
    #qf2q Table
    dfqf = df_final[df_final["LPS_qf2q"].notnull()]
    dfqf_filtered= dfqf[dfqf.FDR_qf2q.le(threshold_qf2q)].reset_index()
    
    # Proteins thar passs the filters 
    listproteins = list(set((list(dfp_filtered["q"])+list(dfpgm_filtered["q"])+list(dfqf_filtered["q"]))))
    logging.info("- proteins that pass the threshold p2qf: " + str(len(set(list(dfp_filtered["q"])))))
    logging.info("- proteins that pass the threshold pgm2p: " + str(len(set(list(dfpgm_filtered["q"])))))
    logging.info("- proteins that pass the threshold qf2q: " + str(len(set(list(dfqf_filtered["q"])))))
    
    df_p2qf,df_p2qf_filtered = p2qfMaker(dfp, listproteins,threshold_p2qf)
    dfqf_2,dfqf_2_filtered= qf2qMaker(dfqf, listproteins,threshold_qf2q)
    
    return df_p2qf,df_p2qf_filtered, dfqf_2,dfqf_2_filtered,dfpgm,dfpgm_filtered,listproteins

def plotting (df_p2qf_filtered,dfpgm_filtered,dfqf_2_filtered,group_path,listproteins,font_size,grid,plot_width,plot_height):

    listafail = []
    dfpgm_filtered["n"].astype('int')
    c  =0
    for prot in listproteins: 
        c = c+1
        q = prot
        df1= df_p2qf_filtered[df_p2qf_filtered.q.eq(prot)]
        df1pgm= dfpgm_filtered[dfpgm_filtered.q.eq(prot)]
        df1pgm["pgmFreq"] = df1pgm["pgmFreq"].astype(float)

        dfqf_3= dfqf_2_filtered[dfqf_2_filtered.q.eq(prot)]
        dfqf_3["qfFreq"] = dfqf_3["qfFreq"].astype(float)
        list1 = list(df1pgm["n"])
        list1 = [int(x) for x in list1]
        dfqf_3["q"] = dfqf_3["q"].replace(prot,"qf2q")
        listw= list1+list(dfqf_3["position"])+list(df1["position"])
        listw.sort()
        df1pgm["a_g"] = df1pgm["a_g_d"] + "<br>" + "p= "+ df1pgm["p"]
        try: 
            fig1 = px.scatter(df1, x="position", y="LPS_p2qf",
                        size="pFreq", color="Missing_Cleavages",
                                 hover_name="Missing_Cleavages", size_max=8, opacity=1, title=list(df1["description"])[0],color_discrete_map={"DT": "lightgreen", "DP": "black"}, width=400, height=400)
            fig1.update_traces(
                marker=dict(symbol="square", line=dict(width=0, color="DarkSlateGrey")),
                selector=dict(mode="markers"),)
        except: 
            fig1= "false"

        try:
            fig2 = px.scatter(df1pgm, x="n", y="New_LPS",
                    size="pgmFreq", color="g",
                             hover_name="g", size_max=90,text="a_g", hover_data={"a_g": True},color_discrete_map={"NM": "orchid", 'Mod' :"red"}, width=400, height=400)
            fig2.for_each_trace(lambda t: t.update(textfont_color=t.marker.color))
            fig2.for_each_trace(lambda t: t.update(textfont_color='rgba(0,0,0,0)'))
        except: 
            fig2 ="false"

        try:
            fig3 = px.scatter(dfqf_3, x="position", y="LPS_qf2q", size="qfFreq",color = 'q',hover_name= 'q', opacity=1, size_max=9,color_discrete_map={"qf2q": "orange"}, width=400, height=400)
            fig3.update_traces(
                marker=dict( symbol="square", line=dict(width=0, color="orange")),
                selector=dict(mode="markers"),)
        except:
            fig3 = "false"


        if fig1 =="false" and fig2!="false" and fig3!= "false":
            fig = go.Figure(data = fig3.data  + fig2.data)
        elif fig1 =="false" and fig2=="false" and fig3!= "false":
            fig = go.Figure(data = fig3.data)
        elif fig1 =="false" and fig2!="false" and fig3== "false":
            fig = go.Figure(data = fig2.data)
        elif fig1 !="false" and fig2=="false" and fig3== "false":     
            fig = go.Figure(data = fig1.data)
        elif fig1 !="false" and fig2=="false" and fig3!= "false":
            fig = go.Figure(data = fig3.data  + fig1.data)  
        elif fig1!="false" and fig2!="false" and fig3== "false":
            fig = go.Figure(data = fig1.data  + fig2.data) 
        else: 
            fig = go.Figure(data = fig3.data + fig1.data + fig2.data)
        try: 
            fig.update_xaxes(range=[0,(listw[-1]+100)])
        except: 
            listafail.append(prot)
            pass

        try : 
            fig.update_layout(title_text=list(dfqf_3["description"])[0])
        except: 
            try : 
                fig.update_layout(title_text=list(df1["description"])[0])
            except: 
                fig.update_layout(title_text=list(df1pgm["description"])[0])
        
        fig.update_traces(textfont_size=1)
        fig.update_layout(yaxis = dict(tickfont = dict(size=font_size)))
        fig.update_layout(xaxis = dict(tickfont = dict(size=font_size)))
        fig.update_layout(paper_bgcolor = "rgba(0,0,0,0)",plot_bgcolor = "rgba(0,0,0,0)")
        if grid == "No": 
            fig.update_yaxes(showline=True, linewidth=5, linecolor='black', gridcolor='white', gridwidth=0,zeroline=True, zerolinewidth=5, zerolinecolor='black')
        else: 
            fig.update_yaxes(showline=True, linewidth=3, linecolor='black', gridcolor='black', gridwidth=1,zeroline=True, zerolinewidth=5, zerolinecolor='black',)
            fig.update_xaxes(showline=True, linewidth=3, linecolor='black', 
                  gridcolor='black', gridwidth=1)
            
        fig.update_xaxes(tickangle=-90, showline=True, linewidth=5, linecolor='black')
        fig.update_layout(width=plot_width, height=plot_height, xaxis=dict(ticks="outside",ticklen=15, tickwidth=5), yaxis=dict(ticks="outside",ticklen=15, tickwidth=5))
        fig.write_html(group_path + "/" + q + ".html", config={'toImageButtonOptions': {'format': 'svg', 'filename': q, 'height': plot_height, 'width': plot_width, 'scale': 1}})



###################
# Main functions #
###################
def main(config):

    """
    Reading configuration file
    """
    logging.info("Reading PTMap configuration file")

    # extract required column parameters
    pgm_first_header, pgm_second_header = config.get("pgm_column_name", ["", ""])
    pgm = f"{pgm_first_header}_{pgm_second_header}"
    g_first_header, g_second_header = config.get("g_column_name", ["", ""])
    g = f"{g_first_header}_{g_second_header}"
    a_first_header, a_second_header = config.get("a_column_name", ["", ""])
    a = f"{a_first_header}_{a_second_header}"
    n_first_header, n_second_header = config.get("n_column_name", ["", ""])
    n = f"{n_first_header}_{n_second_header}"
    e_first_header, e_second_header = config.get("e_column_name", ["", ""])
    e = f"{e_first_header}_{e_second_header}"
    p_first_header, p_second_header = config.get("p_column_name", ["", ""])
    p = f"{p_first_header}_{p_second_header}"
    q_first_header, q_second_header = config.get("q_column_name", ["", ""])
    q = f"{q_first_header}_{q_second_header}"
    d_first_header, d_second_header = config.get("d_column_name", ["", ""])
    d = f"{d_first_header}_{d_second_header}"
    qf_first_header, qf_second_header = config.get("qf_column_name", ["", ""])
    qf = f"{qf_first_header}_{qf_second_header}"
    pFreq_first_header, pFreq_second_header = config.get("pFreq_column_name", ["", ""])
    pFreq = f"{pFreq_first_header}_{pFreq_second_header}"
    qfFreq_first_header, qfFreq_second_header = config.get("qfFreq_column_name", ["", ""])
    qfFreq = f"{qfFreq_first_header}_{qfFreq_second_header}"
    pgmFreq_first_header, pgmFreq_second_header = config.get("pgmFreq_column_name", ["", ""])
    pgmFreq = f"{pgmFreq_first_header}_{pgmFreq_second_header}"
    first_b_first_header, first_b_second_header = config.get("first_b_column_name", ["", ""])
    first_b = f"{first_b_first_header}_{first_b_second_header}"
    description_first_header, description_second_header = config.get("description_column_name", ["", ""])
    description = f"{description_first_header}_{description_second_header}"
    Missing_Cleavages_first_header, Missing_Cleavages_second_header = config.get("Missing_Cleavages_column_name", ["", ""])
    Missing_Cleavages = f"{Missing_Cleavages_first_header}_{Missing_Cleavages_second_header}"

    NM = config.get("NM")

    # extract threshold parameters
    threshold_pgm2p_NM = config.get("threshold_pgm2p_NM") 
    threshold_pgm2p = config.get("threshold_pgm2p") 
    threshold_p2qf = config.get("threshold_p2qf") 
    threshold_qf2q = config.get("threshold_qf2q")     
    pgmFreqThreshold = config.get("pgmFreqThreshold")

    # extract map parameters
    font_size = config.get("font_size")
    grid= config.get("grid")
    plot_width= config.get("plot_width")
    plot_height= config.get("plot_height")

    # extract folders to save the maps
    path_plots_FDR = config.get("path_plots_with_threshold")
    path_plots = config.get("path_plots_Without_threshold")

    # extract the groups
    groups = config.get("groups")

    logging.info(f"Reading input file: {args.infile}...")
    df = readInfile(args.infile,pgm)

    logging.info(f'Processing by groups: {groups}')
    for grp in groups:
        logging.info(f"Preparing workspace for '{grp}'...")
        # prepare workspaces
        group_path_FDR = os.path.join(args.outdir, path_plots_FDR, grp)
        group_path = os.path.join(args.outdir, path_plots, grp)
        if not os.path.exists(group_path):
            os.makedirs(group_path, exist_ok=False)
        if not os.path.exists(group_path_FDR):
            os.makedirs(group_path_FDR, exist_ok=False)

        logging.info(f'- preparing parameters...')

        # read LPS column mappings
        LPS_mappings = config.get("LPS_ColumnNames", {})
        LPS_p2qf = f"{LPS_mappings['p2qc'][0]}_{grp}_{LPS_mappings['p2qc'][1]}"
        LPS_qf2q = f"{LPS_mappings['qc2q'][0]}_{grp}_{LPS_mappings['qc2q'][1]}"
        LPS_pgm2p = f"{LPS_mappings['pgm2p'][0]}_{grp}_{LPS_mappings['pgm2p'][1]}"
        LPS_pgm2p_NM = f"{LPS_mappings['pgm2p_NM'][0]}_{grp}_{LPS_mappings['pgm2p_NM'][1]}"

        # read NM column mappings
        NM_mappings = config.get("NM_ColumnNames", {})
        FDR_pgm = f"{NM_mappings['pgm2p'][0]}_{grp}_{NM_mappings['pgm2p'][1]}"
        FDR_NM = f"{NM_mappings['pgm2p_NM'][0]}_{grp}_{NM_mappings['pgm2p_NM'][1]}"

        # read Filter column mappings
        Filter_mappings = config.get("Filter_ColumnNames", {})
        FDR_p2qf = f"{Filter_mappings['p2qc'][0]}_{grp}_{Filter_mappings['p2qc'][1]}"
        FDR_qf2q = f"{Filter_mappings['qc2q'][0]}_{grp}_{Filter_mappings['qc2q'][1]}"

        logging.info("- obtaining group data...")
        df_final = obtaindf (df,"New_FDR",g,a,n,first_b,LPS_pgm2p,LPS_pgm2p_NM,FDR_NM,FDR_pgm,FDR_p2qf,FDR_qf2q,Missing_Cleavages,LPS_p2qf,LPS_qf2q,e,description, p,q,qf,pFreq,pgmFreq, qfFreq,d,NM,threshold_pgm2p,pgmFreqThreshold)

        logging.info("- preparing data...")
        df_p2qf,df_p2qf_filtered, dfqf_2,dfqf_2_filtered,dfpgm, dfpgm_filtered, listproteins= TablesMaker (df_final,threshold_p2qf,NM,"New_LPS","New_FDR",threshold_pgm2p,FDR_qf2q,threshold_qf2q)
        logging.info("- plotting filtered data...")
        plotting(df_p2qf_filtered,dfpgm_filtered,dfqf_2_filtered,group_path_FDR,listproteins,font_size, grid, plot_width,plot_height)
        logging.info("- plotting all data...")
        plotting(df_p2qf,dfpgm,dfqf_2,group_path,listproteins,font_size, grid, plot_width,plot_height)
    


if __name__ == '__main__':

    # parse arguments
    parser = argparse.ArgumentParser(
        description='PTMaps',
        epilog='''
        Example:
            python PTMaps.py
        ''')
      
    # default PTMaps configuration file
    defaultconfig = os.path.join(os.path.dirname(__file__), 'PTMMap.yml')
    parser.add_argument('-i', '--infile', required=True, help='Path to input file')
    parser.add_argument('-o', '--outdir', required=True, help='Output directory. Will be created if it does not exist')
    parser.add_argument('-c', '--config', default=defaultconfig, help='Path to custom config.ini file')
    parser.add_argument('-v', dest='verbose', action='store_true', help="Increase output verbosity")
    args = parser.parse_args()
   
   # parse config
    with open(args.config) as file:
        config = yaml.load(file, yaml.FullLoader)        
        # get the config section
        config = config.get('PTMMap')

    # prepare workspace
    outdir = args.outdir
    if outdir and not os.path.exists(outdir):
        os.makedirs(outdir, exist_ok=False)

    # logging debug level. By default, info level
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG,
                            format='%(asctime)s - %(levelname)s - %(message)s',
                            datefmt='%m/%d/%Y %I:%M:%S %p',
                            handlers=[logging.FileHandler(os.path.join(outdir, 'PTMMap.log')),
                                      logging.StreamHandler()])
    else:
        logging.basicConfig(level=logging.INFO,
                            format='%(asctime)s - %(levelname)s - %(message)s',
                            datefmt='%m/%d/%Y %I:%M:%S %p',
                            handlers=[logging.FileHandler(os.path.join(outdir, 'PTMMap_debug.log')),
                                      logging.StreamHandler()])

    # start main function
    logging.info('start script: '+"{0}".format(" ".join([x for x in sys.argv])))
    main(config)
    logging.info("end of the script")
