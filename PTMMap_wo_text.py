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
import configparser
import warnings
warnings.filterwarnings("ignore")



###################
# Local functions #
###################
def readInfile(infile):
    '''    
    Read input file to dataframe.
    '''
    df = pd.read_csv(infile,sep="\t", encoding='latin',header=[0,1])
    df.columns = df.columns.map('_'.join)
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
    if row[g] ==NM or row[a] == "U": 
        row[a] = "U"
        dm = ""
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


def obtaindf (infile,pgm,New_FDR, g,a,n,first_b,LPS_pgm2p,LPS_pgm2p_NM,FDR_NM,FDR_pgm,FDR_p2qf,FDR_qf2q,Missing_Cleavages,LPS_p2qf,LPS_qf2q,e,description, p,q,qf,pFreq,pgmFreq, qfFreq,d,NM,threshold_pgm2p,pgmFreqThreshold):
    #list_selectedproteins=["P24752","P50440"]
    df = readInfile(infile)
   

   
    df = df.replace('#NUM!',np.nan)
    df = df.drop_duplicates(subset=[pgm]) 
    
    

    df = df.loc[:,[g, a,n,first_b,LPS_pgm2p,LPS_pgm2p_NM, FDR_NM, FDR_pgm,FDR_p2qf, FDR_qf2q,Missing_Cleavages,LPS_p2qf, LPS_qf2q, e , description,p, q, qf, pFreq, qfFreq, pgmFreq,d]]

    df = df.rename(columns={New_FDR:"New_FDR",FDR_pgm:"FDR_pgm",g:"g", a: "a",n:"n",first_b : "first_b", LPS_pgm2p: "LPS_pgm2p",
                        LPS_pgm2p_NM: "LPS_pgm2p_NM", FDR_NM: "FDR_NM",FDR_p2qf : "FDR_p2qf", FDR_qf2q:"FDR_qf2q",
                      Missing_Cleavages: "Missing_Cleavages",LPS_p2qf:"LPS_p2qf", LPS_qf2q:"LPS_qf2q", e : "e", description: "description",
                       p :"p", q:"q", qf:"qf", pFreq:"pFreq", qfFreq: "qfFreq", pgmFreq: "pgmFreq",d :"d"})
    #df =df [df  ['q'].isin (list_selectedproteins)]

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



    logging.info("Proteins that pass the threshold p2qf: " + str(len(set(list(dfp_filtered["q"])))))
    logging.info("Proteins that pass the threshold pgm2p: " + str(len(set(list(dfpgm_filtered["q"])))))
    logging.info("Proteins that pass the threshold qf2q: " + str(len(set(list(dfqf_filtered["q"])))))
    
    df_p2qf,df_p2qf_filtered = p2qfMaker(dfp, listproteins,threshold_p2qf)
    dfqf_2,dfqf_2_filtered= qf2qMaker(dfqf, listproteins,threshold_qf2q)
    
    return df_p2qf,df_p2qf_filtered, dfqf_2,dfqf_2_filtered,dfpgm,dfpgm_filtered,listproteins


def plotting (df_p2qf_filtered,dfpgm_filtered,dfqf_2_filtered,group_path,listproteins):
    listafail = []
    dfpgm_filtered["n"].astype('int')
    c  =0
    ###listproteins=["K7GLT8","P80021", "Q29554", "A0A4X1VDD8", "Q95283", "F1RRT2"]
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
                             hover_name="g", size_max=90,text='a_g_d',color_discrete_map={"NM": "orchid", 'Mod' :"red"}, width=400, height=400)
            fig2.for_each_trace(lambda t: t.update(textfont_color=t.marker.color))
            
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
        fig.update_layout(yaxis = dict(tickfont = dict(size=50)))
        fig.update_layout(xaxis = dict(tickfont = dict(size=50)))
        fig.update_layout(paper_bgcolor = "rgba(0,0,0,0)",plot_bgcolor = "rgba(0,0,0,0)")
        fig.update_yaxes(showline=True, linewidth=5, linecolor='black', gridcolor='white', gridwidth=0,zeroline=True, zerolinewidth=5, zerolinecolor='black')
        fig.update_xaxes(tickangle=-90, showline=True, linewidth=5, linecolor='black', 
                  gridcolor='white', gridwidth=0)
        #fig.update_layout(width=1000, height=1000, xaxis=dict(ticks="outside",ticklen=15, tickwidth=5), yaxis=dict(ticks="outside",ticklen=15, tickwidth=5))
   

        fig.write_html(group_path+"/"+q+".html")




###################
# Main functions #
###################
def main(args):

    """
    Reading configuration file
    """
    config = configparser.ConfigParser(inline_comment_prefixes='#')
    config.read(args.config)
 
    logging.info("Reading PTMap configuration file") 

    

    pgm_first_header = config["PTMap_Parameters"].get("pgm_column_name") 
    g_first_header = config["PTMap_Parameters"].get("g_column_name") 
    a_first_header = config["PTMap_Parameters"].get("a_column_name") 
    n_first_header = config["PTMap_Parameters"].get("n_column_name") 
    e_first_header = config["PTMap_Parameters"].get("e_column_name") 
    p_first_header = config["PTMap_Parameters"].get("p_column_name") 
    q_first_header = config["PTMap_Parameters"].get("q_column_name") 
    d_first_header = config["PTMap_Parameters"].get("d_column_name")
    qf_first_header = config["PTMap_Parameters"].get("qf_column_name") 
    pFreq_first_header = config["PTMap_Parameters"].get("pFreq_column_name")
    qfFreq_first_header = config["PTMap_Parameters"].get("qfFreq_column_name")
    pgmFreq_first_header = config["PTMap_Parameters"].get("pgmFreq_column_name")
    first_b_first_header  = config["PTMap_Parameters"].get("first_b_column_name") 
    description_first_header = config["PTMap_Parameters"].get("description_column_name")
    Missing_Cleavages_first_header = config["PTMap_Parameters"].get("Missing_Cleavages_column_name")    

    pgm_second_header = config["PTMap_Parameters"].get("pgm_second_header_column_name") 
    g_second_header = config["PTMap_Parameters"].get("g_second_header_column_name") 
    a_second_header = config["PTMap_Parameters"].get("a_second_header_column_name") 
    n_second_header = config["PTMap_Parameters"].get("n_second_header_column_name") 
    e_second_header = config["PTMap_Parameters"].get("e_second_header_column_name") 
    p_second_header = config["PTMap_Parameters"].get("p_second_header_column_name") 
    q_second_header = config["PTMap_Parameters"].get("q_second_header_column_name") 
    d_second_header = config["PTMap_Parameters"].get("d_second_header_column_name")
    qf_second_header = config["PTMap_Parameters"].get("qf_second_header_column_name") 
    pFreq_second_header = config["PTMap_Parameters"].get("pFreq_second_header_column_name")
    qfFreq_second_header = config["PTMap_Parameters"].get("qfFreq_second_header_column_name")
    pgmFreq_second_header = config["PTMap_Parameters"].get("pgmFreq_second_header_column_name")
    first_b_second_header  = config["PTMap_Parameters"].get("first_b_second_header_column_name") 
    description_second_header = config["PTMap_Parameters"].get("description_second_header_column_name")
    Missing_Cleavages_second_header = config["PTMap_Parameters"].get("Missing_Cleavages_second_header_column_name")

    pgm =  pgm_first_header+"_"+pgm_second_header
    g =  g_first_header+"_"+g_second_header
    a =  a_first_header+"_"+a_second_header
    n =  n_first_header+"_"+n_second_header
    e =  e_first_header+"_"+e_second_header
    p =  p_first_header+"_"+p_second_header
    q =  q_first_header+"_"+q_second_header
    d = d_first_header+"_"+d_second_header
    qf = qf_first_header+"_"+qf_second_header
    pFreq = pFreq_first_header+"_"+pFreq_second_header
    qfFreq = qfFreq_first_header+"_"+qfFreq_second_header
    pgmFreq = pgmFreq_first_header+"_"+pgmFreq_second_header
    first_b  = first_b_first_header+"_"+first_b_second_header 
    description = description_first_header+"_"+description_second_header
    Missing_Cleavages=  Missing_Cleavages_first_header+"_"+Missing_Cleavages_second_header
    
    LPS_p2qf_first_header = config["PTMap_Parameters"].get("LPS_p2qf_column_name") 
    LPS_qf2q_first_header = config["PTMap_Parameters"].get("LPS_qf2q_column_name") 
    LPS_pgm2p_first_header = config["PTMap_Parameters"].get("LPS_pgm2p_column_name") 
    LPS_pgm2p_NM_first_header = config["PTMap_Parameters"].get("LPS_pgm2p_NM_column_name") 


    LPS_p2qf_second_header = config["PTMap_Parameters"].get("LPS_p2qf_second_header_column_name") 
    LPS_qf2q_second_header = config["PTMap_Parameters"].get("LPS_qf2q_second_header_column_name") 
    LPS_pgm2p_second_header = config["PTMap_Parameters"].get("LPS_pgm2p_second_header_column_name") 
    LPS_pgm2p_NM_second_header = config["PTMap_Parameters"].get("LPS_pgm2p_NM_second_header_column_name")

    LPS_p2qf = LPS_p2qf_first_header+"_"+LPS_p2qf_second_header
    LPS_qf2q = LPS_qf2q_first_header+"_"+LPS_qf2q_second_header
    LPS_pgm2p =  LPS_pgm2p_first_header+"_"+LPS_pgm2p_second_header
    LPS_pgm2p_NM = LPS_pgm2p_NM_first_header+"_"+LPS_pgm2p_NM_second_header
   

    FDR_pgm_first_header  = config["PTMap_Parameters"].get("Filter_pgm2p_NM_column_name") 
    FDR_NM_first_header  = config["PTMap_Parameters"].get("Filter_pgm2p_column_name") 
    FDR_p2qf_first_header = config["PTMap_Parameters"].get("Filter_p2qf_column_name") 
    FDR_qf2q_first_header = config["PTMap_Parameters"].get("Filter_qf2q_column_name")  

    FDR_pgm_second_header = config["PTMap_Parameters"].get("Filter_pgm2p_NM_second_header_column_name") 
    FDR_NM_second_header  = config["PTMap_Parameters"].get("Filter_pgm2p_second_header_column_name") 
    FDR_p2qf_second_header = config["PTMap_Parameters"].get("Filter_p2qf_second_header_column_name") 
    FDR_qf2q_second_header = config["PTMap_Parameters"].get("Filter_qf2q_second_header_column_name") 
    
    FDR_pgm =FDR_pgm_first_header+"_"+FDR_pgm_second_header
    FDR_NM  =FDR_NM_first_header+"_"+FDR_NM_second_header
    FDR_p2qf =FDR_p2qf_first_header+"_"+FDR_p2qf_second_header
    FDR_qf2q  =FDR_qf2q_first_header+"_"+FDR_qf2q_second_header


    
    threshold_pgm2p_NM = config["PTMap_Parameters"].getfloat("threshold_pgm2p_NM") 
    threshold_pgm2p = config["PTMap_Parameters"].getfloat("threshold_pgm2p") 
    threshold_p2qf = config["PTMap_Parameters"].getfloat("threshold_p2qf") 
    threshold_qf2q = config["PTMap_Parameters"].getfloat("threshold_qf2q") 
    
    pgmFreqThreshold = config["PTMap_Parameters"].getint("pgmFreqThreshold") 
    NM = config["PTMap_Parameters"].get("NM")
                                        
    path_plots_FDR = config["PTMap_Parameters"].get("path_plots_with_threshold")
    path_plots = config["PTMap_Parameters"].get("path_plots_Without_threshold")                                                                  
    group_path_FDR = os.path.join(args.output, path_plots_FDR)
    group_path = os.path.join(args.output, path_plots)
    if not os.path.exists(args.output):
        os.mkdir(args.output)

    if not os.path.exists(group_path):
        os.mkdir(group_path)  
    #else: 
        #os.remove(group_path)
        #os.mkdir(group_path) 
    if not os.path.exists(group_path_FDR):
        os.mkdir(group_path_FDR)
    #else: 
        #os.remove(group_path_FDR)
        #os.mkdir(group_path_FDR) 


    logging.info("Processing input file")
    df_final = obtaindf (args.infile,pgm,"New_FDR",g,a,n,first_b,LPS_pgm2p,LPS_pgm2p_NM,FDR_NM,FDR_pgm,FDR_p2qf,FDR_qf2q,Missing_Cleavages,LPS_p2qf,LPS_qf2q,e,description, p,q,qf,pFreq,pgmFreq, qfFreq,d,NM,threshold_pgm2p,pgmFreqThreshold)

    
    logging.info("Preparing data")
    df_p2qf,df_p2qf_filtered, dfqf_2,dfqf_2_filtered,dfpgm, dfpgm_filtered, listproteins= TablesMaker (df_final,threshold_p2qf,NM,"New_LPS","New_FDR",threshold_pgm2p,FDR_qf2q,threshold_qf2q)



    logging.info("Plotting")
    plotting(df_p2qf_filtered,dfpgm_filtered,dfqf_2_filtered,group_path_FDR,listproteins)
    plotting(df_p2qf,dfpgm,dfqf_2,group_path,listproteins)
    
    logging.info("end of the script")





if __name__ == '__main__':

    try:
        remove('Maps.ini')
    except:
        None

    # parse arguments
    parser = argparse.ArgumentParser(
        description='PTMaps',
        epilog='''
        Example:
            python PTMaps.py
        ''')
      
    # default PTMaps configuration file
    defaultconfig = os.path.join(os.path.dirname(__file__), "config/Solver.ini")
    parser.add_argument('-i', '--infile', required=True, help='Path to input file')
    parser.add_argument('-o',  '--output', required=True, help='Output directory. Will be created if it does not exist')
    parser.add_argument('-c', '--config', default=defaultconfig, help='Path to custom config.ini file')
    parser.add_argument('-v', dest='verbose', action='store_true', help="Increase output verbosity")
    args = parser.parse_args()
   
    # parse config
    config = configparser.ConfigParser(inline_comment_prefixes='#')
    config.read(args.config)      
    # logging debug level. By default, info level
    log_file = outfile = args.infile[:-4] + 'PTMas_log.txt'
    log_file_debug = outfile = args.infile[:-4] + 'PTMas_log_debug.txt'
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG,
                            format='%(asctime)s - %(levelname)s - %(message)s',
                            datefmt='%m/%d/%Y %I:%M:%S %p',
                            handlers=[logging.FileHandler(log_file_debug),
                                      logging.StreamHandler()])
    else:
        logging.basicConfig(level=logging.INFO,
                            format='%(asctime)s - %(levelname)s - %(message)s',
                            datefmt='%m/%d/%Y %I:%M:%S %p',
                            handlers=[logging.FileHandler(log_file),
                                      logging.StreamHandler()])

    #start main function
    logging.info('start script: '+"{0}".format(" ".join([x for x in sys.argv])))
    
    
    main(args)

