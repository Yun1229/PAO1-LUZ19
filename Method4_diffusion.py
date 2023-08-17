# -*- coding: utf-8 -*-
"""
Created on Mon Feb 20 15:48:37 2023

@author: jenny
"""

import pandas as pd

class Node:
    def __init__(self, name,targets,weight,levelWritten):
        self.name = name
        self.targets = []
        self.weight = weight
        self.levelWritten=0

class PaeNode(Node):
    def __init__(self, name,targets,weight,levelWritten):
        super().__init__(name,targets,weight,levelWritten)
        self.node_type = "pae"
        self.targets = targets
        self.weight = weight
        self.levelWritten =levelWritten
        
class CpdNode(Node):
    def __init__(self,name,targets,weight,changed,levelWritten):
        super().__init__(name,targets,weight,levelWritten)
        self.node_type = "cpd"
        self.targets = targets
        self.weight = weight
        self.changed = changed
        self.levelWritten =levelWritten
        

network = pd.read_excel("pae_network.xlsx")
len(network) #7256

network[network['target'].str.contains("0195.1")]
network.loc[1670]['target']='pae:PA01951'
network.loc[1674]['target']='pae:PA01951'
network[network['source'].str.contains("0195.1")]
network.loc[6541]['source']='pae:PA01951'
network.loc[6545]['source']='pae:PA01951'
network[network['source'].str.contains("1552.1")]
network.loc[4640]['source']='pae:PA15521'
network[network['source'].str.contains("1555.1")]
network.loc[4646]['source']='pae:PA15551'
network[network['target'].str.contains("1552.1")]
network.loc[4627]['target']='pae:PA15521'
network[network['target'].str.contains("1555.1")]
network.loc[4633]['target']='pae:PA15551'
sourcelist=[]
targetlist=[]
for i in range(len(network.index)):
    sourcelist.append(network.iloc[i]['source'].split(":")[1])
for i in range(len(network.index)):
    targetlist.append(network.iloc[i]['target'].split(":")[1])
renetworkdict = {
    'source':sourcelist,
    'target':targetlist,
    }

renetwork= pd.DataFrame(renetworkdict)
alllist=sourcelist+targetlist
alllist=list(set(alllist))
paelist=[]
cpdlist=[]



# %% import metabolomics 
Metabolomics = pd.read_excel("Metabolomics.xlsx", sheet_name="LUZ19")
M_5 = Metabolomics[((Metabolomics["FC time=2"]>=0.5) & (Metabolomics["AdjP time=2"]<0.05)) | ((Metabolomics["FC time=2"]<=-0.5) & (Metabolomics["AdjP time=2"]<0.05))]
M_5 = [cpd for cpd in M_5['Compound ID']]
M_5_string = []
for cpd in M_5:
    cpds = cpd.strip().split(';')
    for cpd in cpds:
        if cpd:
            M_5_string.append(cpd.strip())
len(M_5_string) #33


M_10 = Metabolomics[((Metabolomics["FC time=3"]>=0.5) & (Metabolomics["AdjP time=3"]<0.05)) | ((Metabolomics["FC time=3"]<=-0.5) & (Metabolomics["AdjP time=3"]<0.05))]
M_10 = [cpd for cpd in M_10['Compound ID']]
M_10_string = []
for cpd in M_10:
    cpds = cpd.strip().split(';')
    for cpd in cpds:
        if cpd:
            M_10_string.append(cpd.strip())
len(M_10_string) #76


M_15 = Metabolomics[((Metabolomics["FC time=4"]>=0.5) & (Metabolomics["AdjP time=4"]<0.05)) | ((Metabolomics["FC time=4"]<=-0.5) & (Metabolomics["AdjP time=4"]<0.05))]
M_15 = [cpd for cpd in M_15['Compound ID']]
M_15_string = []
for cpd in M_15:
    if not isinstance(cpd, float):
        cpds = cpd.strip().split(';')
        for cpd in cpds:
            if cpd:
                M_15_string.append(cpd.strip())
len(M_15_string) #105


# %% create Node objects 
#for M_5
for ele in alllist:
    if 'PA' in ele:
        exec(ele+ "=" + "PaeNode(" + "ele" + ",[],0,0)")
    if 'C' in ele:
        if ele not in M_5_string:
            exec(ele+ "=" + "CpdNode(" + "ele" + ",[],0,False,0)")
        else:
            exec(ele+ "=" + "CpdNode(" + "ele" + ",[],0,True,0)")
#for M_10            
for ele in alllist:
    if 'PA' in ele:
        exec(ele+ "=" + "PaeNode(" + "ele" + ",[],0,0)")
    if 'C' in ele:
        if ele not in M_10_string:
            exec(ele+ "=" + "CpdNode(" + "ele" + ",[],0,False,0)")
        else:
            exec(ele+ "=" + "CpdNode(" + "ele" + ",[],0,True,0)")
#for M_15            
for ele in alllist:
    if 'PA' in ele:
        exec(ele+ "=" + "PaeNode(" + "ele" + ",[],0,0)")
    if 'C' in ele:
        if ele not in M_15_string:
            exec(ele+ "=" + "CpdNode(" + "ele" + ",[],0,False,0)")
        else:
            exec(ele+ "=" + "CpdNode(" + "ele" + ",[],0,True,0)")            
            

            

# %% add targets to nodes            
for ele in alllist:
    if 'PA' in ele:
        exec("paelist.append("+ele+")")
        for ele2 in renetwork.loc[renetwork['source']==ele,'target'].values.tolist():
            exec(ele+ ".targets.append(" + ele2 + ")")        
    if 'C' in ele:
        exec("cpdlist.append("+ele+")")
        for ele2 in renetwork.loc[renetwork['source']==ele,'target'].values.tolist():
            exec(ele+ ".targets.append(" + ele2 + ")")
            
# len(cpdlist)+len(paelist)==len(alllist)
# 1231 + 1310 = 2541

# %% import transcriptomics
T_5 = pd.read_excel("LUZ19_LB_5min_vs_0min_DESeq.xlsx", sheet_name="padj < 0.01")
T_5_string = list(T_5["Gene"])
len(T_5_string) #168
T_5_string = [elem.replace('.', '') if '.' in elem else elem for elem in T_5_string]
T_5_objects = []
T_5_objects = [gene for gene in set(paelist) if gene.name in T_5_string]
len(T_5_objects) #54 

T_10 = pd.read_excel("LUZ19_LB_10min_vs_0min_DESeq.xlsx", sheet_name="padj < 0.01")
T_10_string = list(T_10["Gene"])
len(T_10_string) #283
T_10_string = [elem.replace('.', '') if '.' in elem else elem for elem in T_10_string]
T_10_objects = []
T_10_objects = [gene for gene in set(paelist) if gene.name in T_10_string]
len(T_10_objects) #85

T_15 = pd.read_excel("LUZ19_LB_15min_vs_0min_DESeq.xlsx", sheet_name="padj < 0.01")
T_15_string = list(T_15["Gene"])
len(T_15_string) #113
T_15_string = [elem.replace('.', '') if '.' in elem else elem for elem in T_15_string]
T_15_objects = []
T_15_objects = [gene for gene in set(paelist) if gene.name in T_15_string]
len(T_15_objects) #25

# %% create M_objects
M_5_objects = []
M_5_objects = [meta for meta in set(cpdlist) if meta.name in M_5_string]
len(M_5_string) #33
len(M_5_objects) #30 

M_10_objects = []
M_10_objects = [meta for meta in set(cpdlist) if meta.name in M_10_string]
len(M_10_string) #76
len(M_10_objects) #67

M_15_objects = []
M_15_objects = [meta for meta in set(cpdlist) if meta.name in M_15_string]
len(M_15_string) #105
len(M_15_objects) #92 

# %% assign weight to cpd 
  
# for cpd in cpdlist:
#     cpd.weight = 0


def findTargetInList(nodeList,level=0):    
    if level > 9 or not nodeList:
        return
    level+=1
    TargetList=[]
    for node in nodeList:
        if node not in linked_nodes:
            node.weight=1
            linked_nodes.append(node)  
        for target in node.targets:       
            value=node.weight/len(node.targets)
            if target.levelWritten==level or target.levelWritten==0:
                target.weight+=value
                target.levelWritten=level
                if target not in linked_nodes:
                    linked_nodes.append(target)
            if target not in TargetList:
                TargetList.append(target)
    findTargetInList(TargetList,level)



#for TP=5
T_5_impact_level = {}
impact_level=0
for gene in T_5_objects:
    T_5_impact_level[gene.name] = 0
    linked_nodes=[]
    for node in cpdlist:
        node.levelWritten = 0
        node.weight = 0      
    for node in paelist:
        node.levelWritten = 0
        node.weight = 0
    findTargetInList([gene])
    impact_level=0
    for node in linked_nodes:       
        if isinstance(node, CpdNode) and node.changed:
            impact_level += node.weight
            print([gene.name,node.name,node.weight,impact_level])              
    T_5_impact_level[gene.name] = impact_level
sorted_dict = {k: v for k, v in sorted(T_5_impact_level.items(), key=lambda item: item[1], reverse=True)}
df = pd.DataFrame(list(sorted_dict.items()), columns=["KEGG_ID", "impact_level"])
#df.to_excel("5min_diffusion.xlsx", index=False)
df_filtered = df[df['impact_level'] != 0]['KEGG_ID']
# df_filtered.to_csv("5min_diffusion_impact.tsv",sep='\t', index=False,header=False)
df_filtered = df_filtered.str.replace('PA', 'pae:PA')
df_filtered.to_excel("5min_diffusion_impact.xlsx", index=False) 




#for TP=10
T_10_impact_level = {}
impact_level=0
for gene in T_10_objects:
    T_10_impact_level[gene.name] = 0
    linked_nodes=[]
    for node in cpdlist:
        node.levelWritten = 0
        node.weight = 0      
    for node in paelist:
        node.levelWritten = 0
        node.weight = 0
    findTargetInList([gene])
    impact_level=0
    for node in linked_nodes:       
        if isinstance(node, CpdNode) and node.changed:
            impact_level += node.weight
            print([gene.name,node.name,node.weight,impact_level])               
    T_10_impact_level[gene.name] = impact_level
sorted_dict = {k: v for k, v in sorted(T_10_impact_level.items(), key=lambda item: item[1], reverse=True)}
df = pd.DataFrame(list(sorted_dict.items()), columns=["KEGG_ID", "impact_level"])
#df.to_excel("10min_diffusion.xlsx", index=False)
df_filtered = df[df['impact_level'] != 0]['KEGG_ID']
# df_filtered.to_csv("5min_diffusion_impact.tsv",sep='\t', index=False,header=False)
df_filtered = df_filtered.str.replace('PA', 'pae:PA')
df_filtered.to_excel("10min_diffusion_impact.xlsx", index=False) 




#for TP=15
T_15_impact_level = {}
impact_level=0
for gene in T_15_objects:
    T_15_impact_level[gene.name] = 0
    linked_nodes=[]
    for node in cpdlist:
        node.levelWritten = 0
        node.weight = 0      
    for node in paelist:
        node.levelWritten = 0
        node.weight = 0
    findTargetInList([gene])
    impact_level=0
    for node in linked_nodes:       
        if isinstance(node, CpdNode) and node.changed:
            impact_level += node.weight
            print([gene.name,node.name,node.weight,impact_level])               
    T_15_impact_level[gene.name] = impact_level
sorted_dict = {k: v for k, v in sorted(T_15_impact_level.items(), key=lambda item: item[1], reverse=True)}
df = pd.DataFrame(list(sorted_dict.items()), columns=["KEGG_ID", "impact_level"])
df.to_excel("15min_diffusion.xlsx", index=False)
df_filtered = df[df['impact_level'] != 0]['KEGG_ID']
# df_filtered.to_csv("5min_diffusion_impact.tsv",sep='\t', index=False,header=False)
df_filtered = df_filtered.str.replace('PA', 'pae:PA')
df_filtered.to_excel("15min_diffusion_impact.xlsx", index=False) 
