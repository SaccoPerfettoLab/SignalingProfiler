'''
import importlib

def check_and_install_package(package_name):
    try:
        # Try to import the package
        importlib.import_module(package_name)
        #print(f"{package_name} is already installed.")
    except ModuleNotFoundError:
        # If the package is not installed, try to install it
        install_package(package_name)

def install_package(package_name):
    import subprocess
    import sys

    #print(f"Installing {package_name}...")
    subprocess.check_call([sys.executable, "-m", "pip", "install", package_name])
    print(f"{package_name} has been successfully installed.")

# Example: Check and install required package
check_and_install_package('requests')
check_and_install_package('csv')
check_and_install_package('os')
check_and_install_package('ssl')
check_and_install_package('requests')
check_and_install_package('sys')
check_and_install_package('pandas')
check_and_install_package('itertools')
check_and_install_package('datetime')
'''

# After checking and installing, import the necessary libraries
import csv
import os
import ssl
import requests, sys
import pandas as pd
from itertools import *
from datetime import date

analysis_type = 'PHENOTYPES'

accepted_list_source=['protein','complex','proteinfamily', 'smallmolecule']
accepted_list_target=['protein','complex','proteinfamily','phenotype', 'smallmolecule']

### Veronica insert filter_list
home_dir = os.path.expanduser('~') 

phospho_path = home_dir + '/' + 'phosphoproteomics.tsv'
prot_path = home_dir + '/' + 'proteomics.tsv'
#tr_path = home_dir + '/' + 'transcriptomics.tsv'

phosphoproteomics = pd.read_table(phospho_path, sep = '\t')
proteomics = pd.read_table(prot_path, sep = '\t')
#transcript = pd.read_table(tr_path, sep = '\t')

phosphoproteomics_genes = phosphoproteomics['gene_name']
proteomics_genes = proteomics['gene_name']
#transcript_genes = transcript['Gene']

# Concatenate the elements and take the unique ones
#filter_list = pd.concat([phosphoproteomics_genes, proteomics_genes, transcript_genes]).unique()
filter_list = pd.concat([phosphoproteomics_genes, proteomics_genes]).unique()
#len(filter_list)
len(filter_list)

# Download SIGNOR
requestURL='https://signor.uniroma2.it/getData.php?'
r=requests.get(requestURL)
if not r.ok:
    r.raise_for_status()
    sys.exit()
signor=r.text

### create dict to convert SIGNOR EFFECTs into shapes
diz_effect_clean = {}
diz_effect_integer = {}
diz_shape= {}

for el in signor.split('\n'):
    col=el.split('\t')
    if col[0] != '':
        effect = col[8]
        if 'down-' in effect:
            diz_effect_clean[effect]= 'inhibition'
        elif 'up-' in effect:
            diz_effect_clean[effect]= 'activation'
        elif 'form' in effect:
            diz_effect_clean[effect]= 'binding'
        else:
            diz_effect_clean[effect]= 'unknown'
diz_effect_integer['inhibition']= -1
diz_effect_integer['activation']= 1
diz_effect_integer['binding']= 1
diz_effect_integer['unknown']= 0

diz_shape['inhibition']= '--|'
diz_shape['activation']= '-->'
diz_shape['binding']= '--[]'
diz_shape['unknown']= '--?'

### extract SIGNOR interactions
all_entities_uni = []
diz_degree = {}
diz_out_degree = {}
diz_in_degree = {}
diz_genename2uniprot = {}
diz_uniprot2genename = {}
diz_input = {}
diz_output = {}
#diz_distance = {}
diz_effect = {}
diz_signor_id = {}


for el in signor.split('\n'):
    col=el.split('\t')
    if col[0] != '':
        name1=col[0].upper().replace(' ','_')
        name2=col[4].upper().replace(' ','_')
        type1=col[1]
        type2=col[5]
        uni1=col[2]
        uni2=col[6]
        effect = diz_effect_clean[col[8]]
        direct = col[22]
        score= col[27]
        all_entities_uni.append(uni1)
        all_entities_uni.append(uni2)
        if score != '' and (name1 in filter_list or type1 != 'protein') and (name2 in filter_list or type2 == 'phenotype'): ### Veronica
            #print(name1, name2)
        #if score != '' and uni1 in filter_list and uni2 in filter_list:
            #distance = 1-float(score)
            signor_id= col[26]
            diz_genename2uniprot[name1]=uni1
            diz_genename2uniprot[name2]=uni2
            diz_uniprot2genename[uni1]=name1
            diz_uniprot2genename[uni2]=name2
            if (direct == 't' or type2 == 'phenotype') and (effect != 'unknown') and (type1 in accepted_list_source and type2 in accepted_list_target) :
                try:
                    diz_input[uni2].append(uni1)
                except:
                    diz_input[uni2] = []
                    diz_input[uni2].append(uni1)
                try:
                    diz_output[uni1].append(uni2)
                except:
                    diz_output[uni1] = []
                    diz_output[uni1].append(uni2)
                #diz_distance[uni1+'|'+uni2] = distance
                diz_signor_id[uni1+'|'+uni2] = signor_id
                try:
                    diz_effect[uni1+'|'+uni2].append(effect)
                except:
                    diz_effect[uni1+'|'+uni2] = []
                    diz_effect[uni1+'|'+uni2].append(effect)
            
for entity in set(all_entities_uni):
    try:
        out_degree = len(set(diz_output[entity]))
    except:
        out_degree = 0 
    try:
        in_degree = len(set(diz_input[entity]))
    except: 
        in_degree = 0
    degree = in_degree + out_degree
    diz_degree [entity] = degree
    diz_out_degree [entity] = out_degree
    diz_in_degree [entity] = in_degree
    
diz_distance = {}
parameter_1 = 1
parameter_2 = 0.001


for el in signor.split('\n'):
    col=el.split('\t')
    if col[0] != '':
        name1=col[0].upper().replace(' ','_')
        name2=col[4].upper().replace(' ','_')
        type1=col[1]
        type2=col[5]
        uni1=col[2]
        uni2=col[6]
        effect = diz_effect_clean[col[8]]
        direct = col[22]
        score= col[27]
        if score != '' and (name1 in filter_list or type1 != 'protein') and (name2 in filter_list or type2 == 'phenotype'): ### Veronica 
            distance = 1-float(score)
            corrected_distance = parameter_1*distance + parameter_2* diz_in_degree[uni2] + parameter_2* diz_out_degree[uni1]
            if (direct == 't' or type2 == 'phenotype') and (effect != 'unknown') and (type1 in accepted_list_source and type2 in accepted_list_target) :
                diz_distance[uni1+'|'+uni2] = corrected_distance
            #print(uni1, uni2, distance, corrected_distance)
            
            
output_list_pathway=[]

for el in signor.split('\n')[1:]:
    col=el.split('\t')
    if col[0] != '':
        name1=col[0].upper().replace(' ','_')
        name2=col[4].upper().replace(' ','_')
        type1=col[1]
        type2=col[5]
        uni1=col[2]
        uni2=col[6]
        if type1 == 'phenotype':
            output_list_pathway.append(uni1)
        if type2 == 'phenotype':
            output_list_pathway.append(uni2)

output_list_pathway= set(output_list_pathway)


diz_effect_boolean = {}
diz_effect_final = {}
diz_effect_shape = {}

for pair in diz_effect:
    if len(set(diz_effect[pair])) < 2:
        diz_effect_final[pair] = diz_effect[pair][0]
        diz_effect_boolean[pair] = diz_effect_integer[diz_effect[pair][0]]
        diz_effect_shape [pair] = diz_shape [diz_effect[pair][0]]
            #print (set(diz_effect[pair]))
    elif len(set(diz_effect[pair])) > 2 and 'inhibition' not in diz_effect[pair] and ('activation' in diz_effect[pair] and 'binding' in diz_effect[pair]):
        diz_effect_final[pair]= 'binding'
        diz_effect_boolean[pair]= 1
        diz_effect_shape [pair] = '--[]'
    else:
        diz_effect_final[pair]= 'unknown'
        diz_effect_boolean[pair]= 0
        diz_effect_shape [pair] = '--?'
        
###prepare input protein
input_list_clean=[]
for el_gn in diz_output: ## use this for full analysis  - entire proteome in SIGNOR
#for el_gn in input_list:  ## use this as test list 
    try:
        el= diz_genename2uniprot[el_gn.upper().replace(' ','_')]
        if el in diz_output:
            input_list_clean.append(el)
    except:
        el=el_gn
        if el in diz_output:
            input_list_clean.append(el)        

###prepare output
output_list_clean=[]
for el_gn in output_list_pathway:## use this for full analysis  - entire proteome in SIGNOR
    
#for el_gn in output_list: ## use this as test list 
    try:
        el= diz_genename2uniprot[el_gn.upper().replace(' ','_')]
        if el in diz_input:
            output_list_clean.append(el)
    except:
        #print(el)
        el=el_gn
        
        if el in diz_input:
            output_list_clean.append(el)
            
##this chunk takes approx 2 mins, it retrieves all the paths linking input entities to output entities

diz_path_list={}
n=0
for start_node in input_list_clean:
    score=0
    m=0
    tot_neigh= []
    path_list = []
    if start_node in diz_output:
        for degree_1_neigh in set(diz_output[start_node]):
            if degree_1_neigh in output_list_clean:
                end_node= degree_1_neigh
                path_list.append(start_node+'|'+ end_node)
                n=n+1
                m=m+1
            if degree_1_neigh in diz_output and degree_1_neigh!= start_node:
                for degree_2_neigh in set(diz_output[degree_1_neigh]):
                    if  degree_2_neigh in output_list_clean:
                        end_node= degree_2_neigh
                        path_list.append(start_node+'|'+degree_1_neigh+ '|'+ end_node)
                        n=n+1
                        m=m+1
                    if degree_2_neigh in diz_output and degree_2_neigh!= start_node and degree_2_neigh!= degree_1_neigh:
                        for degree_3_neigh in set(diz_output[degree_2_neigh]):
                            if  degree_3_neigh in output_list_clean:
                                end_node= degree_3_neigh
                                path_list.append(start_node+'|'+degree_1_neigh+ '|'+ degree_2_neigh+ '|'+ end_node)
                                n=n+1
                                m=m+1
                                    
                            if degree_3_neigh in diz_output and degree_3_neigh!= start_node and degree_3_neigh!= degree_1_neigh and degree_3_neigh!= degree_2_neigh:
                                for degree_4_neigh in set(diz_output[degree_3_neigh]):
                                    if  degree_4_neigh in output_list_clean:
                                        end_node= degree_4_neigh
                                        path_list.append(start_node+'|'+degree_1_neigh+ '|'+ degree_2_neigh+ '|'+ degree_3_neigh+ '|'+ end_node)
                                        n=n+1
                                        m=m+1

    #print(start_node, ' paths: ', m)
    diz_path_list[start_node] = path_list
    
##this chunk takes approx mins, it parses info for all the paths linking input entities to output entities
## it writes the results in a file: Global_result_final_table_minimized.txt

file_output= open(home_dir+'/Global_result_final_table_minimized.txt','w')
header= 'QueryNode\tPath_String\tEndNode\trelations_path\tPath_Score\tPath_Length\tFinal_Effect\tEndPathways\tEndNode_score\n'
file_output.write(header)
#file_output2= open(path2+'Global_result_final_table_minimized.txt','w')
#file_output2.write(header)

for start_node in diz_path_list:
    for path in diz_path_list[start_node]:
        end_node = path.split('|')[-1]
        p=0
        distance_score=0
        path_effect= 1
        path_string= ''
        relations_path = []
        while p < len(path.split('|'))-1:
            distance_score = distance_score + diz_distance[(path.split('|')[p]+'|' + path.split('|')[p+1])]
            path_effect = path_effect * diz_effect_boolean[(path.split('|')[p]+'|' + path.split('|')[p+1])]
            path_string = path_string + diz_uniprot2genename[path.split('|')[p]]+ diz_effect_shape [(path.split('|')[p]+'|' + path.split('|')[p+1])]
            relations_path.append(diz_signor_id[(path.split('|')[p]+'|' + path.split('|')[p+1])])
            p=p+1
        path_string = path_string + diz_uniprot2genename[end_node]
        path_lenght = path.count('|')
        to_write=[diz_uniprot2genename[start_node], 
                  path_string, 
                  diz_uniprot2genename[end_node] , 
                  ';'.join(relations_path),
                  str(round(distance_score,3) ), 
                  str(path_lenght),
                  str(path_effect), 
                  diz_uniprot2genename[end_node] ,
                  '-']
        file_output.write('\t'.join(to_write) + '\n')
        #file_output2.write('\t'.join(to_write) + '\n')
file_output.close() 
#file_output2.close()  
