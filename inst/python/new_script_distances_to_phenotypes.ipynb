{
 "cells": [
  {
   "cell_type": "code",
   "execution_count": 166,
   "metadata": {},
   "outputs": [],
   "source": [
    "import csv\n",
    "import os\n",
    "import ssl\n",
    "import requests, sys\n",
    "import pandas as pd\n",
    "from itertools import *\n",
    "from datetime import date"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 167,
   "metadata": {},
   "outputs": [],
   "source": [
    "analysis_type = 'PHENOTYPES'"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 168,
   "metadata": {},
   "outputs": [],
   "source": [
    "accepted_list_source=['protein','complex','proteinfamily', 'smallmolecule']\n",
    "accepted_list_target=['protein','complex','proteinfamily','phenotype', 'smallmolecule']"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 228,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "7658"
      ]
     },
     "execution_count": 228,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "### Veronica insert filter_list\n",
    "home_dir = os.path.expanduser('~') \n",
    "\n",
    "phospho_path = home_dir + '/' + 'phosphoproteomics.tsv'\n",
    "prot_path = home_dir + '/' + 'proteomics.tsv'\n",
    "#tr_path = home_dir + '/' + 'transcriptomics.tsv'\n",
    "\n",
    "phosphoproteomics = pd.read_table(phospho_path, sep = '\\t')\n",
    "proteomics = pd.read_table(prot_path, sep = '\\t')\n",
    "#transcript = pd.read_table(tr_path, sep = '\\t')\n",
    "\n",
    "phosphoproteomics_genes = phosphoproteomics['gene_name']\n",
    "proteomics_genes = proteomics['gene_name']\n",
    "#transcript_genes = transcript['Gene']\n",
    "\n",
    "# Concatenate the elements and take the unique ones\n",
    "#filter_list = pd.concat([phosphoproteomics_genes, proteomics_genes, transcript_genes]).unique()\n",
    "filter_list = pd.concat([phosphoproteomics_genes, proteomics_genes]).unique()\n",
    "#len(filter_list)\n",
    "len(filter_list)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 212,
   "metadata": {},
   "outputs": [],
   "source": [
    "# Download SIGNOR\n",
    "requestURL='https://signor.uniroma2.it/getData.php?'\n",
    "r=requests.get(requestURL)\n",
    "if not r.ok:\n",
    "    r.raise_for_status()\n",
    "    sys.exit()\n",
    "signor=r.text"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 213,
   "metadata": {},
   "outputs": [],
   "source": [
    "### create dict to convert SIGNOR EFFECTs into shapes\n",
    "diz_effect_clean = {}\n",
    "diz_effect_integer = {}\n",
    "diz_shape= {}\n",
    "\n",
    "for el in signor.split('\\n'):\n",
    "    col=el.split('\\t')\n",
    "    if col[0] != '':\n",
    "        effect = col[8]\n",
    "        if 'down-' in effect:\n",
    "            diz_effect_clean[effect]= 'inhibition'\n",
    "        elif 'up-' in effect:\n",
    "            diz_effect_clean[effect]= 'activation'\n",
    "        elif 'form' in effect:\n",
    "            diz_effect_clean[effect]= 'binding'\n",
    "        else:\n",
    "            diz_effect_clean[effect]= 'unknown'\n",
    "diz_effect_integer['inhibition']= -1\n",
    "diz_effect_integer['activation']= 1\n",
    "diz_effect_integer['binding']= 1\n",
    "diz_effect_integer['unknown']= 0\n",
    "\n",
    "diz_shape['inhibition']= '--|'\n",
    "diz_shape['activation']= '-->'\n",
    "diz_shape['binding']= '--[]'\n",
    "diz_shape['unknown']= '--?'"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 65,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "{'inhibition': '--|', 'activation': '-->', 'binding': '--[]', 'unknown': '--?'}"
      ]
     },
     "execution_count": 65,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "diz_shape"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "### extract SIGNOR interactions\n",
    "all_entities_uni = []\n",
    "diz_degree = {}\n",
    "diz_out_degree = {}\n",
    "diz_in_degree = {}\n",
    "diz_genename2uniprot = {}\n",
    "diz_uniprot2genename = {}\n",
    "diz_input = {}\n",
    "diz_output = {}\n",
    "#diz_distance = {}\n",
    "diz_effect = {}\n",
    "diz_signor_id = {}\n",
    "\n",
    "\n",
    "for el in signor.split('\\n'):\n",
    "    col=el.split('\\t')\n",
    "    if col[0] != '':\n",
    "        name1=col[0].upper().replace(' ','_')\n",
    "        name2=col[4].upper().replace(' ','_')\n",
    "        type1=col[1]\n",
    "        type2=col[5]\n",
    "        uni1=col[2]\n",
    "        uni2=col[6]\n",
    "        effect = diz_effect_clean[col[8]]\n",
    "        direct = col[22]\n",
    "        score= col[27]\n",
    "        all_entities_uni.append(uni1)\n",
    "        all_entities_uni.append(uni2)\n",
    "        if score != '' and name1 in filter_list and (name2 in filter_list or type2 == 'phenotype'): ### Veronica\n",
    "            #print(name1, name2)\n",
    "        #if score != '' and uni1 in filter_list and uni2 in filter_list:\n",
    "            #distance = 1-float(score)\n",
    "            signor_id= col[26]\n",
    "            diz_genename2uniprot[name1]=uni1\n",
    "            diz_genename2uniprot[name2]=uni2\n",
    "            diz_uniprot2genename[uni1]=name1\n",
    "            diz_uniprot2genename[uni2]=name2\n",
    "            if (direct == 't' or type2 == 'phenotype') and (effect != 'unknown') and (type1 in accepted_list_source and type2 in accepted_list_target) :\n",
    "                try:\n",
    "                    diz_input[uni2].append(uni1)\n",
    "                except:\n",
    "                    diz_input[uni2] = []\n",
    "                    diz_input[uni2].append(uni1)\n",
    "                try:\n",
    "                    diz_output[uni1].append(uni2)\n",
    "                except:\n",
    "                    diz_output[uni1] = []\n",
    "                    diz_output[uni1].append(uni2)\n",
    "                #diz_distance[uni1+'|'+uni2] = distance\n",
    "                diz_signor_id[uni1+'|'+uni2] = signor_id\n",
    "                try:\n",
    "                    diz_effect[uni1+'|'+uni2].append(effect)\n",
    "                except:\n",
    "                    diz_effect[uni1+'|'+uni2] = []\n",
    "                    diz_effect[uni1+'|'+uni2].append(effect)\n",
    "            \n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 216,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "(3705, 3781)"
      ]
     },
     "execution_count": 216,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "len(diz_output), len(diz_input)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 230,
   "metadata": {},
   "outputs": [],
   "source": [
    "for entity in set(all_entities_uni):\n",
    "    try:\n",
    "        out_degree = len(set(diz_output[entity]))\n",
    "    except:\n",
    "        out_degree = 0 \n",
    "    try:\n",
    "        in_degree = len(set(diz_input[entity]))\n",
    "    except: \n",
    "        in_degree = 0\n",
    "    degree = in_degree + out_degree\n",
    "    diz_degree [entity] = degree\n",
    "    diz_out_degree [entity] = out_degree\n",
    "    diz_in_degree [entity] = in_degree"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 231,
   "metadata": {},
   "outputs": [],
   "source": [
    "diz_distance = {}\n",
    "parameter_1 = 1\n",
    "parameter_2 = 0.001\n",
    "\n",
    "\n",
    "for el in signor.split('\\n'):\n",
    "    col=el.split('\\t')\n",
    "    if col[0] != '':\n",
    "        name1=col[0].upper().replace(' ','_')\n",
    "        name2=col[4].upper().replace(' ','_')\n",
    "        type1=col[1]\n",
    "        type2=col[5]\n",
    "        uni1=col[2]\n",
    "        uni2=col[6]\n",
    "        effect = diz_effect_clean[col[8]]\n",
    "        direct = col[22]\n",
    "        score= col[27]\n",
    "        if score != '' and name1 in filter_list and (name2 in filter_list or type2 == 'phenotype'): ### Veronica \n",
    "            distance = 1-float(score)\n",
    "            corrected_distance = parameter_1*distance + parameter_2* diz_in_degree[uni2] + parameter_2* diz_out_degree[uni1]\n",
    "            if (direct == 't' or type2 == 'phenotype') and (effect != 'unknown') and (type1 in accepted_list_source and type2 in accepted_list_target) :\n",
    "                diz_distance[uni1+'|'+uni2] = corrected_distance\n",
    "            #print(uni1, uni2, distance, corrected_distance)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 232,
   "metadata": {},
   "outputs": [],
   "source": [
    "output_list_pathway=[]\n",
    "\n",
    "for el in signor.split('\\n')[1:]:\n",
    "    col=el.split('\\t')\n",
    "    if col[0] != '':\n",
    "        name1=col[0].upper().replace(' ','_')\n",
    "        name2=col[4].upper().replace(' ','_')\n",
    "        type1=col[1]\n",
    "        type2=col[5]\n",
    "        uni1=col[2]\n",
    "        uni2=col[6]\n",
    "        if type1 == 'phenotype':\n",
    "            output_list_pathway.append(uni1)\n",
    "        if type2 == 'phenotype':\n",
    "            output_list_pathway.append(uni2)\n",
    "\n",
    "output_list_pathway= set(output_list_pathway)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 233,
   "metadata": {},
   "outputs": [],
   "source": [
    "diz_effect_boolean = {}\n",
    "diz_effect_final = {}\n",
    "diz_effect_shape = {}\n",
    "\n",
    "for pair in diz_effect:\n",
    "    if len(set(diz_effect[pair])) < 2:\n",
    "        diz_effect_final[pair] = diz_effect[pair][0]\n",
    "        diz_effect_boolean[pair] = diz_effect_integer[diz_effect[pair][0]]\n",
    "        diz_effect_shape [pair] = diz_shape [diz_effect[pair][0]]\n",
    "            #print (set(diz_effect[pair]))\n",
    "    elif len(set(diz_effect[pair])) > 2 and 'inhibition' not in diz_effect[pair] and ('activation' in diz_effect[pair] and 'binding' in diz_effect[pair]):\n",
    "        diz_effect_final[pair]= 'binding'\n",
    "        diz_effect_boolean[pair]= 1\n",
    "        diz_effect_shape [pair] = '--[]'\n",
    "    else:\n",
    "        diz_effect_final[pair]= 'unknown'\n",
    "        diz_effect_boolean[pair]= 0\n",
    "        diz_effect_shape [pair] = '--?'"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 234,
   "metadata": {},
   "outputs": [],
   "source": [
    "###prepare input protein\n",
    "input_list_clean=[]\n",
    "for el_gn in diz_output: ## use this for full analysis  - entire proteome in SIGNOR\n",
    "#for el_gn in input_list:  ## use this as test list \n",
    "    try:\n",
    "        el= diz_genename2uniprot[el_gn.upper().replace(' ','_')]\n",
    "        if el in diz_output:\n",
    "            input_list_clean.append(el)\n",
    "    except:\n",
    "        el=el_gn\n",
    "        if el in diz_output:\n",
    "            input_list_clean.append(el)        \n",
    "\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 235,
   "metadata": {},
   "outputs": [],
   "source": [
    "###prepare output\n",
    "output_list_clean=[]\n",
    "for el_gn in output_list_pathway:## use this for full analysis  - entire proteome in SIGNOR\n",
    "    \n",
    "#for el_gn in output_list: ## use this as test list \n",
    "    try:\n",
    "        el= diz_genename2uniprot[el_gn.upper().replace(' ','_')]\n",
    "        if el in diz_input:\n",
    "            output_list_clean.append(el)\n",
    "    except:\n",
    "        #print(el)\n",
    "        el=el_gn\n",
    "        \n",
    "        if el in diz_input:\n",
    "            output_list_clean.append(el)\n",
    "            \n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 236,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "35467"
      ]
     },
     "execution_count": 236,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "len(signor.split('\\n'))"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 237,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "1566"
      ]
     },
     "execution_count": 237,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "len(input_list_clean)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 238,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "127"
      ]
     },
     "execution_count": 238,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "len(output_list_clean)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 239,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "93672"
      ]
     },
     "execution_count": 239,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "##this chunk takes approx 2 mins, it retrieves all the paths linking input entities to output entities\n",
    "\n",
    "diz_path_list={}\n",
    "n=0\n",
    "for start_node in input_list_clean:\n",
    "    score=0\n",
    "    m=0\n",
    "    tot_neigh= []\n",
    "    path_list = []\n",
    "    if start_node in diz_output:\n",
    "        for degree_1_neigh in set(diz_output[start_node]):\n",
    "            if degree_1_neigh in output_list_clean:\n",
    "                end_node= degree_1_neigh\n",
    "                path_list.append(start_node+'|'+ end_node)\n",
    "                n=n+1\n",
    "                m=m+1\n",
    "            if degree_1_neigh in diz_output and degree_1_neigh!= start_node:\n",
    "                for degree_2_neigh in set(diz_output[degree_1_neigh]):\n",
    "                    if  degree_2_neigh in output_list_clean:\n",
    "                        end_node= degree_2_neigh\n",
    "                        path_list.append(start_node+'|'+degree_1_neigh+ '|'+ end_node)\n",
    "                        n=n+1\n",
    "                        m=m+1\n",
    "                    if degree_2_neigh in diz_output and degree_2_neigh!= start_node and degree_2_neigh!= degree_1_neigh:\n",
    "                        for degree_3_neigh in set(diz_output[degree_2_neigh]):\n",
    "                            if  degree_3_neigh in output_list_clean:\n",
    "                                end_node= degree_3_neigh\n",
    "                                path_list.append(start_node+'|'+degree_1_neigh+ '|'+ degree_2_neigh+ '|'+ end_node)\n",
    "                                n=n+1\n",
    "                                m=m+1\n",
    "                                    \n",
    "                            if degree_3_neigh in diz_output and degree_3_neigh!= start_node and degree_3_neigh!= degree_1_neigh and degree_3_neigh!= degree_2_neigh:\n",
    "                                for degree_4_neigh in set(diz_output[degree_3_neigh]):\n",
    "                                    if  degree_4_neigh in output_list_clean:\n",
    "                                        end_node= degree_4_neigh\n",
    "                                        path_list.append(start_node+'|'+degree_1_neigh+ '|'+ degree_2_neigh+ '|'+ degree_3_neigh+ '|'+ end_node)\n",
    "                                        n=n+1\n",
    "                                        m=m+1\n",
    "\n",
    "    #print(start_node, ' paths: ', m)\n",
    "    diz_path_list[start_node] = path_list\n",
    "    \n",
    "n\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 240,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "93672"
      ]
     },
     "execution_count": 240,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 241,
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Today's date: 20230803\n",
      "LastAnalysis_PHENOTYPES/results_output_20230803/\n"
     ]
    }
   ],
   "source": [
    "'''\n",
    "today = date.today()\n",
    "today = str(today).replace('-','')\n",
    "print(\"Today's date:\", today)\n",
    "path = 'LastAnalysis_'+ analysis_type +'/results_output_' + today+'/'\n",
    "path2 = 'LastAnalysis_'+ analysis_type +'/results_output_general/'\n",
    "print (path)\n",
    "# Check whether the specified path exists or not\n",
    "isExist = os.path.exists(path)\n",
    "if not isExist:\n",
    "\n",
    "   # Create a new directory because it does not exist\n",
    "   os.makedirs(path)\n",
    "   print(\"The new directory is created!\")\n",
    "isExist = os.path.exists(path2)\n",
    "if not isExist:\n",
    "\n",
    "   # Create a new directory because it does not exist\n",
    "   os.makedirs(path2)\n",
    "   print(\"The new directory is created!\")\n",
    "'''\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 244,
   "metadata": {},
   "outputs": [],
   "source": [
    "##this chunk takes approx mins, it parses info for all the paths linking input entities to output entities\n",
    "## it writes the results in a file: Global_result_final_table_minimized.txt\n",
    "\n",
    "file_output= open(home_dir+'/Global_result_final_table_minimized.txt','w')\n",
    "header= 'QueryNode\\tPath_String\\tEndNode\\trelations_path\\tPath_Score\\tPath_Length\\tFinal_Effect\\tEndPathways\\tEndNode_score\\n'\n",
    "file_output.write(header)\n",
    "#file_output2= open(path2+'Global_result_final_table_minimized.txt','w')\n",
    "#file_output2.write(header)\n",
    "\n",
    "for start_node in diz_path_list:\n",
    "    for path in diz_path_list[start_node]:\n",
    "        end_node = path.split('|')[-1]\n",
    "        p=0\n",
    "        distance_score=0\n",
    "        path_effect= 1\n",
    "        path_string= ''\n",
    "        relations_path = []\n",
    "        while p < len(path.split('|'))-1:\n",
    "            distance_score = distance_score + diz_distance[(path.split('|')[p]+'|' + path.split('|')[p+1])]\n",
    "            path_effect = path_effect * diz_effect_boolean[(path.split('|')[p]+'|' + path.split('|')[p+1])]\n",
    "            path_string = path_string + diz_uniprot2genename[path.split('|')[p]]+ diz_effect_shape [(path.split('|')[p]+'|' + path.split('|')[p+1])]\n",
    "            relations_path.append(diz_signor_id[(path.split('|')[p]+'|' + path.split('|')[p+1])])\n",
    "            p=p+1\n",
    "        path_string = path_string + diz_uniprot2genename[end_node]\n",
    "        path_lenght = path.count('|')\n",
    "        to_write=[diz_uniprot2genename[start_node], \n",
    "                  path_string, \n",
    "                  diz_uniprot2genename[end_node] , \n",
    "                  ';'.join(relations_path),\n",
    "                  str(round(distance_score,3) ), \n",
    "                  str(path_lenght),\n",
    "                  str(path_effect), \n",
    "                  diz_uniprot2genename[end_node] ,\n",
    "                  '-']\n",
    "        file_output.write('\\t'.join(to_write) + '\\n')\n",
    "        #file_output2.write('\\t'.join(to_write) + '\\n')\n",
    "file_output.close() \n",
    "#file_output2.close()  \n",
    "\n"
   ]
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3 (ipykernel)",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.9.13"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 4
}
