#!/usr/bin/python3
"""
Author: Aldrin Yim
Version 0.1 (alpha)
Date: Mar 2019
Usage: python3 circos-data-generation.py --help

"""
import sys, os, subprocess, argparse, pysam, pandas as pd, seaborn as sns, matplotlib.pyplot as plt, math, statistics, numpy as np
from collections import ChainMap, defaultdict
from itertools import combinations
import multiprocessing, tqdm, itertools
from genes import Genes

current_wd = os.getcwd()
print('working directory: ',current_wd)

parser = argparse.ArgumentParser()
parser.add_argument("-x", "--excel", help=" ESSENTIAL Input excel file for pairwise comparison")
args = parser.parse_args()

if args.excel == None:
    parser.print_help()
    sys.exit()

size = 1000
colors = ['group1','group2','group3','group4','group5','group6']
def get_length_for_kary(ddict,cond1):
    return len(ddict[cond1])*size-1

def get_common_length_for_kary(ddict,cond1):
    return len(ddict[cond1])*size

def compare_conditions(ddict,cond1,cond2):
    return 0

def get_highlight_coors(kary_ord, mtag_dict, target_gene, target_cond):
    mtag = mtag_dict[target_cond]
    start = str((kary_ord[target_cond].index(target_gene)+1)*size-size+1)
    end = str((kary_ord[target_cond].index(target_gene)+1)*size-1)
    return (mtag, start, end)

conditions_dict = defaultdict(list)
all_genes = []

with open(args.excel) as f:
    for line in f:
        readin = line.rstrip()
        data = readin.split('\t')
        conditions_dict[data[0]]=list(set(data[1:]))
        all_genes.append(list(set(data[1:])))
f.close()

mc_tag_dict = {}
i = 1
for k,v in conditions_dict.items():
    mc_tag_dict[k] = 'mc'+str(i)
    i+=1
print(mc_tag_dict)

unique_genes = list(set(list(itertools.chain.from_iterable(all_genes))))
#print(len(unique_genes))
#print(unique_genes)

gconnect = Genes(unique_genes,conditions_dict)

for k,v in conditions_dict.items():
    gconnect.put_connection(k,v)

gconnect.cal_connectivity()
gconnect.specify_conditions(['PNS-specific','CNS-specific'])
sorted_gene_list = gconnect.get_sorted_marginal()
apoe_related_conds = gconnect.get_gene_conditions('Spp1')
print(gconnect.get_gene_connectivity('Spp1'))
print(apoe_related_conds)
print('reporter = ', gconnect.is_gene_specified('St14') )

#Construct karyotype
#koutf = open('mac_karyotype.txt','w')
#for k,v in conditions_dict.items():
#    koutf.write('chr - '+mc_tag_dict[k]+' '+k+' '+'0 '+str(get_length_for_kary(conditions_dict,k))+' red\n')
#koutf.close()


#for k,v in mc_tag_dict.items():
#    if k == "CNS-specific":
#        sorted_genes = gconnect.get_cond_specific_sorted_marginal(k)
#        for gname,count in sorted_genes:
#            if count > 1:
#                print(gname,"|",count)

#Construct order of karyotype based on connectivity
karyotype_order = defaultdict(list)
for k,v in mc_tag_dict.items():
    sorted_genes = gconnect.get_cond_specific_sorted_marginal(k)
    for gname,count in sorted_genes:
        #if k != "PNS-specific" and k != "CNS-specific":
        if count >= 2:
                #print(count)
            karyotype_order[k].append(gname)
        #else:
        #    karyotype_order[k].append(gname)

print(karyotype_order['PNS-specific'])
print(karyotype_order['CNS-specific'])
#Construct karyotype
koutf = open('mac_karyotype.txt','w')
for k,v in conditions_dict.items():
    koutf.write('chr - '+mc_tag_dict[k]+' '+k+' '+'0 '+str(get_length_for_kary(karyotype_order,k))+' white\n')

band_start = 0
for k,v in karyotype_order.items():
    for i in range(0,len(v)):
        band_start = i*size
        band_end = band_start+size-1
        color = ''
        if i % 2 == 0:
            color = 'grey'
        else:
            color = 'white'
        koutf.write('band '+mc_tag_dict[k]+' band1 band1 '+str(band_start)+' '+str(band_end)+' '+color+'\n')

koutf.close()

all_connected_genes = gconnect.get_all_genes_with_connection()

ideohloutf = open('ideogram_highlight.txt','w')
linkoutf = open('connection_link.txt','w')
textoutf = open('genes_label.txt','w')
transparency = ['_a3','_a3','_a2','_a2','_a1','_a1','_a1']
tdict={}
for ind_gene in all_connected_genes:
    conditions = gconnect.get_gene_conditions(ind_gene)
    if len(conditions) >= 2:
        tuplist = []

        if ind_gene not in tdict.keys():
            tdict[ind_gene]=conditions

        for tcond in conditions:
            tuplist.append(get_highlight_coors(karyotype_order,mc_tag_dict,ind_gene,tcond))
        #exhaust combinations

        all_combos = list(combinations(tuplist,2))
        for each_iter in all_combos:
            (condA_abr, condA_start, condA_end) = each_iter[0]
            (condB_abr, condB_start, condB_end) = each_iter[1]
            tcolor = 'grey'
            treporter = gconnect.is_gene_specified(ind_gene)
            if treporter['PNS-specific'] == 1:
                tcolor = 'red_a3'
                #tcolor = 'red'+transparency[len(conditions)-2]
                linkoutf.write(condA_abr+' '+condA_start+' '+condA_end+' '+condB_abr+' '+condB_start+' '+condB_end+' color='+tcolor+'\n')
                #if len(conditions)>=2:
                    #print(ind_gene,' ',len(conditions))

            elif treporter['CNS-specific'] == 1:
                tcolor = 'blue_a3'
                #tcolor = 'blue'+transparency[len(conditions)-2]
                linkoutf.write(condA_abr+' '+condA_start+' '+condA_end+' '+condB_abr+' '+condB_start+' '+condB_end+' color='+tcolor+'\n')

            #linkoutf.write(condA_abr+' '+condA_start+' '+condA_end+' '+condB_abr+' '+condB_start+' '+condB_end+' color='+colors[len(conditions)-2]+'\n')
            #linkoutf.write(condA_abr+' '+condA_start+' '+condA_end+' '+condB_abr+' '+condB_start+' '+condB_end+' color='+tcolor+'\n')
            #print(condA_abr+' '+condA_start+' '+condA_end+' '+condB_abr+' '+condB_start+' '+condB_end+'\n')

gene_conditionsf = open('gene_conds.txt','w')
for i,v in tdict.items():
    gene_conditionsf.write(i+'\t')
    for j in v:
        gene_conditionsf.write(j+'\t')
    gene_conditionsf.write('\n')
gene_conditionsf.close()

for ind_gene in all_connected_genes:
    if gconnect.has_condition('PNS-specific',ind_gene):
        if gconnect.get_gene_connectivity(ind_gene) > 3:
            (mt,tstart,tend)=get_highlight_coors(karyotype_order,mc_tag_dict,ind_gene,'PNS-specific')
            textoutf.write(mt+' '+tstart+' '+tend+' '+ind_gene+' color='+colors[gconnect.get_gene_connectivity(ind_gene)-2]+'n\n')

for k,v in karyotype_order.items():
    for i in range(0,len(v)):
        counts = len(gconnect.get_gene_conditions(v[i]))
        band_start = i*size
        band_end = band_start+size-1
        if counts >= 2:
            ideohloutf.write(mc_tag_dict[k]+' '+str(band_start)+' '+str(band_end)+' fill_color='+colors[counts-2]+'\n')

textoutf.close()
linkoutf.close()
ideohloutf.close()
