#!/usr/bin/python3
"""
Author: Aldrin Yim
Version 0.1 (alpha)
Date: Mar 2019
Usage: python3 expression_visualization.py --help

"""
import matplotlib
matplotlib.use('Agg')
import sys, os, subprocess, argparse, pysam, pandas as pd, seaborn as sns, matplotlib.pyplot as plt, math, statistics, numpy as np
from collections import ChainMap, defaultdict
from itertools import combinations
import multiprocessing, tqdm, itertools, csv

current_wd = os.getcwd()
print('working directory: ',current_wd)

parser = argparse.ArgumentParser()
parser.add_argument("-g", "--go", help=" ESSENTIAL GO results from GOzilla")
args = parser.parse_args()

if args.go == None:
    parser.print_help()
    sys.exit()

initi=0
x = []
y = []

color_map = []
color_coors = defaultdict(list)

with open(args.go) as f:
    for line in f:
        readin = line.rstrip()
        data=readin.split('\t')
        if initi == 0:
            initi = 1
        else:
            description = data[1]
            pvalue = float(data[2])
            enrichment = float(data[4])
            color = data[6]
            mlog10pvalue = -math.log(pvalue,10)
            log2enrichment = math.log(enrichment,2)

            x.append(mlog10pvalue)
            #y.append(log2enrichment)
            y.append(enrichment)

            if color != "grey":
                tcx = color+'x'
                tcy = color+'y'
                color_map.append(color)
                color_coors[tcx].append(mlog10pvalue)
                color_coors[tcy].append(enrichment)
                #color_coors[tcy].append(log2enrichment)

print(x,y)

fig = plt.figure(figsize=(5,5))
ax = fig.add_subplot(111)
ax.axis([0,math.ceil(max(x))+1,0,max(y)+0.03])
#ax.axis([-15,15,-12,12])
dot_size = 100
plt.grid(b=True, which='major', color='k',alpha=0.3, linestyle='--')
#ax.grid(False,which='both',color='grey',alpha=0.3)
#ax.grid(False)
ax.axhline(y=0,alpha=0.8,color='k',linewidth=1)
ax.axvline(x=0,alpha=0.8,color='k',linewidth=1)
ax.set_facecolor('white')
plt.scatter(x,y, s=dot_size, alpha=0.3, c='grey')

for c in color_map:
    tcx = c+'x'
    tcy = c+'y'
    plt.scatter(color_coors[tcx],color_coors[tcy], s=100, alpha=0.9, c=c)
#, marker='c'

#plt.show()
plt.savefig('PNS-unique-v2.pdf')
