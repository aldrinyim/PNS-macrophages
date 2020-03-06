#!/usr/bin/python3
"""
Author: Aldrin Yim
Version 0.1 (alpha)
Date: Oct 2018
Usage: python3 edgeR_foldchange_plotting.py --help

"""

import sys, os, subprocess, argparse, glob, scipy, numpy as np, math, matplotlib.pyplot as plt, seaborn as sns
import multiprocessing, tqdm
import scipy as sp
from joblib import Parallel, delayed
from scipy.stats import stats
from collections import ChainMap, defaultdict

parser = argparse.ArgumentParser()
parser.add_argument("-c", "--classes", help=" ESSENTIAL Classification")
parser.add_argument("-s", "--subclass", help=" ESSENTIAL Subclassification")
parser.add_argument("-e", "--deg", help=" ESSENTIAL DEG matrix from edgeR")
parser.add_argument("-d", "--dir", help=" ESSENTIAL EdgeR diff exp output directory")
parser.add_argument("-o", "--out", help=" ESSENTIAL Output file name")
parser.add_argument("-m", "--matrix", help=" TEST TMM matrix or any normalized matrix")
args = parser.parse_args()

if args.subclass == None or args.out == None or args.classes == None:
    parser.print_help()
    sys.exit()

#file name = counts.matrix.Brain_vs_Fascia.edgeR.DE_results
prefix = 'counts.matrix.'
de_files = []
samples = defaultdict(list)
conditions = defaultdict(list)

def obtain_file(cond1,cond2):
    sample_in_cond1 = conditions[cond1]
    sample_in_cond2 = conditions[cond2]
    cond1_file = []
    cond2_file = []

    for file in de_files:
        for sam in sample_in_cond1:
            if sam in file:
                cond1_file.append(file)

    for file in cond1_file:
        for sam2 in sample_in_cond2:
            if sam2 in file:
                cond2_file.append(file)
    return cond2_file

with open(args.classes) as f:
    for line in f:
        readin = line.strip()
        data = readin.split('\t')
        samples[data[0]].append(data[1])
f.close()

print(samples)

with open(args.subclass) as f:
    for line in f:
        readin = line.strip()
        data = readin.split('\t')
        conditions[data[1]].append(data[0])
f.close()

print(conditions)

#Data structure of diff_exp_matrix
#diff_exp_matrix[gene][sample]=expression_value

initialize = 0
header_str = ''
header = []
diff_exp_matrix = defaultdict(dict)
gene_list = []
with open(args.matrix) as f:
    for line in f:
        readin = line.rstrip()
        data = readin.split('\t')
        if initialize == 0:
            initialize = 1
            header_str = readin
            header = data
        else:
            for i in range(1,len(data)):
                diff_exp_matrix[data[0]][header[i]]=float(data[i])
            gene_list.append(data[0])

for i in os.listdir(args.dir):
    if 'counts' in i:
        de_files.append(i)
print('number of files : ',len(de_files))

#annotated_genes = ("P2ry12","Cx3cr1","Tgfbr1","Trem2","Siglech","Ms4a4a","Apoe","Hivep2","Sall1","Tmem119","Gata6","Ccl2","Il1rl1","Cybb","Arg1","Gpnmb","Cp","Icam2","Spic","Cd74","Itgax","F11r","Ms4a7","Cd11a"\
#         ,"Pecam1","Cd93","Ccr2")

#annotated_genes = ("Spic","Vcam1","Pecam1","Itgal","Epcam","Clec7a","Ms4a6a")

annotated_genes = ("Cx3cr1","Sall1","Apoe","Pecam1","Siglech","Trem2","Arg1","Ccr2")

pns_can_defiles = obtain_file('PNS','Canonical')
print(len(pns_can_defiles))
print('--------------------------')
print(pns_can_defiles)
print('--------------------------')
pns_can_upreg = defaultdict(list)
can_pns_upreg = defaultdict(list)

for file in pns_can_defiles:
    initialize = 0
    with open(args.dir+file) as f:
        for line in f:
            readin = line.strip()
            data = readin.split('\t')
            if initialize == 0:
                initialize = 1
            else:
                gene_name = data[0]
                sam1 = data[1]
                sam2 = data[2]
                fc = float(data[3])
                pvalue = float(data[5])
                sam1_group = [key for key, value in conditions.items() if sam1 in value][0]
                #print(sam1,'|',sam1_group)

                if pvalue < 0.0001:
                    if sam1_group == 'PNS' and fc > 0:
                        pns_can_upreg[gene_name].append({sam1:[fc,pvalue]})
                    elif sam1_group == 'PNS' and fc < 0:
                        can_pns_upreg[gene_name].append({sam2:[-fc,pvalue]})
                    elif sam1_group == 'Canonical' and fc > 0:
                        can_pns_upreg[gene_name].append({sam1:[fc,pvalue]})
                    elif sam1_group == 'Canonical' and fc < 0:
                        pns_can_upreg[gene_name].append({sam2:[-fc,pvalue]})

nil = list(pns_can_upreg.keys())
for ni in nil:
    if can_pns_upreg.get(ni) != None:
        pns_can_upreg.pop(ni,None)
        can_pns_upreg.pop(ni,None)

cns_can_defiles = obtain_file('CNS','Canonical')
print(len(cns_can_defiles))
print('--------------------------')
print(cns_can_defiles)
print('--------------------------')
cns_can_upreg = defaultdict(list)
can_cns_upreg = defaultdict(list)

for file in cns_can_defiles:
    initialize = 0
    with open(args.dir+file) as f:
        for line in f:
            readin = line.strip()
            data = readin.split('\t')
            if initialize == 0:
                initialize = 1
            else:
                gene_name = data[0]
                sam1 = data[1]
                sam2 = data[2]
                fc = float(data[3])
                pvalue = float(data[5])
                sam1_group = [key for key, value in conditions.items() if sam1 in value][0]
                #print(sam1,'|',sam1_group)

                if pvalue < 0.0001:
                    if sam1_group == 'CNS' and fc > 0:
                        cns_can_upreg[gene_name].append({sam1:[fc,pvalue]})
                    elif sam1_group == 'CNS' and fc < 0:
                        can_cns_upreg[gene_name].append({sam2:[-fc,pvalue]})
                    elif sam1_group == 'Canonical' and fc > 0:
                        can_cns_upreg[gene_name].append({sam1:[fc,pvalue]})
                    elif sam1_group == 'Canonical' and fc < 0:
                        cns_can_upreg[gene_name].append({sam2:[-fc,pvalue]})

nil = list(cns_can_upreg.keys())
for ni in nil:
    if can_cns_upreg.get(ni) != None:
        cns_can_upreg.pop(ni,None)
        can_cns_upreg.pop(ni,None)

cns_pns_defiles = obtain_file('CNS','Canonical')
print(len(cns_pns_defiles))
print('--------------------------')
print(cns_pns_defiles)
print('--------------------------')
cns_pns_upreg = defaultdict(list)
pns_cns_upreg = defaultdict(list)

for file in cns_pns_defiles:
    initialize = 0
    with open(args.dir+file) as f:
        for line in f:
            readin = line.strip()
            data = readin.split('\t')
            if initialize == 0:
                initialize = 1
            else:
                gene_name = data[0]
                sam1 = data[1]
                sam2 = data[2]
                fc = float(data[3])
                pvalue = float(data[5])
                sam1_group = [key for key, value in conditions.items() if sam1 in value][0]
                #print(sam1,'|',sam1_group)
                if pvalue < 0.0001:
                    if sam1_group == 'CNS' and fc > 0:
                        cns_pns_upreg[gene_name].append({sam1:[fc,pvalue]})
                    elif sam1_group == 'CNS' and fc < 0:
                        pns_cns_upreg[gene_name].append({sam2:[-fc,pvalue]})
                    elif sam1_group == 'PNS' and fc > 0:
                        pns_cns_upreg[gene_name].append({sam1:[fc,pvalue]})
                    elif sam1_group == 'PNS' and fc < 0:
                        cns_pns_upreg[gene_name].append({sam2:[-fc,pvalue]})

cns_pns_shared_genes = defaultdict(list)
nil = list(cns_can_upreg.keys())
for ni in nil:
    if pns_can_upreg.get(ni) != None:
        cns_pns_shared_genes[ni].append(pns_can_upreg[ni])
        cns_can_upreg.pop(ni,None)
        pns_can_upreg.pop(ni,None)

cns_pns_codown_genes = defaultdict(list)
nil = list(can_cns_upreg.keys())
for ni in nil:
    if can_pns_upreg.get(ni) != None:
        cns_pns_codown_genes[ni].append(can_cns_upreg[ni])
        can_cns_upreg.pop(ni,None)
        can_pns_upreg.pop(ni,None)

annotated_order = []
annotated_genes_coors_x = []
annotated_genes_coors_y = []

pns_cns_shared_coors_x = []
pns_cns_shared_coors_y = []

pns_cns_codown_coors_x = []
pns_cns_codown_coors_y = []

pns_specific_coors_x = []
pns_specific_coors_y = []

can_specific_coors_x = []
can_specific_coors_y = []

cns_specific_coors_x = []
cns_specific_coors_y = []

all_gene_coors_x = []
all_gene_coors_y = []

fit_gene_coors_x = []
fit_gene_coors_y = []

remain_gene_coors_x = []
remain_gene_coors_y = []

low_gene_coors_x = []
low_gene_coors_y = []

# 6 files to be written
'''
cns_spcific_outf = open('cns_specific_genes.txt','w')
can_specific_outf = open('can_specific_genes.txt','w')
pns_specific_outf = open('pns_specific_genes.txt','w')
drg_specific_outf = open('drg_specific_genes.txt','w')

nerves_cns_outf = open('nerves_cns_shared_genes.txt','w')
nerves_drg_outf = open('nerves_drg_shared_genes.txt','w')
'''

def write_expression_matrix(headline,dematrix,genelist,outname):
    outf = open(outname,'w')
    outf.write(headline+'\n')
    header = headline.split('\t')
    for g in genelist:
        if g in dematrix:
            outf.write(g+'\t')
            for i in range(1,len(header)):
                outf.write(str(np.log2(dematrix[g][header[i]]))+'\t')
            outf.write('\n')
    outf.close()

def append_exp(gene,group):
    exp_list = []
    exp_count = 0
    for i in group:
        if diff_exp_matrix[gene][i] >= 1:
            if diff_exp_matrix[gene][i] > 5:
                exp_count+=1
            exp_list.append(diff_exp_matrix[gene][i])

    expressivity = exp_count/len(group)

    return (exp_list,expressivity)

def plot_ci_manual(t, s_err, n, xs, x2, y2, ax=None):
    """Return an axes of confidence bands using a simple approach.

    Notes
    -----
    .. math:: \left| \: \hat{\mu}_{y|x0} - \mu_{y|x0} \: \right| \; \leq \; T_{n-2}^{.975} \; \hat{\sigma} \; \sqrt{\frac{1}{n}+\frac{(x_0-\bar{x})^2}{\sum_{i=1}^n{(x_i-\bar{x})^2}}}
    .. math:: \hat{\sigma} = \sqrt{\sum_{i=1}^n{\frac{(y_i-\hat{y})^2}{n-2}}}

    References
    ----------
    .. [1]: M. Duarte.  "Curve fitting," JUpyter Notebook.
       http://nbviewer.ipython.org/github/demotu/BMC/blob/master/notebooks/CurveFitting.ipynb

    """
    #if ax is None:
    #    ax = plt.gca()

    ci = t*s_err*np.sqrt(1/n + (x2-np.mean(xs))**2/np.sum((x-np.mean(xs))**2))
    print('ci value : ', ci)
    ax.fill_between(x2, y2+ci, y2-ci, alpha=0.3, color='blue', edgecolor="")

    return ax

#Testing
"""
write_expression_matrix(header_str,diff_exp_matrix,cns_can_upreg,'plot3_cns_specific_genes.txt')
write_expression_matrix(header_str,diff_exp_matrix,pns_can_upreg,'plot3_pns_specific_genes.txt')
write_expression_matrix(header_str,diff_exp_matrix,cns_pns_shared_genes,'plot3_shared_cns_pns_genes.txt')


write_expression_matrix(header_str,diff_exp_matrix,cns_upreg,'plot2_cns_specific_genes.txt')
write_expression_matrix(header_str,diff_exp_matrix,nerves_upreg,'plot2_nerves_specific_genes.txt')
write_expression_matrix(header_str,diff_exp_matrix,can_upreg,'plot2_can_specific_genes.txt')
write_expression_matrix(header_str,diff_exp_matrix,drg_upreg,'plot2_drg_specific_genes.txt')
write_expression_matrix(header_str,diff_exp_matrix,nerve_cns_common_genes,'plot2_nerves_cns_shared_genes.txt')
write_expression_matrix(header_str,diff_exp_matrix,nerve_drg_common_genes,'plot2_nerves_drg_shared_genes.txt')
"""

for g in gene_list:
    nerves_exp = []
    can_exp = []
    cns_exp = []
    drg_exp = []

    pns_group = [y for x in [samples[j] for j in conditions['PNS']] for y in x]
    can_group = [y for x in [samples[j] for j in conditions['Canonical']] for y in x]
    cns_group = [y for x in [samples[j] for j in conditions['CNS']] for y in x]

    pns_exp,pns_expressivity = append_exp(g,pns_group)
    can_exp,can_expressivity = append_exp(g,can_group)
    cns_exp,cns_expressivity = append_exp(g,cns_group)

    pns = np.mean(pns_exp)
    can = np.mean(can_exp)
    cns = np.mean(cns_exp)

    if pns_expressivity < 0.9:
        if g in pns_can_upreg.keys():
            pns_can_upreg.pop(g,None)
        if g in cns_pns_shared_genes.keys():
            cns_can_upreg[g].append(cns_pns_shared_genes[g])
            cns_pns_shared_genes.pop(g,None)
    if pns/cns < 2 and pns/can < 2:
        pns_can_upreg.pop(g,None)

    if cns_expressivity < 0.9:
        if g in cns_can_upreg.keys():
            cns_can_upreg.pop(g,None)
        if g in cns_pns_shared_genes.keys():
            pns_can_upreg[g].append(cns_pns_shared_genes[g])
            cns_pns_shared_genes.pop(g,None)
    if cns/pns < 2 or cns/can < 2:
        cns_can_upreg.pop(g,None)

    if cns < 5:
        if g in cns_can_upreg.keys():
            cns_can_upreg.pop(g,None)
        if g in cns_pns_shared_genes.keys():
            pns_can_upreg[g].append(cns_pns_shared_genes[g])
            cns_pns_shared_genes.pop(g,None)

    if pns < 5:
        if g in pns_can_upreg.keys():
            pns_can_upreg.pop(g,None)
        if g in cns_pns_shared_genes.keys():
            cns_can_upreg[g].append(cns_pns_shared_genes[g])
            cns_pns_shared_genes.pop(g,None)

    """
    for i in nerves_group:
        nerves_exp.append(diff_exp_matrix[g][i])
    for i in can_group:
        can_exp.append(diff_exp_matrix[g][i])
    for i in cns_group:
        cns_exp.append(diff_exp_matrix[g][i])
    for i in drg_group:
        drg_exp.append(diff_exp_matrix[g][i])
    """


    #print(nerves,'|',can,'|',cns,'|',drg)

    #TWO comparisons will be performed
    if pns > 0 and can > 0 and cns > 0:
        #print(g)
        pns_can_ratio = np.log2(pns/can)
        cns_can_ratio = np.log2(cns/can)

        #pns_can_ratio = pns/can
        #cns_can_ratio = cns/can

        all_gene_coors_x.append(cns_can_ratio)
        all_gene_coors_y.append(pns_can_ratio)
        #categorizing dots

        if abs(pns_can_ratio) > 1 or abs(cns_can_ratio) > 1:
            if g in cns_can_upreg.keys():
                cns_specific_coors_x.append(cns_can_ratio)
                cns_specific_coors_y.append(pns_can_ratio)
            elif g in pns_can_upreg.keys():
                pns_specific_coors_x.append(cns_can_ratio)
                pns_specific_coors_y.append(pns_can_ratio)
            elif g in cns_pns_shared_genes.keys():
                pns_cns_shared_coors_x.append(cns_can_ratio)
                pns_cns_shared_coors_y.append(pns_can_ratio)
            elif g in cns_pns_codown_genes.keys():
                pns_cns_codown_coors_x.append(cns_can_ratio)
                pns_cns_codown_coors_y.append(pns_can_ratio)
            else:
                remain_gene_coors_x.append(cns_can_ratio)
                remain_gene_coors_y.append(pns_can_ratio)
        else:
            low_gene_coors_x.append(cns_can_ratio)
            low_gene_coors_y.append(pns_can_ratio)

        fit_gene_coors_x.append(cns_can_ratio)
        fit_gene_coors_y.append(pns_can_ratio)

    if g in annotated_genes:
        annotated_order.append(g)
        annotated_genes_coors_x.append(cns_can_ratio)
        annotated_genes_coors_y.append(pns_can_ratio)

        if g == 'Cx3cr1':
            print('pns exp : ', pns_exp)
            print('can exp : ', can_exp)
            print('cns exp : ', cns_exp)
            print(g,'|',cns_can_ratio,'|',pns_can_ratio)
        if g == 'Trem2':
            print('pns exp : ', pns_exp)
            print('can exp : ', can_exp)
            print('cns exp : ', cns_exp)
            print(g,'|',cns_can_ratio,'|',pns_can_ratio)
        if g == 'Apoe':
            print('pns exp : ', pns_exp)
            print('can exp : ', can_exp)
            print('cns exp : ', cns_exp)
            print(g,'|',cns_can_ratio,'|',pns_can_ratio)


x = pns_cns_shared_coors_x+pns_cns_codown_coors_x
y = pns_cns_shared_coors_y+pns_cns_codown_coors_y

# Modeling with Numpy
p, cov = np.polyfit(x, y, 1, cov=True)        # parameters and covariance from of the fit
y_model = np.polyval(p, x)                    # model using the fit parameters; NOTE: parameters here are coefficients

# Statistics
n = len(x)                              # number of observations
m = len(p)                                    # number of parameters
DF = n - m                                    # degrees of freedom
t = sp.stats.t.ppf(0.95, n - m)                  # used for CI and PI bands

# Estimates of Error in Data/Model
resid = y - y_model
chi2 = np.sum((resid/y_model)**2)             # chi-squared; estimates error in data
chi2_red = chi2/(DF)                          # reduced chi-squared; measures goodness of fit
s_err = np.sqrt(np.sum(resid**2)/(DF))        # standard deviation of the error


x2 = np.linspace(np.min(x), np.max(x), 100)
y2 = np.linspace(np.min(y_model), np.max(y_model), 100)


fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111)
ax.axis([math.floor(min(x))-0.5,math.ceil(max(x))+0.5,math.floor(min(y))-0.5,math.ceil(max(y))+0.5])
#ax.axis([-15,15,-12,12])
dot_size = 1.5
#ax.grid(False,which='both',color='grey',alpha=0.3)
ax.grid(False)
ax.axhline(y=0,alpha=0.8,color='k',linewidth=1)
ax.axvline(x=0,alpha=0.8,color='k',linewidth=1)
ax.set_facecolor('white')
plt.scatter(pns_specific_coors_x,pns_specific_coors_y, s=dot_size, alpha=0.9, c='purple')
plt.scatter(cns_specific_coors_x,cns_specific_coors_y, s=dot_size, alpha=0.7, c='green')
plt.scatter(pns_cns_shared_coors_x,pns_cns_shared_coors_y, s=dot_size, alpha=0.7, c='red')
plt.scatter(pns_cns_codown_coors_x,pns_cns_codown_coors_y, s=dot_size, alpha=0.7, c='blue')
plt.scatter(remain_gene_coors_x,remain_gene_coors_y,s=dot_size, alpha=0.3,c='black')
plt.scatter(low_gene_coors_x,low_gene_coors_y, s=1, alpha=.8, c='black')
#plt.scatter(gene_coors_x,gene_coors_y, s=dot_size, alpha=.3, c='black')
for i in range(0,len(annotated_order)):
    ax.annotate(annotated_order[i],(annotated_genes_coors_x[i],annotated_genes_coors_y[i]),fontsize=12)


#ax.plot(x,y_model,"-", color="0.1", linewidth=1.5, alpha=0.5, label="Fit")

# Confidence Interval (select one)
#plot_ci_manual(t, s_err, n, x, x2, y2, ax=ax)
#plot_ci_bootstrap(n, x, y, resid, ax=ax)

# Prediction Interval
#pi = t*s_err*np.sqrt(1+1/n+(x2-np.mean(x))**2/np.sum((x-np.mean(x))**2))
#ax.fill_between(x2, y2+pi, y2-pi, color="None", linestyle="--")
#ax.plot(x2, y2-pi, "--", color="0.5", label="95% Prediction Limits")
#ax.plot(x2, y2+pi, "--", color="0.5")


#plt.show()
plt.savefig('biplot_figs_prediction_interval_expressivity.jpg')

write_expression_matrix(header_str,diff_exp_matrix,cns_can_upreg,'plot_expressivity_cns_specific_genes.txt')
write_expression_matrix(header_str,diff_exp_matrix,pns_can_upreg,'plot_expressivity_pns_specific_genes.txt')
write_expression_matrix(header_str,diff_exp_matrix,cns_pns_shared_genes,'plot_expressivity_shared_cns_pns_genes.txt')
