#!/bin/env python3
import os
import sys
import pickle
from argparse import ArgumentParser
import numpy as np
import pandas as pd
import seaborn as sns
sns.set_style('darkgrid')
import matplotlib.pyplot as plt
import pysam
import pysamstats
from Bio import SeqIO
from itertools import permutations


## file loading
def readVcf(vcf):
    vcf_cols=['CHROM',
            'POS',
            'ID',
            'REF',
            'ALT',
            'QUAL',
            'FILTER',
            'INFO',
            'FORMAT',
            'SAMPLE']
    df=pd.read_csv(vcf,
        comment='#',
        sep='\t',
        names=vcf_cols)
    return df

def addPysamstats(vcf,bam,ref):
    mybam = pysam.AlignmentFile(bam)
    l=[]
    pos=vcf['POS'].unique()
    for rec in pysamstats.stat_variation_strand(mybam, ref):
        if rec['pos'] in pos:
            l.append(rec)
    df=pd.DataFrame(l)
    vcf=vcf.merge(df,left_on=['CHROM','POS'],right_on=['chrom','pos'],how='left')
    return vcf

def addFeatures(df):
    df['ALT_len']=df.ALT.map(len)
    df['REF_len']=df.REF.map(len)
    df['INDEL length']=df['ALT_len'] -  df['REF_len']

    df['5prime proximity']=df['POS']-df['POS'].shift(1)
    df['3prime proximity']=df['POS']-df['POS'].shift(-1)
    df['3prime proximity']=df['3prime proximity'].abs()
    df['proximty']=df[['5prime proximity','3prime proximity']].min(axis=1)
    df['proximty'].fillna(5000,inplace=True)
    df['alphabeta']=df.apply(alphaBeta,axis=1)
    df=df[df.REF_len==1]
    df=df[df.ALT_len==1]

    if len(df) > 0:
        df['baseChange']=df.apply(addBaseChangeN,axis=1)

    df['total']=df['A']+df['C']+df['G']+df['T']
    bases=['A','T','C','G','insertions','deletions']
    for b in bases:
        df['{0} %'.format(b)]=(df[b]/df['reads_all'])*100
        df['{0}_PS'.format(b)]=(df[b]/df['total'])
    df['reads_all']=df[['A','T','C','G']].sum(axis=1)
    df['top_base'] = df[['A','T','C','G']].max(axis=1)
    df['top_base_seq'] = df[['A','T','C','G']].idxmax(axis=1)
    df['majority base %'] = (df['top_base'] / df['reads_all'])
    df['Top Base matches VC'] = np.where(df.ALT == df.top_base_seq,1,0)
    df=df[df['majority base %'].notna()]
    if len(df)>0:
        df['ps']=df.apply(ps, axis=1)
        df['strand bias']=df.apply(strandBias,axis=1)
    return df

def addTruth(df, truth):
    seqs={}
    with open(truth,'rt') as inf:
        for seq in SeqIO.parse(inf,'fasta'):
            seqs[str(seq.id)]=str(seq.seq)
    df['SNP validation']=df.apply(checkRef,axis=1,seqs=seqs)
    return df


def checkRef(row,seqs=None):
    chrom=str(row['CHROM'])
    pos  =int(row['POS']) -1 
    ref_base=seqs[chrom][pos]
    if ref_base==row['ALT'][0]:
        return True
    else:
        return False



def totalReads(row):
    l=row['INFO'].split(';')
    return int(l[2].replace('TotalReads=',''))

def supportFraction(row):
    l=row['INFO'].split(';')
    if len(l)>2:
        return float(l[4].replace('SupportFraction=',''))
    else:
        return None

def BaseCalledReadsWithVariant(row):
    l=row['INFO'].split(';')
    if len(l)>2:
        return float(l[0].replace('BaseCalledReadsWithVariant=',''))
    else:
        return None

def BaseCalledFraction(row):
    l=row['INFO'].split(';')
    if len(l)>2:
        return float(l[1].replace('BaseCalledFraction=',''))
    else:
        return None

basesN=['A', 'C', 'G', 'T']
perm = permutations(basesN,2)
N,perms=0,{}
for p in perm:
    perms.setdefault(p[0],{}).setdefault(p[1],N)
    N+=1

def addBaseChangeN(row):
    r=row['REF']
    a=row['ALT']
    p=perms[r][a]
    return p

def alphaBeta(row):
    if row['REF'] == 'A' and row['ALT'] == 'G':
        return 0
    elif row['REF'] == 'A' and row['ALT'] == 'T':
        return 1
    elif row['REF'] == 'A' and row['ALT'] == 'C':
        return 1
    elif row['REF'] == 'C' and row['ALT'] == 'A':
        return 1
    elif row['REF'] == 'C' and row['ALT'] == 'G':
        return 1
    elif row['REF'] == 'C' and row['ALT'] == 'T':
        return 0
    elif row['REF'] == 'G' and row['ALT'] == 'A':
        return 0
    elif row['REF'] == 'G' and row['ALT'] == 'C':
        return 1
    elif row['REF'] == 'G' and row['ALT'] == 'T':
        return 1
    elif row['REF'] == 'T' and row['ALT'] == 'A':
        return 1
    elif row['REF'] == 'T' and row['ALT'] == 'C':
        return 0
    elif row['REF'] == 'T' and row['ALT'] == 'G':
        return 1

def ps(row):
    if row['ALT']=='N':
        x=None
    else:
        x = float(row[row['ALT']+'_PS'])
    return x

def strandBias(r):
    if r['reads_rev']==0:
        rs=0
    else:
        rs=r['reads_rev']/r['reads_all']
    if r['reads_fwd']==0:
        fs=0
    else:
        fs=r['reads_fwd']/r['reads_all']
    return 1*fs




## plots

def plotFeatureImportances(importances,std,indices,features,feat):
    plt.figure()
    plt.title("Relative importance of features")
    df=pd.DataFrame({'Importance':importances,'std':std,'indices':indices,'feature':features})
    df.to_csv('{0}_feat_importance.csv'.format(feat.replace(' ','_')))
    g = sns.barplot(y='feature', x="Importance", xerr=df['std'], capsize=.2, data=df)
#    plt.bar(importances[indices],features
#            color="r", yerr=std[indices], align="center")
    #plt.xticks(features, indices)
    #plt.xlim([-1, len(features)])
    plt.tight_layout()
    plt.savefig('{0}_feat_importance.pdf'.format(feat.replace(' ','_')))
    plt.savefig('{0}_feat_importance.png'.format(feat.replace(' ','_')))
    plt.savefig('{0}_feat_importance.svg'.format(feat.replace(' ','_')))
    plt.clf()

def plotRocCurve(d):
    for i in d:
        plt.plot(d[i]['fpr'], d[i]['tpr'], label='{0}, AUC:{1:.2f}'.format(i,d[i]['AUC']))
    plt.plot([0, 1], [0, 1], color='darkblue', linestyle='--')
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver Operating Characteristic (ROC) Curve')
    plt.legend()
    plt.tight_layout()
    plt.savefig('ROC.pdf')
    plt.savefig('ROC.png')
    plt.savefig('ROC.svg')
    #plt.show()
    plt.clf()

def plotRecallPrecision(d):
    for i in d:
        average_precision = average_precision_score(d[i]['y_test'], d[i]['y_score'])
        precision, recall, _ = precision_recall_curve(d[i]['y_test'], d[i]['y_score'])
#        step_kwargs = ({'step': 'post'}
#               if 'step' in signature(plt.fill_between).parameters
#               else {})
        plt.step(recall, precision, where='post', label='{0} AP:{1:.2f}'.format(i,average_precision))
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.ylim([0.0, 1.05])
    plt.xlim([0.0, 1.0])
    plt.title('Precision-Recall curve')
    plt.legend()
    plt.tight_layout()
    plt.savefig('PR.pdf')
    plt.savefig('PR.png')
    plt.savefig('PR.svg')

    plt.clf()


def saveTree(model,features):
    ## Extract single tree
    estimator = model.estimators_[5]
    from sklearn.tree import export_graphviz
    # Export as dot file
    export_graphviz(estimator, out_file='tree.dot',
                feature_names = features,
                class_names = ['False','True'],
                rounded = True, proportion = False,
                precision = 2, filled = True)

def plotRecallDepth(df,prefix='unfiltered'):
    g=sns.scatterplot('depth','Recall',hue='Sample',data=df)
    #g.legend(loc='bottom left', bbox_to_anchor=(1.25, 0.5), ncol=1)
    #plt.show()
    plt.savefig('{0}_depth_recall.png'.format(prefix))
    plt.savefig('{0}_depth_recall.pdf'.format(prefix))
    plt.savefig('{0}_depth_recall.svg'.format(prefix))
    plt.clf()

def plot_FP_depth(df,prefix='unfiltered'):
    g=sns.scatterplot('depth','FP',hue='Model',style='Strain',data=df)
    #g.legend(loc='bottom left', bbox_to_anchor=(1.25, 0.5), ncol=1)
    #plt.show()
    plt.savefig('{0}_depth_FP.png'.format(prefix))
    plt.savefig('{0}_depth_FP.pdf'.format(prefix))
    plt.savefig('{0}_depth_FP.svg'.format(prefix))
    plt.clf()


def plot_recall_depth_comp(df,prefix='unfiltered'):
    g=sns.scatterplot('depth','Recall',hue='Model',style='Strain',data=df,s=50)
    g.legend(loc='lower right',ncol=2)
    g.set(ylim=(0, 1))
    g.set(xlim=(0, 130))
    plt.savefig('{0}_depth_recall.png'.format(prefix))
    plt.savefig('{0}_depth_recall.pdf'.format(prefix))
    plt.savefig('{0}_depth_recall.svg'.format(prefix))
    #plt.show()
    plt.clf()

def plot_TN_depth_comp(df,prefix='unfiltered'):
    g=sns.scatterplot('depth','TN',hue='Model',style='Strain',data=df,s=50)
    #g.legend(loc='upper left',bbox_to_anchor=(1.04,1), ncol=1)
    g.set(ylim=(0, None))
    g.set(xlim=(0, 130))
    #plt.tight_layout()
    plt.savefig('{0}_TN_recall.png'.format(prefix))
    plt.savefig('{0}_TN_recall.pdf'.format(prefix))
    plt.savefig('{0}_TN_recall.svg'.format(prefix))

