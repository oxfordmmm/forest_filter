#!/bin/env python3
import warnings
warnings.filterwarnings("ignore")
import os
import sys
import pickle
from argparse import ArgumentParser
import numpy as np
import pandas as pd
import seaborn as sns
sns.set_style('darkgrid')
import matplotlib.pyplot as plt
from functools import reduce
import pysam
import pysamstats
from modules.utils import readVcf, addPysamstats, addFeatures
# roc curve and auc score
from sklearn.neighbors import KNeighborsClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_curve
from sklearn.metrics import precision_recall_curve
#from sklearn.utils.fixes import signature
from sklearn.metrics import average_precision_score
from sklearn.metrics import roc_auc_score
from sklearn import svm
from sklearn.model_selection import StratifiedKFold
from itertools import permutations



### train #####

class train:
    def getData(self):
        dfs=[]
        for vcf,bam,ref in zip(self.vcfs,self.bams,self.refs):
            df = readVcf(vcf)
            df = addPysamstats(df,bam,ref)
            dfs.append(df)
            df = addFeatures(df)

        df=pd.concat(dfs)
        print(df)
        return df,dfs

    def runTrain(self):
        df,allFrames=self.getData()
        #df=df[~df['SNP validation'].isna()]
        #d={}
        #dfs=[]
        #models={}
        #for feat in feature_combinations:
        #    features=feature_combinations[feat]
        #    d[feat]=slef.trainTest(features,df,feat)
        #    models[feat]=d[feat]['model']
    
        #plot_roc_curve(d)
        #plot_recall_precision(d)
        #return models,allFrames
    
    def trainTest(self, features,df,feat):
        # train and test
        features.append('SNP validation')
        print(features)
        #l=df['POS'].max()
        #df2=df[df['POS'] < int((l/100)*30)]
        df2=df[features]
        df2=df2.dropna()
        features.remove('SNP validation')
        print(features)
        X=np.array(df2[features])
        Y=np.array(df2['SNP validation'])
        Y=Y.astype('int')
        trainX, testX, trainy, testy = train_test_split(X, Y, test_size=0.3, random_state=1)
        model = RandomForestClassifier()
        #model = GradientBoostingClassifier()
        model.fit(trainX, trainy)
        prob = model.predict_proba(testX)
        probs = prob[:, 1]
        auc = roc_auc_score(testy, probs)
        fpr, tpr, thresholds = roc_curve(testy, probs)
        filename = '{0}_model.sav'.format(feat.replace(' ','_'))
        pickle.dump(model, open(filename, 'wb'))
        saveTree(model,features)
        # classify all
        preds= model.predict(X)
        probs=model.predict_proba(X)
        importances = model.feature_importances_
        std = np.std([tree.feature_importances_ for tree in model.estimators_],axis=0)
        indices = np.argsort(importances)[::-1]
        plot_feature_importances(importances,std,indices,features,feat)
        #probs=probs[:, 1]
        return {'fpr':fpr,'tpr':tpr,'AUC':auc,'probs':probs,'preds':preds,
                'y_test':Y,'y_score':preds,'model':model}
    
    def run(self,opts):
        self.vcfs=opts.vcf_files
        self.bams=opts.bam_files
        self.refs=opts.ref_files
        self.runTrain()

    def getTrainArgs(self,parser):
        parser.add_argument('-v', '--vcf_files', required=True, nargs='+',
                                 help='Input VCF files to be used for training or filtering')
        parser.add_argument('-f', '--ref_files', required=True, nargs='+',
                                 help='reference fasta files for training')
        parser.add_argument('-b', '--bam_files', required=True, nargs='+',
                                 help='pysamstats files')
        return parser
