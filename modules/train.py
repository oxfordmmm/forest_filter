#!/bin/env python3
import os
import sys
import pickle
from argparse import ArgumentParser
import numpy as np
import pandas as pd
from functools import reduce
from itertools import permutations
from modules.utils import readVcf, addPysamstats, addFeatures, addTruth, plotFeatureImportances 
from modules.utils import plotFeatureImportances, plotRocCurve, plotRecallPrecision
from modules.features import feature_combinations
# roc curve and auc score
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_curve
from sklearn.metrics import precision_recall_curve
#from sklearn.utils.fixes import signature
from sklearn.metrics import average_precision_score
from sklearn.metrics import roc_auc_score



### train #####

class train:
    def run(self,opts):
        self.vcfs=opts.vcf_files
        self.bams=opts.bam_files
        self.refs=opts.ref_files
        self.truths=opts.truth_files
        self.prefix=opts.prefix
        self.runTrain()

    def runTrain(self):
        df,allFrames=self.getData()
        df=df[~df['SNP validation'].isna()]
        d={}
        dfs=[]
        models={}
        for feat in feature_combinations:
            features=feature_combinations[feat]
            d[feat]=self.trainTest(features,df,feat)
            models[feat]=d[feat]['model']
    
        plotRocCurve(d)

    def getData(self):
        dfs=[]
        for vcf,bam,ref,truth in zip(self.vcfs,self.bams,self.refs,self.truths):
            df = readVcf(vcf)
            df = addPysamstats(df,bam,ref)
            df = addFeatures(df)
            df = addTruth(df,truth)
            dfs.append(df)

        df=pd.concat(dfs)
        #df.to_csv('test/data/train_test_data.csv')
        return df,dfs

    
    def trainTest(self, features,df,feat):
        # prepare train and test data
        features.append('SNP validation')
        df2=df[features]
        df2=df2.dropna()
        
        features.remove('SNP validation')
        X=np.array(df2[features])
        Y=np.array(df2['SNP validation'])
        Y=Y.astype('int')
        trainX, testX, trainy, testy = train_test_split(X, Y, test_size=0.3, random_state=1)

        # prepare model
        model = RandomForestClassifier()
        model.fit(trainX, trainy)
        prob = model.predict_proba(testX)
        probs = prob[:, 1]
        auc = roc_auc_score(testy, probs)
        fpr, tpr, thresholds = roc_curve(testy, probs)
        filename = '{0}_{1}_model.sav'.format(self.prefix,feat.replace(' ','_'))
        pickle.dump(model, open(filename, 'wb'))
   
        # classify all
        preds= model.predict(X)
        probs=model.predict_proba(X)
        importances = model.feature_importances_
        std = np.std([tree.feature_importances_ for tree in model.estimators_],axis=0)
        indices = np.argsort(importances)[::-1]
        plotFeatureImportances(importances,std,indices,features,feat)
        probs=probs[:, 1]
        return {'fpr':fpr,'tpr':tpr,'AUC':auc,'probs':probs,'preds':preds,
                'y_test':Y,'y_score':preds,'model':model}
    

    def getTrainArgs(self,parser):
        parser.add_argument('-v', '--vcf_files', required=True, nargs='+',
                                 help='Input VCF files to be used for training or filtering')
        parser.add_argument('-r', '--ref_files', required=True, nargs='+',
                                 help='reference fasta files for training')
        parser.add_argument('-t', '--truth_files', required=True, nargs='+',
                                 help='truth fasta files for training')
        parser.add_argument('-b', '--bam_files', required=True, nargs='+',
                                 help='Aligned and sorted bam file with index')
        parser.add_argument('-p', '--prefix', required=False, default=None,
                                 help='Prefix for model')
        return parser
