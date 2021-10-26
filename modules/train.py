#!/env python3
import warnings
warnings.filterwarnings("ignore")
import os
import sys
from argparse import ArgumentParser
import numpy as np
import pandas as pd
import seaborn as sns
sns.set_style('darkgrid')
import matplotlib.pyplot as plt
from functools import reduce
# roc curve and auc score
#from sklearn.neighbors import KNeighborsClassifier
#from sklearn.ensemble import RandomForestClassifier
#from sklearn.ensemble import GradientBoostingClassifier
#from sklearn.model_selection import train_test_split
#from sklearn.metrics import roc_curve
#from sklearn.metrics import precision_recall_curve
##from sklearn.utils.fixes import signature
#from sklearn.metrics import average_precision_score
#from sklearn.metrics import roc_auc_score
#from sklearn import svm
#from sklearn.model_selection import StratifiedKFold
#from itertools import permutations
#import pickle

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

def plot_recall_depth(df,prefix='unfiltered'):
    g=sns.scatterplot('depth','Recall',hue='Sample',data=df)
    #g.legend(loc='bottom left', bbox_to_anchor=(1.25, 0.5), ncol=1)
    #plt.show()
    plt.savefig('figs/{0}_depth_recall.png'.format(prefix))
    plt.savefig('figs/{0}_depth_recall.pdf'.format(prefix))
    plt.savefig('figs/{0}_depth_recall.svg'.format(prefix))
    plt.clf()

def plot_FP_depth(df,prefix='unfiltered'):
    g=sns.scatterplot('depth','FP',hue='Model',style='Strain',data=df)
    #g.legend(loc='bottom left', bbox_to_anchor=(1.25, 0.5), ncol=1)
    #plt.show()
    plt.savefig('figs/{0}_depth_FP.png'.format(prefix))
    plt.savefig('figs/{0}_depth_FP.pdf'.format(prefix))
    plt.savefig('figs/{0}_depth_FP.svg'.format(prefix))
    plt.clf()


def plot_recall_depth_comp(df,prefix='unfiltered'):
    g=sns.scatterplot('depth','Recall',hue='Model',style='Strain',data=df,s=50)
    g.legend(loc='lower right',ncol=2)
    g.set(ylim=(0, 1))
    g.set(xlim=(0, 130))
    plt.savefig('figs/{0}_depth_recall.png'.format(prefix))
    plt.savefig('figs/{0}_depth_recall.pdf'.format(prefix))
    plt.savefig('figs/{0}_depth_recall.svg'.format(prefix))
    #plt.show()
    plt.clf()

def plot_TN_depth_comp(df,prefix='unfiltered'):
    g=sns.scatterplot('depth','TN',hue='Model',style='Strain',data=df,s=50)
    #g.legend(loc='upper left',bbox_to_anchor=(1.04,1), ncol=1)
    g.set(ylim=(0, None))
    g.set(xlim=(0, 130))
    #plt.tight_layout()
    plt.savefig('figs/{0}_TN_recall.png'.format(prefix))
    plt.savefig('figs/{0}_TN_recall.pdf'.format(prefix))
    plt.savefig('figs/{0}_TN_recall.svg'.format(prefix))
    #plt.show()
    plt.clf()

def plot_accuracy_depth_comp(df,prefix='unfiltered'):
    g=sns.scatterplot('depth','Accuracy',hue='Model',style='Strain',data=df,s=50)
    g.legend(loc='lower right',ncol=2)
    g.set(ylim=(0, 1))
    g.set(xlim=(0, 130))
    plt.savefig('figs/{0}_Accuracy_recall.pdf'.format(prefix))
    plt.savefig('figs/{0}_Accuracy_recall.svg'.format(prefix))
    plt.savefig('figs/{0}_Accuracy_recall.png'.format(prefix))
    #plt.show()
    plt.clf()

def plot_F1_depth_comp(df,prefix='unfiltered'):
    g=sns.scatterplot('depth','F1 Score',hue='Model',style='Strain',data=df,s=50)
    g.legend(loc='lower right',ncol=2)
    g.set(ylim=(0, 1))
    g.set(xlim=(0, 130))
    plt.savefig('figs/{0}_F1_depth.png'.format(prefix))
    plt.savefig('figs/{0}_F1_depth.pdf'.format(prefix))
    plt.savefig('figs/{0}_F1_depth.svg'.format(prefix))
    #plt.show()
    plt.clf()

def plot_prob_box(df):
    df=df[df['SNP type'] != 'missed']
    df=df.sample(5000)
    g=sns.FacetGrid(df,row='SNP type',height=2, aspect=4)
    g=g.map(sns.distplot,"composite prob")
    plt.savefig('figs/probs_distplot.png')
    plt.savefig('figs/probs_distplot.pdf')
    plt.savefig('figs/probs_distplot.svg')
    #plt.show()
    plt.clf()

class train:
    def getData(self):
        p='snps/'
        runs=os.listdir(p)
        runs=[r for r in runs if r.endswith('.csv')]
        #runs=[r for r in runs if 'fast' in r]
        dfs=[]
        for run in runs:
            csv='{0}{1}'.format(p,run)
            df=pd.read_csv(csv)
            dfs.append(df)
        df=pd.concat(dfs)
        return df,dfs

    def runTrain(self):
        df,allFrames=self.getData()
        df=df[~df['SNP validation'].isna()]
        d={}
        dfs=[]
        models={}
        for feat in feature_combinations:
            features=feature_combinations[feat]
            d[feat]=slef.trainTest(features,df,feat)
            models[feat]=d[feat]['model']
    
        plot_roc_curve(d)
        plot_recall_precision(d)
        return models,allFrames
    
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
        self.vcfs=opts.vcfs
        self.pysams=opts.pysams
        self.refs=opts.refs
        models,dfs=runTrain()

    def getTrainArgs(self,parser):
        parser.add_argument('-v', '--vcf_files', required=True, nargs='+',
                                 help='Input VCF files to be used for training or filtering')
        parser.add_argument('-f', '--ref_files', required=True, nargs='+',
                                 help='reference fasta files for training')
        parser.add_argument('-p', '--pysamstats_files', required=True, nargs='+',
                                 help='pysamstats files')
        return parser
