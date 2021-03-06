#!/env python3
import warnings
warnings.filterwarnings("ignore")
import os
import sys
from argparse import ArgumentParser
import pickle
import vcf
import numpy as np
import pandas as pd
from modules.utils import readVcf, addPysamstats, addFeatures, addTruth, plotFeatureImportances
from modules.features import feature_combinations
# roc curve and auc score
from sklearn.ensemble import RandomForestClassifier


class classify:
    def getClassifyArgs(self,parser):
        parser.add_argument('-v', '--inVCF', required=True,
                                 help='Input VCF files to be used for training or filtering')
        parser.add_argument('-o', '--outVCF', required=True,
                                 help='Input VCF files to be used for training or filtering')
        parser.add_argument('-r', '--reference', required=True,
                                 help='reference sequence fasta file')
        parser.add_argument('-m', '--model', required=True, 
                                 help='Random forest model to classify SNPs')
        parser.add_argument('-b', '--bam', required=True, 
                                 help='Sorted aligned bam file')
        parser.add_argument('-f', '--probFilt',required=False,default=0,
                             help='probability threshold filter')
        parser.add_argument('-mw', '--maskWeak', required=False,action='store_true',
                             help='mask weakly supported positions from pysam')
        parser.add_argument('-c', '--combination', required=False, default='composite',
                             help='name of the features to use, default=composite')
        return parser

    def run(self, opts):
        self.inVCF=opts.inVCF
        self.outVCF=opts.outVCF
        self.bam=opts.bam
        self.ref=opts.reference
        self.modelFile=opts.model
        self.combination=opts.combination
        self.keep=set()
        self.probFilt=float(opts.probFilt)
        self.maskWeak = opts.maskWeak

        # run
        self.loadModel()
        self.getData()
        self.classify()
        self.filter()

    def getData(self):
        df = readVcf(self.inVCF)
        df = addPysamstats(df,self.bam,self.ref)
        self.SNPs = addFeatures(df)
        self.SNPs.to_csv('test/data/example_classifier_data_full.csv')

    def loadModel(self):
        self.model = pickle.load(open(self.modelFile, 'rb'))

    def maskProbs(self,r):
        if r['probs'] >= self.probFilt:
            return False
        else:
            return True

    def classify(self):
        #??prep data
        features=feature_combinations[self.combination]
        X=np.array(self.SNPs[features])

        # predict
        preds = self.model.predict(X)
        probs = self.model.predict_proba(X)

        # rearrange data
        probs=pd.DataFrame(probs,columns=[True,False])
        p=probs.max(axis=1)
        self.SNPs['preds']=preds
        self.SNPs['probs']=p

        # apply masking
        self.SNPs['mask']=self.SNPs.apply(self.maskProbs,axis=1)
        keep=self.SNPs[self.SNPs.preds==True]
        s=set(keep['POS'])
        self.keep.update(s)
        mask=self.SNPs[self.SNPs['mask']==True]
        psmask=self.SNPs[self.SNPs['ps']<0.8]
        self.mask=set(mask['POS'])
        self.mask.update(list(psmask['POS']))

    def _vcf_reader(self):
        vcf_reader = vcf.Reader(open(self.inVCF, 'r'))
        for record in vcf_reader:
            yield record

    def filter(self):
        vcf_reader = self._vcf_reader()
        vcf_oneread= vcf.Reader(open(self.inVCF, 'r'))
        vcf_writer = vcf.Writer(open(self.outVCF, 'w'), vcf_oneread)
        for record in vcf_reader:
            if record.POS not in self.keep: continue
            else:
                if self.maskWeak == True:
                    if record.POS in self.mask:
                        record.ALT = 'N'
                        print('masking',record.POS)
                vcf_writer.write_record(record)
        vcf_writer.close()
