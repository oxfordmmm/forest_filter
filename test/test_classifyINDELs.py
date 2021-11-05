#!/usr/bin/env python3
import unittest
from modules.classify import classify
from modules.features import feature_combinations
import pandas as pd
import json

class TestCalc(unittest.TestCase):

    def test_getData(self):
        t=classify()
        t.bam='test/data/INDELs/MRSA_r9.10.sorted.bam'
        t.inVCF='test/data/INDELs/MRSA_r9_10.vcf'
        t.ref='test/data/INDELs/MRSA_mut.fasta'
        t.mode='INDELs_only'
        t.combination='compositeINDELs'
        t.getData()
        #t.SNPs.to_csv('test/data/INDELs/example_classifier_data.csv')
        dfExpected=pd.read_csv('test/data/INDELs/example_classifier_data.csv',index_col=0)
        pd._testing.assert_frame_equal(t.SNPs, dfExpected)

    def test_classifier(self):
        t=classify()
        t.combination='compositeINDELs'
        t.mode='INDELs_only'
        t.keep=set()
        t.probFilt=0
        t.maskWeak = False

        t.modelFile='test/data/INDELs/Medaka_INDELs_MRSA_r9_compositeINDELs_model.sav'
        t.SNPs=pd.read_csv('test/data/INDELs/example_classifier_data.csv',index_col=0)
        t.loadModel()
        t.classify()
        #t.SNPs.to_csv('test/data/INDELs/example_classified_data_full.csv')
        dfExpected=pd.read_csv('test/data/INDELs/example_classified_data_full.csv',index_col=0)
        pd._testing.assert_frame_equal(t.SNPs, dfExpected)

if __name__ == '__main__':
    unittest.main()
