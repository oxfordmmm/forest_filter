#!/usr/bin/env python3
import unittest
from modules.train import train
from modules.features import feature_combinations
import pandas as pd
import json

class TestCalc(unittest.TestCase):

    def test_getData(self):
        t=train()
        t.prefix='test'
        t.bams=['test/data/INDELs/MRSA_r9.10.sorted.bam']
        t.vcfs=['test/data/INDELs/MRSA_r9_10.vcf']
        t.refs=['test/data/INDELs/MRSA_mut.fasta']
        t.truths=['test/data/INDELs/MRSA_mut.csv']
        df,dfs=t.getData()
        df.sort_values(by='POS',inplace=True)
        df=df.reset_index(drop=True)
        #df.to_csv('test/data/INDELs/combined_results.csv')
        dfExpected=pd.read_csv('test/data/INDELs/combined_results.csv',index_col=0)
        pd._testing.assert_frame_equal(df, dfExpected)

    def test_train(self):
        t=train()
        t.prefix='test'
        df=pd.read_csv('test/data/INDELs/combined_results.csv',index_col=0)
        feat='compositeINDELs'
        features=feature_combinations[feat]
        d=t.trainTest(features,df,feat)
#        del d['model']
#        d['fpr']=d['fpr'].tolist()
#        d['tpr']=d['tpr'].tolist()
#        d['preds']=d['preds'].tolist()
#        d['probs']=d['probs'].tolist()
#        d['y_test']=d['y_test'].tolist()
#        d['y_score']=d['y_score'].tolist()
#        d['AUC']=float(d['AUC'])
#        #with open('test/data/train_results.json', 'w') as fp:
#        #    json.dump(d, fp)
#
#        with open('test/data/train_results.json', 'r') as fp:
#            exp=json.load(fp)
#
#        self.assertAlmostEqual(d['AUC'], exp['AUC'], delta=0.2)
#        #print(d['preds'])
#        #print(exp['preds'])
#        #self.assertCountEqual(d['preds'],exp['preds'])

if __name__ == '__main__':
    unittest.main()
