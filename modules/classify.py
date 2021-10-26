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


class classify:
    def getClassifyArgs(self,parser):
        parser.add_argument('-v', '--vcf_files', required=True, nargs='+',
                                 help='Input VCF files to be used for training or filtering')
        parser.add_argument('-p', '--pysamstats_files', required=True, nargs='+',
                                 help='pysamstats files')
        parser.add_argument('-m', '--model', required=False, 
                                 help='Random forest model to classify SNPs')
        return parser

    def run():
        ...
