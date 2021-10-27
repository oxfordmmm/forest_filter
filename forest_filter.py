#!/env python3
import sys
import pickle
from argparse import ArgumentParser
from modules.train import train
from modules.classify import classify



if __name__=="__main__":
    # args
    parser = ArgumentParser(description='Train Random forest classifier and filter SNPs')

    subparsers = parser.add_subparsers(dest='subparser')

    # train and test
    trainParse = subparsers.add_parser('train', help='Train Random Forest on VCFs, reference sequences and bam files')
    train=train()
    trainParse = train.getTrainArgs(trainParse)
    trainParse.set_defaults(func=train.run)

    # classify
    classifyParse = subparsers.add_parser('classify', help='Filter VCF files using trained Random Forest model')
    classify=classify()
    classifyParse = classify.getClassifyArgs(classifyParse)
    classifyParse.set_defaults(func=classify.run)

    args = parser.parse_args()
    if hasattr(args,'func'):
        args.func(args)

