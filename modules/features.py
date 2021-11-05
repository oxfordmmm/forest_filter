#!/usr/bin/env python3

feature_combinations={
       'QUAL only':['QUAL'],
       'composite':['QUAL','reads_all','proximty','baseChange','majority base %',
           'Top Base matches VC','deletions %','insertions %', 'strand bias'],
       'compositeINDELs':['QUAL','reads_all','proximty',
           'Top Base matches VC','deletions %','insertions %', 'strand bias',
           'INDEL length']
       }

