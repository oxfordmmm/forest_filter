# forest_filter
Using random forest classfiers to filter SNP calls

# Tests
 ![Tests](https://github.com/oxfordmmm/forest_filter/actions/workflows/test.yml/badge.svg)


## Install

Forest filter can be installed using Conda.

### Conda
To install the prerequesits needed to run forest_filter:
```bash
conda env create -f environment.yml
``` 

## Usage

Filter_forest has two componants, `train` and `classify`.

### Train
A model can be trained using the following command.

```bash
forest_filter.py train \
	-v MRSA_r9_10.vcf \
	-b MRSA_r9.10.sorted.bam \
	-r MRSA252_mut.fasta \
	-t MRSA252.fasta
```

This should improve upon "QUAL" only as shown in this diagram
![ROC](figs/ROC.png?raw=true "ROC")

The features and their importances for training are are.
![FEATURES](figs/composite_feat_importance.png?raw=true "Feature importances")



### Classify
The SNPs can be classified and therefore filtered with this command.

```bash
forest_filter.py classify \
	-v MRSA_r9_10.vcf \
	-b MRSA_r9.10.sorted.bam \
	-r MRSA252_mut.fasta \
	-m r9.4.1_composite_model.sav \
	-o filtered.vcf
```

## Testing

To run the tests, use these commands.

```bash
python3 -m unittest test/test_train.py
python3 -m unittest test/test_classify.py 
```

These tests are also run by github actions.
