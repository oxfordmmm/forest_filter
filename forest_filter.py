#!/env python3
import sys
import pickle
from argparse import ArgumentParser
from modules.train import train
from modules.classify import classify
feature_combinations={
       'QUAL only':['QUAL'],
       'composite':['QUAL','reads_all','proximty','baseChange','majority base %','Top Base matches Nanopolish','deletions %','insertions %']
       }


###################### plots #####################################

def plot_feature_importances(importances,std,indices,features,feat):
    plt.figure()
    plt.title("Relative importance of features")
    df=pd.DataFrame({'Importance':importances,'std':std,'indices':indices,'feature':features})
    df.to_csv('{0}_feat_importance.csv'.format(feat.replace(' ','_')))
    g = sns.barplot(y='feature', x="Importance", xerr=df['std'], capsize=.2, data=df)
#    plt.bar(importances[indices],features
#            color="r", yerr=std[indices], align="center")
    #plt.xticks(features, indices)
    #plt.xlim([-1, len(features)])
    plt.tight_layout()
    plt.savefig('figs/{0}_feat_importance.pdf'.format(feat.replace(' ','_')))
    plt.savefig('figs/{0}_feat_importance.png'.format(feat.replace(' ','_')))
    plt.savefig('figs/{0}_feat_importance.svg'.format(feat.replace(' ','_')))
    plt.clf()

def plot_roc_curve(d):
    for i in d:
        plt.plot(d[i]['fpr'], d[i]['tpr'], label='{0}, AUC:{1:.2f}'.format(i,d[i]['AUC']))
    plt.plot([0, 1], [0, 1], color='darkblue', linestyle='--')
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver Operating Characteristic (ROC) Curve')
    plt.legend()
    plt.tight_layout()
    plt.savefig('figs/ROC.pdf')
    plt.savefig('figs/ROC.png')
    plt.savefig('figs/ROC.svg')
    #plt.show()
    plt.clf()

def plot_recall_precision(d):
    for i in d:
        average_precision = average_precision_score(d[i]['y_test'], d[i]['y_score'])
        precision, recall, _ = precision_recall_curve(d[i]['y_test'], d[i]['y_score'])
#        step_kwargs = ({'step': 'post'}
#               if 'step' in signature(plt.fill_between).parameters
#               else {})
        plt.step(recall, precision, where='post', label='{0} AP:{1:.2f}'.format(i,average_precision))
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.ylim([0.0, 1.05])
    plt.xlim([0.0, 1.0])
    plt.title('Precision-Recall curve')
    plt.legend()
    plt.tight_layout()
    plt.savefig('figs/PR.pdf')
    plt.savefig('figs/PR.png')
    plt.savefig('figs/PR.svg')

    plt.clf()

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



#################### train and classify ##############################

def maskProbs(r,feat='feat'):
    if r['{0} pred'.format(feat)]==False:
        return False
    elif r['{0} pred'.format(feat)] == True and r['{0} prob'.format(feat)] >= 0.95:
        return True
    elif r['{0} pred'.format(feat)] == True and r['{0} prob'.format(feat)] < 0.95:
        return False
    else:
        return None

def terms(r,col='classification'):
    if r['SNP validation'] == True and r[col] == True:
        return 'TP'
    elif r['SNP validation'] == True and r[col] == False:
        return 'FN'
    elif r['SNP validation'] == False and r[col] == False:
        return 'TN'
    elif r['SNP validation'] == False and r[col] == True:
        return 'FP'
    elif r['Origin'] == 'truth' and r['read depth'] < 1:
        return 'no reads'
    elif r['Origin'] == 'truth' and pd.isnull(r['read depth'])==True:
         return 'no reads'
    elif r['Origin'] == 'truth' and pd.isnull(r[col])==True:
        return 'missed'



def getModel(f):
    models=['fast','HAC','HAC-trained']
    for model in models:
        if model in f: m = model
    return m

def getStrain(f):
    strains=['WHOQ','WHOF','WHOX','WHOV','H18-208']
    for strain in strains:
        if strain in f: s = strain
    return s

def table_recall_precision(df,sample,depth,col='SNP type'):
    g=df.groupby([col])['POS'].count()
    g=g.reset_index()
    df2=g.T
    df2.columns = df2.iloc[0]
    df2=df2.reindex(df2.index.drop([col]))
    iterms=['TP','TN','FP','FN','missed','no reads']
    for term in iterms:
        if term not in df2:
            df2[term] = 0
    df2['FN']=df2['FN']+df2['missed']
    df2['Accuracy'] = (df2['TP'] + df2['TN']) / (df2['TP']+df2['FP']+df2['FN']+df2['TN'])
    df2['Precision']= df2['TP'] / (df2['TP']+df2['FP'])
    df2['Recall']   = df2['TP'] / (df2['TP']+df2['FN'])
    df2['F1 Score'] = 2*(df2['Recall'] * df2['Precision']) / (df2['Recall'] + df2['Precision'])
    df2['depth']=depth
    df2['Sample']=sample
    print(df2)
    df2['Strain']=df2['Sample'].map(getStrain)
    df2['Model']=df2['Sample'].map(getModel)
    df2['col']=col
    df2.to_csv('csvs/{0}_{1}_{2}_recall_precision.csv'.format(sample,col.replace(' ','_'),depth))
    return df2

def classifyEach(models,dfs,feat='composite'):
    features=feature_combinations[feat]
    model=models[feat]
    classifications,dfPR,dfPRFilt=[],[],[]
    for df in dfs:
        depth=df['sub depth'].to_list()[0]
        df2=df[features]
        df2=df2.dropna()
        X=np.array(df2)
        preds=model.predict(X)
        probs=model.predict_proba(X)
        probsdf=pd.DataFrame(probs,columns=[True,False])
        p=probsdf.max(axis=1)
        df2['{0} prob'.format(feat)] = p
        df2['{0} prob true'.format(feat)] = probs[:, 1]
        df2['{0} prob false'.format(feat)] = probs[:, 0]
        df2['{0} pred'.format(feat)]=preds
        df2=df2[['{0} prob'.format(feat),'{0} prob true'.format(feat),'{0} prob false'.format(feat),'{0} pred'.format(feat)]]
        df=df.join(df2,lsuffix='perGenome_')
        df['SNP type']=df.apply(terms,axis=1,col='{0} pred'.format(feat))
        df['{0} pred filtered'.format(feat)]=df.apply(maskProbs,axis=1,feat=feat)
        df['SNP type filtered']=df.apply(terms,axis=1,col='{0} pred filtered'.format(feat))
        #depth=int(df.reads_all.median())
        sample=df.Sample.unique()
        sample=[s for s in sample if pd.isnull(s) == False ]
        sample=sample[0]
        ## calculate precision recall etc for unfiltered SNPS
        dfPR.append(table_recall_precision(df,sample,depth))
        ## same for filtered SNPs
        dfPRFilt.append(table_recall_precision(df,sample,depth,col='SNP type filtered'))
        df.to_csv('csvs/{0}_{1}_classifications.csv'.format(sample,depth),index=False)
        classifications.append(df)

    df=pd.concat(dfPR)
    df2=pd.concat(dfPRFilt)
    alldf=pd.concat(classifications)
    alldf.to_csv('csvs/all_classifications.csv',index=False)
    return df,df2,alldf,classifications




############### args ###############


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

