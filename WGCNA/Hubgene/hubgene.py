import pandas as pd
import numpy as np
from scipy.stats import pearsonr
from sklearn.decomposition import PCA
import fastcluster

#pearson corr btwn gene expression & mod eigengene
def calcMember(expressionMatrix, modEigen, modAssignment):
    modMember = {}


    for gene in expressionMatrix.columns:
        module = modAssignment[gene]
        eigen = modEigen[module]
        corr, _ = pearsonr(expressionMatrix[gene], eigen)
        modMember[gene] = corr


    return pd.Series(modMember, name="moduleMember")


#pearson corr btwn expresson & trait vector
def calcSignificance(expressionMatrix, ronaVector):
    significance = {}
   
    for gene in expressionMatrix.columns:
        corr, _ = pearsonr(expressionMatrix[gene], ronaVector)
        significance[gene] = corr
   
    return pd.Series(significance, name="geneSignificance")


def idHub(modMember, significance, memberThresh=0.7, geneThresh=0.05):
    hubGenes = (modMember.abs() > memberThresh) & (significance < geneThresh)
    return hubGenes[hubGenes].index.tolist()


expressionMatrix = 

mod_member = calcMember(expressionMatrix, modEigen, modAssignment)
gene_sig = calcSignificance(expressionMatrix, ronaVector)
hub_genes = idHub(mod_member, gene_sig)


print("Hub Genes: ", hub_genes)


