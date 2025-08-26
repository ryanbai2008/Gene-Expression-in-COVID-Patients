import pandas as pd
from sklearn.decomposition import PCA
from scipy.stats import pearsonr


expr_df = pd.read_csv("data/truncated normal counts.csv", index_col=0)
if expr_df.shape[0] > expr_df.shape[1]:
    expr_df = expr_df.T
TOM_diss_df = pd.read_csv("data/tom_matrix_values.csv", index_col=0)
covid_traits_df = pd.read_csv("data/traits.csv", index_col=0)
dendrogram_df = pd.read_csv("WGCNA/Dendrogram/dendrogram_clusters.csv")

common_samples = expr_df.index.intersection(covid_traits_df.index)
expr_df = expr_df.loc[common_samples]
covid_traits_df = covid_traits_df.loc[common_samples]

module_df = {}
for i in range(1, 751):
    cluster_genes = dendrogram_df[dendrogram_df["Cluster"] == i]["Gene"].tolist()
    cluster_genes = [g for g in cluster_genes if g in expr_df.columns]
    module_df[i] = expr_df[cluster_genes]


pca = PCA(n_components=1)
eigengenes = {}
for module in module_df:
    data = module_df[module]
    if data.shape[1] == 0:
        eigengenes[module] = pd.Series([float('nan')] * data.shape[0], index=data.index)
    elif data.shape[1] == 1:
        eigengenes[module] = data.iloc[:, 0]
    else:
        eig = pca.fit_transform(data)
        eigengenes[module] = pd.Series(eig[:, 0], index=data.index)


COVID_correlation = {}
for module, eigengene in eigengenes.items():
    valid_idx = eigengene.dropna().index.intersection(covid_traits_df.index)
    r, p = pearsonr(eigengene.loc[valid_idx], covid_traits_df.loc[valid_idx, "Status"])
    COVID_correlation[module] = {"correlation": r, "p_value": p}


def calcMember(expressionMatrix, modEigen, modAssignment):
    modMember = {}
    for gene in dendrogram_df["Gene"]:
        if gene not in expressionMatrix.columns:
            continue
        module = modAssignment[gene]
        eigen = modEigen[module]
        corr, _ = pearsonr(expressionMatrix[gene], eigen)
        modMember[gene] = corr
    return pd.Series(modMember, name="moduleMembership")

def calcSignificance(expressionMatrix, traitVector):
    significance = {}
    pvalues = {}
    for gene in dendrogram_df["Gene"]:
        if gene not in expressionMatrix.columns:
            continue
        corr, p = pearsonr(expressionMatrix[gene], traitVector)
        significance[gene] = corr
        pvalues[gene] = p
    return pd.Series(significance, name="geneSignificance"), pd.Series(pvalues, name="geneSigPvalue")


def idHub(modMember, significance, pvalues, memberThresh=0.7, sigThresh=0.2, pvalThresh=0.05):
    hubGenes = (
        (modMember.abs() > memberThresh) &
        (significance.abs() > sigThresh) &
        (pvalues < pvalThresh)
    )
    return hubGenes[hubGenes].index.tolist()

expressionMatrix = expr_df
modAssignment = dict(zip(dendrogram_df["Gene"], dendrogram_df["Cluster"]))

modEigen = {}
for module_id, eig in eigengenes.items():
    modEigen[module_id] = eig

ronaVector = covid_traits_df.loc[expressionMatrix.index, "Status"]

mod_member = calcMember(expressionMatrix, modEigen, modAssignment)
gene_sig, gene_sig_pval = calcSignificance(expressionMatrix, ronaVector)
hub_genes = idHub(mod_member, gene_sig, gene_sig_pval, memberThresh=0.95, sigThresh=0.03, pvalThresh=0.05)

pd.DataFrame(eigengenes).to_csv("module_eigengenes.csv")
pd.DataFrame(COVID_correlation).T.to_csv("module_trait_correlation.csv")
mod_member.to_csv("module_membership.csv")
gene_sig.to_csv("gene_significance.csv")
gene_sig_pval.to_csv("gene_significance_pval.csv")


# Convert TOM dissimilarity to similarity
TOM_sim = 1 - TOM_diss_df.values
genes = TOM_diss_df.columns

threshold = 0.1  # adjust based on network density
edges_list = []
for i in range(len(genes)):
    for j in range(i+1, len(genes)):
        if TOM_sim[i, j] > threshold:
            edges_list.append([genes[i], genes[j], TOM_sim[i, j]])

edges = pd.DataFrame(edges_list, columns=["Gene1", "Gene2", "Weight"])

edges.to_csv("module_network_edges.csv")

print("Hub Genes:", hub_genes)
hub_genes = pd.DataFrame(hub_genes, columns=["Hub Genes"])
hub_genes.to_csv("data/hub_genes.csv")

