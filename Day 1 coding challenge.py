import pandas as pd
import numpy as np


df = pd.read_csv("GSE293745_gene_reads_count.anno.csv")  


# Calculate mean expression across experimental samples
df['mean_patient'] = df[['C_P1', 'C_P2', 'C_P3']].mean(axis=1)

# Calculate mean expression across control samples
df['mean_control'] = df[['C_V1', 'C_V2', 'C_V3']].mean(axis=1)


# Adding a small constant to avoid division by zero
df['log2_fc'] = np.log2((df['mean_patient'] + 1e-9) / (df['mean_control'] + 1e-9))


# Sort by absolute log2 fold change, descending
df_top = df.reindex(df['log2_fc'].abs().sort_values(ascending=False).index)
# Show top 10 genes with highest expression difference
print("----------------------------------------------------------------")
print(df_top[['gene_symbol', 'mean_patient', 'mean_control', 'log2_fc']].head(10))

# Pick a gene from top list
gene = df_top.iloc[0]['gene_symbol']

# Filter row for that gene
gene_data = df[df['gene_symbol'] == gene][['C_P1', 'C_P2', 'C_P3', 'C_V1', 'C_V2', 'C_V3']].iloc[0]



print("----------------------------------------------------------------")
# begining of go enrichment analysis
# Define upregulated and downregulated genes
upregulated = df[df['log2_fc'] > 1]['gene_symbol'].dropna().unique().tolist()
downregulated = df[df['log2_fc'] < -1]['gene_symbol'].dropna().unique().tolist()

from gprofiler import GProfiler

gp = GProfiler(return_dataframe=True)

# Run enrichment on upregulated genes
go_results_up = gp.profile(organism='hsapiens', query=upregulated, sources=['GO:BP', 'GO:MF', 'GO:CC'])

# Run enrichment on downregulated genes
go_results_down = gp.profile(organism='hsapiens', query=downregulated, sources=['GO:BP', 'GO:MF', 'GO:CC'])

# Preview top results
print("Top GO terms in upregulated genes:")
print(go_results_up[['name', 'p_value', 'term_size', 'intersection_size']].head())
print("----------------------------------------------------------------")
print("\nTop GO terms in downregulated genes:")
print(go_results_down[['name', 'p_value', 'term_size', 'intersection_size']].head())
print("----------------------------------------------------------------")


#completing enrichr 
# Threshold for upregulation (e.g., log2_fc > 1)
upregulated = df[df['log2_fc'] > 1].sort_values('log2_fc', ascending=False)

top200_genes = upregulated['gene_symbol'].dropna().drop_duplicates().head(200)

# Save to .txt
top200_genes.to_csv("top200_upregulated.txt", index=False, header=False)
df[df["gene_symbol"] == "EZH2"]
print(df[df['gene_symbol'].str.contains("EZH2", case=False, na=False)])

print("----------------------------------------------------------------")
# Threshold for downregulation (e.g., log2_fc < -1)
downregulated = df[df['log2_fc'] < -1].sort_values('log2_fc')


top200_down_genes = downregulated['gene_symbol'].dropna().drop_duplicates().head(200)

# Save to a txt file
top200_down_genes.to_csv("top200_downregulated.txt", index=False, header=False)

# preview to see if it worked
print("Top downregulated genes:")
print(top200_down_genes.head())

