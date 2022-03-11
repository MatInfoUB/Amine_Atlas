from sklearn.preprocessing import StandardScaler
from utils import rdkit_smiles_from_input_smiles
import pandas as pd
from umap import UMAP
from sklearn.manifold import TSNE
import plotly.express as px

df = pd.read_csv("processed_data/MHFP_0311.csv")
df['RDKIT_SMILES'] = df['SMILES'].map(lambda x: rdkit_smiles_from_input_smiles(x))
df = df.drop(columns=['SMILES']).set_index(["CID", "RDKIT_SMILES", "Type"])

scaler = StandardScaler()
np_scaled = scaler.fit_transform(df)

umap_model = UMAP(n_neighbors=50, min_dist=0.25, random_state=42)
version = "0311"
np_umap_result = umap_model.fit_transform(np_scaled)
df_umap_result = pd.DataFrame(np_umap_result, columns=['UMAP-1', 'UMAP-2'])
df = df.reset_index()
df = df[['CID', 'RDKIT_SMILES', 'Type']]
df_concat = pd.concat([df, df_umap_result], axis=1)
df_concat.to_csv("processed_data/UMAP_"+version+".csv", index=False)

# To do TSNE, uncomment this block
# tsne = TSNE(n_components=2, random_state=42, init='pca', perplexity=30, n_iter=1000)
# version2="0311"
# np_tsne_result = tsne.fit_transform(np_scaled)
# df_tsne_result = pd.DataFrame(np_tsne_result, columns=['TSNE-1', 'TSNE-2'])
# df = df.reset_index()
# features = df[['CID', 'RDKIT_SMILES', 'Type']]
# df_concat = pd.concat([df, df_tsne_result], axis=1)
# df_concat.to_csv("processed_data/TSNE_"+version2+".csv", index=False)


