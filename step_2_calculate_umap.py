from sklearn.preprocessing import StandardScaler
from utils import rdkit_smiles_from_input_smiles
import pandas as pd
from umap import UMAP


df = pd.read_csv("processed_data/final_db_with_MHFP.csv")
df['RDKIT_SMILES'] = df['SMILES'].map(lambda x: rdkit_smiles_from_input_smiles(x))
df = df.drop(columns=['SMILES']).set_index(["CID", "RDKIT_SMILES", "Type"])

scaler = StandardScaler()
np_scaled = scaler.fit_transform(df)

umap_model = UMAP(n_neighbors=50, min_dist=0.25, random_state=42)
np_umap_result = umap_model.fit_transform(np_scaled)
df_umap_result = pd.DataFrame(np_umap_result, columns=['UMAP-1', 'UMAP-2'])
df = df.reset_index()
df = df[['CID', 'RDKIT_SMILES', 'Type']]
df_concat = pd.concat([df, df_umap_result], axis=1)
df_concat.to_csv("processed_data/final_db_with_UMAP.csv", index=False)

