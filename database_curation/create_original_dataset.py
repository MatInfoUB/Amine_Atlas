from utils import get_compounds_from_cids
import pandas as pd


df_amine = pd.read_csv("data/processed_original_set.csv")

df_amine.fillna(-9999, inplace=True)
df_amine['CID'] = df_amine['CID'].astype(int)
print(df_amine.head())
df_amine_output = get_compounds_from_cids(list(df_amine["CID"]))
df_to_concat = df_amine_output[['SMILES', 'AIDs']]
print(df_to_concat.head())
df_cat = pd.concat([df_amine, df_to_concat], axis=1)
print(df_cat.head(30))
df_cat.to_csv("data/original_set_with_SMILES_and_AIDs.csv", index=False)

