import pandas as pd
from rdkit import Chem
from mhfp.encoder import MHFPEncoder


# calculate mhfp descriptors
df = pd.read_csv("input_data/final_perovskite_amines_database.csv", usecols=['CID', 'SMILES', 'Type'])

smiles_list = df['SMILES'].tolist()
mols = [Chem.MolFromSmiles(smiles) for smiles in smiles_list]

mhfp_encoder = MHFPEncoder()
mfhps = [mhfp_encoder.encode_mol(mol) for mol in mols]
df_mfhps = pd.DataFrame(mfhps)
df_concat = pd.concat([df, df_mfhps], axis=1)
df_concat.to_csv("processed_data/final_db_with_MHFP.csv", index=False)



