from sklearn import preprocessing
from sklearn.decomposition import PCA
import pandas as pd


def merge_smiles_and_descriptors(smiles_path, name_column, descriptors_path, column_drop=True, name_drop=True,
                                 gmin_drop=True, rdkit_smiles=True):
    if rdkit_smiles:
        df_smiles = pd.read_csv(smiles_path, usecols=["RDKIT_SMILES", name_column], low_memory=False)
    else:
        df_smiles = pd.read_csv(smiles_path, usecols=["SMILES", name_column], low_memory=False)
    df_descriptors = pd.read_csv(descriptors_path, na_values=["Infinity", "-Infinity"], low_memory=False)
    df_descriptors = df_descriptors.drop(["Name"], axis=1)
    df_merge = pd.concat([df_smiles.reset_index(drop=True), df_descriptors.reset_index(drop=True)], axis=1)
    if name_drop:
        df_merge = df_merge.drop([name_column], axis=1)
    if column_drop:
        df_merge = df_merge.dropna(axis=1)
    if gmin_drop:
        df_merge = df_merge.drop(columns=["gmin"])  # for ammonium NH3
    return df_merge


def pca_onestep(df):
    scaler = preprocessing.StandardScaler()
    df_scaled = pd.DataFrame(scaler.fit_transform(df), columns=df.columns)
    pca = PCA(n_components=3)
    pca_result = pca.fit_transform(df_scaled)
    scores = pd.DataFrame(pca_result, index=df.index, columns=['PC-1', 'PC-2', 'PC-3'])
    return scores

