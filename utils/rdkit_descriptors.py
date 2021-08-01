from rdkit.Chem import Descriptors, MolFromSmiles


def count_heavy_atoms(smiles):
    m = MolFromSmiles(smiles)
    try:
        return Descriptors.HeavyAtomCount(m)
    except:
        return -1


def get_rdkit_descriptors(df):
    df['heavy_atom_count'] = df.apply(lambda x: count_heavy_atoms(x['RDKIT_SMILES']), axis=1)
    return df

