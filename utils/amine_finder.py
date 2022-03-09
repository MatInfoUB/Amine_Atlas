import pandas as pd
from rdkit.Chem import Fragments, PandasTools, MolFromSmiles
import numpy as np


def find_fr_Ar_N(x):  # Number of aromatic nitrogens
    try:
        return Fragments.fr_Ar_N(x['mol'])
    except:
        return -1


def find_fr_NH0(x):  # Number of Tertiary amines
    try:
        return Fragments.fr_NH0(x['mol'])
    except:
        return -1


def find_fr_NH1(x):  # Number of Secondary amines
    try:
        return Fragments.fr_NH1(x['mol'])
    except:
        return -1


def find_fr_NH2(x):  # Number of Primary amines
    try:
        return Fragments.fr_NH2(x['mol'])
    except:
        return -1


def find_fr_aniline(x):  # Number of anilines
    try:
        return Fragments.fr_aniline(x['mol'])
    except:
        return -1


def find_fr_benzene(x): #  Number of benzene rings
    try:
        return Fragments.fr_benzene(x['mol'])
    except:
        return -1


def find_fr_piperzine(x): # Number of piperzine rings
    try:
        return Fragments.fr_piperzine(x['mol'])
    except:
        return -1


def find_fr_thiazole(x): # Number of thiazole rings
    try:
        return Fragments.fr_thiazole(x['mol'])
    except:
        return -1


def find_fr_thiophene(x): # Number of thiophene rings
    try:
        return Fragments.fr_thiophene(x['mol'])
    except:
        return -1


def find_all_amine_types_simple(df_origin):
    PandasTools.AddMoleculeColumnToFrame(df_origin, smilesCol='SMILES', molCol='mol', includeFingerprints=True)
    df_origin['fr_NH0'] = df_origin.apply(lambda x: find_fr_NH0(x), axis=1)
    df_origin['fr_NH1'] = df_origin.apply(lambda x: find_fr_NH1(x), axis=1)
    df_origin['fr_NH2'] = df_origin.apply(lambda x: find_fr_NH2(x), axis=1)
    df_origin['fr_all_N'] = df_origin['fr_NH0'] + df_origin['fr_NH1'] + df_origin['fr_NH2']
    df_origin = df_origin.drop(columns=["mol"])
    return df_origin


def find_benzyl(x):   # newly added on 07/28/2021
    patt = MolFromSmiles('Cc1ccccc1')
    if x['mol'].HasSubstructMatch(patt):
        return 1
    return 0


def find_2phenylpropyl(x):  # newly added on 07/29/2021
    patt = MolFromSmiles('CC(CN)c1ccccc1')
    if x['mol'].HasSubstructMatch(patt):
        return 1
    return 0


def find_1phenylethyl(x):   # newly added on 07/29/2021
    patt = MolFromSmiles('CC(N)c1ccccc1')
    if x['mol'].HasSubstructMatch(patt):
        return 1
    return 0


def find_amphetamine(x):   # newly added on 07/29/2021
    patt = MolFromSmiles('CC(N)Cc1ccccc1')
    if x['mol'].HasSubstructMatch(patt):
        return 1
    return 0


def find_n_alkylbenzyl(x):   # newly added on 07/29/2021
    patt = MolFromSmiles('CNCc1ccccc1')
    if x['mol'].HasSubstructMatch(patt):
        return 1
    return 0


def find_n_methylphenylethyl(x):   # newly added on 07/29/2021
    patt = MolFromSmiles('CNCCc1ccccc1')
    if x['mol'].HasSubstructMatch(patt):
        return 1
    return 0


def find_aromatic_piperazine(x):   # newly added on 07/29/2021
    try:
        return Fragments.fr_piperzine(x['mol'])
    except:
        return -1


def find_phenylmethyl(x):   # newly added on 07/29/2021
    patt = MolFromSmiles('NCc1ccccc1')
    if x['mol'].HasSubstructMatch(patt):
        return 1
    return 0


def find_phenylethyl(x):  # newly added on 07/29/2021
    patt = MolFromSmiles('NCCc1ccccc1')
    if x['mol'].HasSubstructMatch(patt):
        return 1
    return 0


def find_phenylpropyl(x):   # newly added on 07/29/2021
    patt = MolFromSmiles('NCCCc1ccccc1')
    if x['mol'].HasSubstructMatch(patt):
        return 1
    return 0


def find_phenylbutyl(x):   # newly added on 07/29/2021
    patt = MolFromSmiles('NCCCCc1ccccc1')
    if x['mol'].HasSubstructMatch(patt):
        return 1
    return 0


# 2022.03.08，modify classification of hetero aromatic amine
def find_imidazole(x):
    patt = MolFromSmiles('C1=CN=CN1')
    if x['mol'].HasSubstructMatch(patt):
        return 1
    return 0


def find_indole(x):
    patt = MolFromSmiles('C1=CC=C2C(=C1)C=CN2')
    if x['mol'].HasSubstructMatch(patt):
        return 1
    return 0


def find_triazole(x):
    patt1 = MolFromSmiles('C1=NNN=C1')
    patt2 = MolFromSmiles('C1=NC=NN1')
    if x['mol'].HasSubstructMatch(patt1) or x['mol'].HasSubstructMatch(patt2):
        return 1
    return 0


def find_thiophene(x):
    patt = MolFromSmiles('C1=CSC=C1')
    if x['mol'].HasSubstructMatch(patt):
        return 1
    return 0


def find_pyridine(x):
    patt = MolFromSmiles('C1=CC=NC=C1')
    if x['mol'].HasSubstructMatch(patt):
        return 1
    return 0


def find_rdkit_fragments(df):
    PandasTools.AddMoleculeColumnToFrame(df, smilesCol='RDKIT_SMILES', molCol='mol', includeFingerprints=False)
    df['fr_Ar_N'] = df.apply(lambda x: find_fr_Ar_N(x), axis=1)
    df['fr_NH0'] = df.apply(lambda x: find_fr_NH0(x), axis=1)
    df['fr_NH1'] = df.apply(lambda x: find_fr_NH1(x), axis=1)
    df['fr_NH2'] = df.apply(lambda x: find_fr_NH2(x), axis=1)
    df['fr_aniline'] = df.apply(lambda x: find_fr_aniline(x), axis=1)
    df['fr_benzene'] = df.apply(lambda x: find_fr_benzene(x), axis=1)
    df['fr_thiazole'] = df.apply(lambda x: find_fr_thiazole(x), axis=1)
    df['fr_thiophene'] = df.apply(lambda x: find_fr_thiophene(x), axis=1)
    df = df.drop(columns=["mol"])
    return df


def find_benzyl_amines(df):
    PandasTools.AddMoleculeColumnToFrame(df, smilesCol='RDKIT_SMILES', molCol='mol', includeFingerprints=False)
    df['2phenylpropyl'] = df.apply(lambda x: find_2phenylpropyl(x), axis=1)
    df['1phenylethyl'] = df.apply(lambda x: find_1phenylethyl(x), axis=1)
    df['amphetamine'] = df.apply(lambda x: find_amphetamine(x), axis=1)
    df['n_alkylbenzyl'] = df.apply(lambda x: find_n_alkylbenzyl(x), axis=1)
    df['n_methylphenylethyl'] = df.apply(lambda x: find_n_methylphenylethyl(x), axis=1)
    df['aromatic_piperazine'] = df.apply(lambda x: find_aromatic_piperazine(x), axis=1)
    df['phenylmethyl'] = df.apply(lambda x: find_phenylmethyl(x), axis=1)
    df['phenylethyl'] = df.apply(lambda x: find_phenylethyl(x), axis=1)
    df['phenylpropyl'] = df.apply(lambda x: find_phenylpropyl(x), axis=1)
    df['phenylbutyl'] = df.apply(lambda x: find_phenylbutyl(x), axis=1)
    df = df.drop(columns=["mol"])
    return df


def find_heterocyclic(df):
    PandasTools.AddMoleculeColumnToFrame(df, smilesCol='RDKIT_SMILES', molCol='mol', includeFingerprints=False)
    df['imidazole'] = df.apply(lambda x: find_imidazole(x), axis=1)
    df['indole'] = df.apply(lambda x: find_indole(x), axis=1)
    df['triazole'] = df.apply(lambda x: find_triazole(x), axis=1)
    df['thiophene'] = df.apply(lambda x: find_thiophene(x), axis=1)
    df['pyridine'] = df.apply(lambda x: find_pyridine(x), axis=1)
    df = df.drop(columns=["mol"])
    return df




