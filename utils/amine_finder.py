import pandas as pd
from rdkit.Chem import Fragments, PandasTools, MolFromSmiles
import numpy as np


# def find_fr_Al_COO(x):  # Number of aliphatic carboxylic acids
#     try:
#         return Fragments.fr_Al_COO(x['mol'])
#     except:
#         return -1
#
#
# def find_fr_Al_OH(x):  # Number of aliphatic hydroxyl groups
#     try:
#         return Fragments.fr_Al_OH(x['mol'])
#     except:
#         return -1
#
#
# def find_fr_Al_OH_noTert(x):  # Number of aliphatic hydroxyl groups excluding tert-OH
#     try:
#         return Fragments.fr_Al_OH_noTert(x['mol'])
#     except:
#         return -1


# def find_fr_ArN(x):  # Number of N functional groups attached to aromatics
#     try:
#         return Fragments.fr_ArN(x['mol'])
#     except:
#         return -1


# def find_fr_Ar_COO(x):  # Number of Aromatic carboxylic acid
#     try:
#         return Fragments.fr_Ar_COO(x['mol'])
#     except:
#         return -1
#
#
def find_fr_Ar_N(x):  # Number of aromatic nitrogens
    try:
        return Fragments.fr_Ar_N(x['mol'])
    except:
        return -1

#
# def find_fr_Ar_NH(x):  # Number of aromatic amines
#     try:
#         return Fragments.fr_Ar_NH(x['mol'])
#     except:
#         return -1
#
#
# def find_fr_Ar_OH(x):  # Number of aromatic hydroxyl groups
#     try:
#         return Fragments.fr_Ar_OH(x['mol'])
#     except:
#         return -1
#
#
# def find_fr_COO(x):  # Number of carboxylic acids
#     try:
#         return Fragments.fr_COO(x['mol'])
#     except:
#         return -1
#
#
# def find_fr_COO2(x):  # Number of carboxylic acids
#     try:
#         return Fragments.fr_COO2(x['mol'])
#     except:
#         return -1
#
#
# def find_fr_C_O(x):  # Number of carbonyl O
#     try:
#         return Fragments.fr_C_O(x['mol'])
#     except:
#         return -1
#
#
# def find_fr_C_O_no_COO(x):  # Number of carbonyl O, excluding COOH
#     try:
#         return Fragments.fr_C_O_noCOO(x['mol'])
#     except:
#         return -1
#
#
# def find_fr_C_S(x):  # Number of thiocarbonyl
#     try:
#         return Fragments.fr_C_S(x['mol'])
#     except:
#         return -1
#
#
# def find_fr_HOCCN(x):  # Number of C(OH)CCN-Ctert-alkyl or C(OH)CCNcyclic
#     try:
#         return Fragments.fr_HOCCN(x['mol'])
#     except:
#         return -1
#
#
# def find_fr_Imine(x):  # Number of Imines
#     try:
#         return Fragments.fr_Imine(x['mol'])
#     except:
#         return -1
#
#
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
#
#
# def find_fr_N_O(x):  # Number of hydroxylamine groups
#     try:
#         return Fragments.fr_N_O(x['mol'])
#     except:
#         return -1
#
#
# def find_fr_Ndealkylation1(x):  # Number of XCCNR groups
#     try:
#         return Fragments.fr_Ndealkylation1(x['mol'])
#     except:
#         return -1
#
#
# def find_fr_Ndealkylation2(x):  # Number of tert-alicyclic amines (no heteroatoms, not quinine-like bridged N)
#     try:
#         return Fragments.fr_Ndealkylation2(x['mol'])
#     except:
#         return -1
#
#
# def find_fr_Nhpyrrole(x):  # Number of H-pyrrole nitrogens
#     try:
#         return Fragments.fr_Nhpyrrole(x['mol'])
#     except:
#         return -1
#
#
# def find_fr_SH(x):  # Number of thiol groups
#     try:
#         return Fragments.fr_SH(x['mol'])
#     except:
#         return -1
#
#
# def find_fr_aldehyde(x):  # Number of aldehydes
#     try:
#         return Fragments.fr_aldehyde(x['mol'])
#     except:
#         return -1
#
#
# def find_fr_alkyl_carbamate(x):  # Number of alkyl carbamates (subject to hydrolysis)
#     try:
#         return Fragments.fr_alkyl_carbamate(x['mol'])
#     except:
#         return -1
#
#
# def find_fr_alkyl_halide(x):  # Number of alkyl halides
#     try:
#         return Fragments.fr_alkyl_halide(x['mol'])
#     except:
#         return -1
#
#
# def find_fr_allylic_oxid(x):  # Number of allylic oxidation sites excluding steroid dienone
#     try:
#         return Fragments.fr_allylic_oxid(x['mol'])
#     except:
#         return -1
#
#
# def find_fr_amide(x):  # Number of amides
#     try:
#         return Fragments.fr_amide(x['mol'])
#     except:
#         return -1
#
#
# def find_fr_amidine(x):  # Number of amidine groups
#     try:
#         return Fragments.fr_amidine(x['mol'])
#     except:
#         return -1


def find_fr_aniline(x):  # Number of anilines
    try:
        return Fragments.fr_aniline(x['mol'])
    except:
        return -1


# def find_fr_aryl_methyl(x):  # Number of aryl methyl sites for hydroxylation
#     try:
#         return Fragments.fr_aryl_methyl(x['mol'])
#     except:
#         return -1
#
#
# def find_fr_azide(x):  # Number of azide groups
#     try:
#         return Fragments.fr_azide(x['mol'])
#     except:
#         return -1
#
#
# def find_fr_azo(x):  # Number of azo groups
#     try:
#         return Fragments.fr_azo(x['mol'])
#     except:
#         return -1
#
#
# def find_fr_barbitur(x):  #  Number of barbiturate groups
#     try:
#         return Fragments.fr_barbitur(x['mol'])
#     except:
#         return -1


def find_fr_benzene(x): #  Number of benzene rings
    try:
        return Fragments.fr_benzene(x['mol'])
    except:
        return -1


# def find_fr_benzodiazepine(x):  # Number of benzodiazepines with no additional fused rings
#     try:
#         return Fragments.fr_benzodiazepine(x['mol'])
#     except:
#         return -1
#
#
# def find_fr_bicyclic(x): # Bicyclic
#     try:
#         return Fragments.fr_bicyclic(x['mol'])
#     except:
#         return -1
#
#
# def find_fr_diazo(x): #  Number of diazo groups
#     try:
#         return Fragments.fr_diazo(x['mol'])
#     except:
#         return -1
#
#
# def find_fr_dihydropyridine(x):  # Number of dihydropyridines
#     try:
#         return Fragments.fr_dihydropyridine(x['mol'])
#     except:
#         return -1
#
#
# def find_fr_epoxide(x): #  Number of epoxide rings
#     try:
#         return Fragments.fr_epoxide(x['mol'])
#     except:
#         return -1
#
#
# def find_fr_ester(x): # Number of esters
#     try:
#         return Fragments.fr_ester(x['mol'])
#     except:
#         return -1
#
#
# def find_fr_ether(x): # Number of ether oxygens (including phenoxy)
#     try:
#         return Fragments.fr_ether(x['mol'])
#     except:
#         return -1
#
#
# def find_fr_furan(x): # Number of furan rings
#     try:
#         return Fragments.fr_furan(x['mol'])
#     except:
#         return -1
#
#
# def find_fr_guanido(x): # Number of guanidine groups
#     try:
#         return Fragments.fr_guanido(x['mol'])
#     except:
#         return -1
#
#
# def find_fr_halogen(x): # Number of halogens
#     try:
#         return Fragments.fr_halogen(x['mol'])
#     except:
#         return -1
#
#
# def find_fr_hdrzine(x): # Number of hydrazine groups
#     try:
#         return Fragments.fr_hdrzine(x['mol'])
#     except:
#         return -1
#
#
# def find_fr_hdrzone(x): # Number of hydrazone groups
#     try:
#         return Fragments.fr_hdrzone(x['mol'])
#     except:
#         return -1
#
#
# def find_fr_imidazole(x): # Number of imidazole rings
#     try:
#         return Fragments.fr_imidazole(x['mol'])
#     except:
#         return -1
#
#
# def find_fr_imide(x): # Number of imide groups
#     try:
#         return Fragments.fr_imide(x['mol'])
#     except:
#         return -1
#
#
# def find_fr_isocyan(x): # Number of isocyanates
#     try:
#         return Fragments.fr_isocyan(x['mol'])
#     except:
#         return -1
#
#
# def find_fr_isothiocyan(x): # Number of isothiocyanates
#     try:
#         return Fragments.fr_isothiocyan(x['mol'])
#     except:
#         return -1
#
#
# def find_fr_ketone(x): # Number of ketones
#     try:
#         return Fragments.fr_ketone(x['mol'])
#     except:
#         return -1
#
#
# def find_fr_ketone_Topliss(x): # Number of ketones excluding diaryl, a,b-unsat. dienones, heteroatom on Calpha
#     try:
#         return Fragments.fr_ketone_Topliss(x['mol'])
#     except:
#         return -1
#
#
# def find_fr_lactam(x): # Number of beta lactams
#     try:
#         return Fragments.fr_lactam(x['mol'])
#     except:
#         return -1
#
#
# def find_fr_lactone(x): # Number of cyclic esters (lactones)
#     try:
#         return Fragments.fr_lactone(x['mol'])
#     except:
#         return -1
#
#
# def find_fr_methoxy(x): # Number of methoxy groups -OCH3
#     try:
#         return Fragments.fr_methoxy(x['mol'])
#     except:
#         return -1
#
#
# def find_fr_morpholine(x): # Number of morpholine rings
#     try:
#         return Fragments.fr_morpholine(x['mol'])
#     except:
#         return -1
#
#
# def find_fr_nitrile(x): # Number of nitriles
#     try:
#         return Fragments.fr_nitrile(x['mol'])
#     except:
#         return -1
#
#
# def find_fr_nitro(x): # Number of nitro groups
#     try:
#         return Fragments.fr_nitro(x['mol'])
#     except:
#         return -1
#
#
# def find_fr_nitro_arom(x): # Number of nitro benzene ring substituents
#     try:
#         return Fragments.fr_nitro_arom(x['mol'])
#     except:
#         return -1
#
#
# def find_fr_nitro_arom_nonortho(x): # Number of non-ortho nitro benzene ring substituents
#     try:
#         return Fragments.fr_nitro_arom_nonortho(x['mol'])
#     except:
#         return -1
#
#
# def find_fr_nitroso(x): # Number of nitroso groups, excluding NO2
#     try:
#         return Fragments.fr_nitroso(x['mol'])
#     except:
#         return -1
#
#
# def find_fr_oxazole(x): # Number of oxazole rings
#     try:
#         return Fragments.fr_oxazole(x['mol'])
#     except:
#         return -1
#
#
# def find_fr_oxime(x): # Number of oxime groups
#     try:
#         return Fragments.fr_oxime(x['mol'])
#     except:
#         return -1
#
#
# def find_fr_para_hydroxylation(x): # Number of para-hydroxylation sites
#     try:
#         return Fragments.fr_para_hydroxylation(x['mol'])
#     except:
#         return -1
#
#
# def find_fr_phenol(x): # Number of phenols
#     try:
#         return Fragments.fr_phenol(x['mol'])
#     except:
#         return -1
#
#
# def find_fr_phenol_noOrthoHbond(x): # Number of phenolic OH excluding ortho intramolecular Hbond substituents
#     try:
#         return Fragments.fr_phenol_noOrthoHbond(x['mol'])
#     except:
#         return -1
#
#
# def find_fr_phos_acid(x): # Number of phosphoric acid groups
#     try:
#         return Fragments.fr_phos_acid(x['mol'])
#     except:
#         return -1
#
#
# def find_fr_phos_ester(x): # Number of phosphoric ester groups
#     try:
#         return Fragments.fr_phos_ester(x['mol'])
#     except:
#         return -1
#
#
# def find_fr_piperdine(x): # Number of piperdine rings
#     try:
#         return Fragments.fr_piperdine(x['mol'])
#     except:
#         return -1
#
#
def find_fr_piperzine(x): # Number of piperzine rings
    try:
        return Fragments.fr_piperzine(x['mol'])
    except:
        return -1
#
#
# def find_fr_priamide(x): # Number of primary amides
#     try:
#         return Fragments.fr_priamide(x['mol'])
#     except:
#         return -1
#
#
# def find_fr_prisulfonamd(x): # Number of primary sulfonamides
#     try:
#         return Fragments.fr_prisulfonamd(x['mol'])
#     except:
#         return -1
#
#
# def find_fr_pyridine(x): # Number of pyridine rings
#     try:
#         return Fragments.fr_pyridine(x['mol'])
#     except:
#         return -1
#
#
# def find_fr_quatN(x): # Number of quarternary nitrogens
#     try:
#         return Fragments.fr_quatN(x['mol'])
#     except:
#         return -1
#
#
# def find_fr_sulfide(x): # Number of thioether
#     try:
#         return Fragments.fr_sulfide(x['mol'])
#     except:
#         return -1
#
#
# def find_fr_sulfonamd(x): # Number of sulfonamides
#     try:
#         return Fragments.fr_sulfonamd(x['mol'])
#     except:
#         return -1
#
#
# def find_fr_sulfone(x): # Number of sulfone groups
#     try:
#         return Fragments.fr_sulfone(x['mol'])
#     except:
#         return -1
#
#
# def find_fr_term_acetylene(x): # Number of terminal acetylenes
#     try:
#         return Fragments.fr_term_acetylene(x['mol'])
#     except:
#         return -1
#
#
# def find_fr_tetrazole(x): # Number of tetrazole rings
#     try:
#         return Fragments.fr_tetrazole(x['mol'])
#     except:
#         return -1
#
#
def find_fr_thiazole(x): # Number of thiazole rings
    try:
        return Fragments.fr_thiazole(x['mol'])
    except:
        return -1
#
#
# def find_fr_thiocyan(x): # Number of thiocyanates
#     try:
#         return Fragments.fr_thiocyan(x['mol'])
#     except:
#         return -1


def find_fr_thiophene(x): # Number of thiophene rings
    try:
        return Fragments.fr_thiophene(x['mol'])
    except:
        return -1


# def find_fr_unbrch_alkane(x): # Number of unbranched alkanes of at least 4 members (excludes halogenated alkanes)
#     try:
#         return Fragments.fr_unbrch_alkane(x['mol'])
#     except:
#         return -1
#
#
# def find_fr_urea(x): # Number of urea groups
#     try:
#         return Fragments.fr_urea(x['mol'])
#     except:
#         return -1


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


# def find_pyrimidine(x):
#     patt = MolFromSmiles('C1=CN=CN=C1')
#     if x['mol'].HasSubstructMatch(patt):
#         return 1
#     return 0


def find_rdkit_fragments(df):
    PandasTools.AddMoleculeColumnToFrame(df, smilesCol='RDKIT_SMILES', molCol='mol', includeFingerprints=False)
    # df['fr_Al_COO'] = df.apply(lambda x: find_fr_Al_COO(x), axis=1)
    # df['fr_Al_OH'] = df.apply(lambda x: find_fr_Al_OH(x), axis=1)
    # df['fr_Al_OH_noTert'] = df.apply(lambda x: find_fr_Al_OH_noTert(x), axis=1)
    # df['fr_ArN'] = df.apply(lambda x: find_fr_ArN(x), axis=1)
    # df['fr_Ar_COO'] = df.apply(lambda x: find_fr_Ar_COO(x), axis=1)
    df['fr_Ar_N'] = df.apply(lambda x: find_fr_Ar_N(x), axis=1)
    # df['fr_Ar_NH'] = df.apply(lambda x: find_fr_Ar_NH(x), axis=1)
    # df['fr_Ar_OH'] = df.apply(lambda x: find_fr_Ar_OH(x), axis=1)
    # df['fr_COO'] = df.apply(lambda x: find_fr_COO(x), axis=1)
    # df['fr_COO2'] = df.apply(lambda x: find_fr_COO2(x), axis=1)
    # df['fr_C_O'] = df.apply(lambda x: find_fr_C_O(x), axis=1)
    # df['fr_C_O_noCOO'] = df.apply(lambda x: find_fr_C_O_no_COO(x), axis=1)
    # df['fr_C_S'] = df.apply(lambda x: find_fr_C_S(x), axis=1)
    # df['fr_HOCCN'] = df.apply(lambda x: find_fr_HOCCN(x), axis=1)
    # df['fr_Imine'] = df.apply(lambda x: find_fr_Imine(x), axis=1)
    df['fr_NH0'] = df.apply(lambda x: find_fr_NH0(x), axis=1)
    df['fr_NH1'] = df.apply(lambda x: find_fr_NH1(x), axis=1)
    df['fr_NH2'] = df.apply(lambda x: find_fr_NH2(x), axis=1)
    # df['fr_N_O'] = df.apply(lambda x: find_fr_N_O(x), axis=1)
    # df['fr_Ndealkylation1'] = df.apply(lambda x: find_fr_Ndealkylation1(x), axis=1)
    # df['fr_Ndealkylation2'] = df.apply(lambda x: find_fr_Ndealkylation2(x), axis=1)
    # df['fr_Nhpyrrole'] = df.apply(lambda x: find_fr_Nhpyrrole(x), axis=1)
    # df['fr_SH'] = df.apply(lambda x: find_fr_SH(x), axis=1)
    # df['fr_aldehyde'] = df.apply(lambda x: find_fr_aldehyde(x), axis=1)
    # df['fr_alkyl_carbamate'] = df.apply(lambda x: find_fr_alkyl_carbamate(x), axis=1)
    # df['fr_alkyl_halide'] = df.apply(lambda x: find_fr_alkyl_halide(x), axis=1)
    # df['fr_allylic_oxid'] = df.apply(lambda x: find_fr_allylic_oxid(x), axis=1)
    # df['fr_amide'] = df.apply(lambda x: find_fr_amide(x), axis=1)
    # df['fr_amidine'] = df.apply(lambda x: find_fr_amidine(x), axis=1)
    df['fr_aniline'] = df.apply(lambda x: find_fr_aniline(x), axis=1)
    # df['fr_aryl_methyl'] = df.apply(lambda x: find_fr_aryl_methyl(x), axis=1)
    # df['fr_azide'] = df.apply(lambda x: find_fr_azide(x), axis=1)
    # df['fr_azo'] = df.apply(lambda x: find_fr_azo(x), axis=1)
    # df['fr_barbitur'] = df.apply(lambda x: find_fr_barbitur(x), axis=1)
    df['fr_benzene'] = df.apply(lambda x: find_fr_benzene(x), axis=1)
    # df['fr_benzodiazepine'] = df.apply(lambda x: find_fr_benzodiazepine(x), axis=1)
    # df['fr_bicyclic'] = df.apply(lambda x: find_fr_bicyclic(x), axis=1)
    # df['fr_diazo'] = df.apply(lambda x: find_fr_diazo(x), axis=1)
    # df['fr_dihydropyridine'] = df.apply(lambda x: find_fr_dihydropyridine(x), axis=1)
    # df['fr_epoxide'] = df.apply(lambda x: find_fr_epoxide(x), axis=1)
    # df['fr_ester'] = df.apply(lambda x: find_fr_ester(x), axis=1)
    # df['fr_ether'] = df.apply(lambda x: find_fr_ether(x), axis=1)
    # df['fr_furan'] = df.apply(lambda x: find_fr_furan(x), axis=1)
    # df['fr_guanido'] = df.apply(lambda x: find_fr_guanido(x), axis=1)
    # df['fr_halogen'] = df.apply(lambda x: find_fr_halogen(x), axis=1)
    # df['fr_hdrzine'] = df.apply(lambda x: find_fr_hdrzine(x), axis=1)
    # df['fr_hdrzone'] = df.apply(lambda x: find_fr_hdrzone(x), axis=1)
    # df['fr_imidazole'] = df.apply(lambda x: find_fr_imidazole(x), axis=1)
    # df['fr_imide'] = df.apply(lambda x: find_fr_imide(x), axis=1)
    # df['fr_isocyan'] = df.apply(lambda x: find_fr_isocyan(x), axis=1)
    # df['fr_isothiocyan'] = df.apply(lambda x: find_fr_isothiocyan(x), axis=1)
    # df['fr_ketone'] = df.apply(lambda x: find_fr_ketone(x), axis=1)
    # df['fr_ketone_Topliss'] = df.apply(lambda x: find_fr_ketone_Topliss(x), axis=1)
    # df['fr_lactam'] = df.apply(lambda x: find_fr_lactam(x), axis=1)
    # df['fr_lactone'] = df.apply(lambda x: find_fr_lactone(x), axis=1)
    # df['fr_methoxy'] = df.apply(lambda x: find_fr_methoxy(x), axis=1)
    # df['fr_morpholine'] = df.apply(lambda x: find_fr_morpholine(x), axis=1)
    # df['fr_nitrile'] = df.apply(lambda x: find_fr_nitrile(x), axis=1)
    # df['fr_nitro'] = df.apply(lambda x: find_fr_nitro(x), axis=1)
    # df['fr_nitro_arom'] = df.apply(lambda x: find_fr_nitro_arom(x), axis=1)
    # df['fr_nitro_arom_nonortho'] = df.apply(lambda x: find_fr_nitro_arom_nonortho(x), axis=1)
    # df['fr_nitroso'] = df.apply(lambda x: find_fr_nitroso(x), axis=1)
    # df['fr_oxazole'] = df.apply(lambda x: find_fr_oxazole(x), axis=1)
    # df['fr_oxime'] = df.apply(lambda x: find_fr_oxime(x), axis=1)
    # df['fr_para_hydroxylation'] = df.apply(lambda x: find_fr_para_hydroxylation(x), axis=1)
    # df['fr_phenol'] = df.apply(lambda x: find_fr_phenol(x), axis=1)
    # df['fr_phenol_noOrthoHbond'] = df.apply(lambda x: find_fr_phenol_noOrthoHbond(x), axis=1)
    # df['fr_phos_acid'] = df.apply(lambda x: find_fr_phos_acid(x), axis=1)
    # df['fr_phos_ester'] = df.apply(lambda x: find_fr_phos_ester(x), axis=1)
    # df['fr_piperdine'] = df.apply(lambda x: find_fr_piperdine(x), axis=1)
    # df['fr_piperzine'] = df.apply(lambda x: find_fr_piperzine(x), axis=1)
    # df['fr_priamide'] = df.apply(lambda x: find_fr_priamide(x), axis=1)
    # df['fr_prisulfonamd'] = df.apply(lambda x: find_fr_prisulfonamd(x), axis=1)
    # df['fr_pyridine'] = df.apply(lambda x: find_fr_pyridine(x), axis=1)
    # df['fr_quatN'] = df.apply(lambda x: find_fr_quatN(x), axis=1)
    # df['fr_sulfide'] = df.apply(lambda x: find_fr_sulfide(x), axis=1)
    # df['fr_sulfonamd'] = df.apply(lambda x: find_fr_sulfonamd(x), axis=1)
    # df['fr_sulfone'] = df.apply(lambda x: find_fr_sulfone(x), axis=1)
    # df['fr_term_acetylene'] = df.apply(lambda x: find_fr_term_acetylene(x), axis=1)
    # df['fr_tetrazole'] = df.apply(lambda x: find_fr_tetrazole(x), axis=1)
    df['fr_thiazole'] = df.apply(lambda x: find_fr_thiazole(x), axis=1)
    # df['fr_thiocyan'] = df.apply(lambda x: find_fr_thiocyan(x), axis=1)
    df['fr_thiophene'] = df.apply(lambda x: find_fr_thiophene(x), axis=1)
    # df['fr_unbrch_alkane'] = df.apply(lambda x: find_fr_unbrch_alkane(x), axis=1)
    # df['fr_urea'] = df.apply(lambda x: find_fr_urea(x), axis=1)
    # df['fr_benzyl'] = df.apply(lambda x: find_benzyl(x), axis=1)  # newly added on 07/28/2021
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


# def find_nitramine(x):
#     patt = MolFromSmiles('N[N+](=O)O')
#     if x['mol'].HasSubstructMatch(patt):
#         return 1
#     return 0


# def find_nitrosamine(df):
#     PandasTools.AddMoleculeColumnToFrame(df, smilesCol='SMILES', molCol='mol', includeFingerprints=False)
#     df['fr_NH0'] = df.apply(lambda x: find_fr_NH0(x), axis=1)
#     df['fr_NH1'] = df.apply(lambda x: find_fr_NH1(x), axis=1)
#     df['fr_NH2'] = df.apply(lambda x: find_fr_NH2(x), axis=1)
#     df['fr_all_N'] = df['fr_NH0'] + df['fr_NH1'] + df['fr_NH2']
#     # df['fr_nitroso'] = df.apply(lambda x: find_fr_nitroso(x), axis=1)
#     # df['nitramine'] = df.apply(lambda x: find_nitramine(x), axis=1)
#     df = df.drop(columns=["mol"])
#     return df



