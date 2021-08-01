import pandas as pd
from utils import mol_from_smiles, rdkit_smiles_from_mol


class Amine:

    def __init__(self, input_smiles, df_descriptors):

        self.mol = mol_from_smiles(input_smiles)
        if self.mol is not None:
            self.rdkit_smiles = rdkit_smiles_from_mol(self.mol)
            self.descriptors = fetch_padel_descriptors(self.rdkit_smiles, df_descriptors)
        else:
            self.rdkit_smiles = ""
            self.descriptors = {}

    def get_mol(self):
        return self.mol

    def get_rdkit_smiles(self):
        return self.rdkit_smiles

    def get_descriptors(self):
        return self.descriptors

    def get_nc(self):
        if len(self.descriptors) != 0:
            return int(float(self.descriptors["nC"]))
        return -1

    def get_nf(self):
        if len(self.descriptors) != 0:
            return int(float(self.descriptors["nF"]))
        return -1

    def get_no(self):
        if len(self.descriptors) != 0:
            return int(float(self.descriptors["nO"]))
        return -1

    def get_nh(self):
        if len(self.descriptors) != 0:
            return int(float(self.descriptors["nH"]))
        return -1

    def get_nn(self):
        if len(self.descriptors) != 0:
            return int(float(self.descriptors["nN"]))
        return -1

    def get_ns(self):
        if len(self.descriptors) != 0:
            return int(float(self.descriptors["nS"]))
        return -1

    def get_np(self):
        if len(self.descriptors) != 0:
            return int(float(self.descriptors["nP"]))
        return -1

    def get_ni(self):
        if len(self.descriptors) != 0:
            return int(float(self.descriptors["nI"]))
        return -1

    def get_atom_num(self):
        if len(self.descriptors) != 0:
            return int(float(self.descriptors["nAtom"]))
        return -1

    def get_heavy_atom_num(self):
        if len(self.descriptors) != 0:
            return int(float(self.descriptors["nHeavyAtom"]))
        return -1

    def get_aromatic_bond_num(self):
        if len(self.descriptors) != 0:
            return int(float(self.descriptors["nAromBond"]))
        return -1

    def get_ring_num(self):
        if len(self.descriptors) != 0:
            return int(float(self.descriptors["nRing"]))
        return -1

    def get_hetero_ring_num(self):
        if len(self.descriptors) != 0:
            return int(float(self.descriptors["nHeteroRing"]))
        return -1

    def get_longest_aliphatic_atom_num(self):
        if len(self.descriptors) != 0:
            return int(float(self.descriptors["nAtomLAC"]))
        return -1

    def get_longest_chain_atom_num(self):
        if len(self.descriptors) != 0:
            return int(float(self.descriptors["nAtomLC"]))
        return -1


def fetch_padel_descriptors(smiles, df):
    df_descriptors = df[df['RDKIT_SMILES'] == smiles].copy()
    # df_descriptors.drop(columns=['SMILES', 'RDKIT_SMILES'], inplace=True)
    df_descriptors.drop(columns=['RDKIT_SMILES'], inplace=True)
    return df_descriptors.to_dict('records')[0]
