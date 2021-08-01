from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import SaltRemover, MolToSmiles, MolFromSmiles, RemoveHs
import numpy as np


SALTS_FILE = "salts/Salts_from_epa_list.txt"

""" Source: http://www.mail-archive.com/rdkit-discuss@lists.sourceforge.net/msg02669.html """
""" contribution from Hans de Winter """


def _InitialiseNeutralisationReactions():
    patts= (
        # Imidazoles
        ('[n+;H]','n'),
        # Amines
        ('[N+;!H0]','N'),
        # Carboxylic acids and alcohols
        ('[$([O-]);!$([O-][#7])]','O'),
        # Thiols
        ('[S-;X1]','S'),
        # Sulfonamides
        ('[$([N-;X2]S(=O)=O)]','N'),
        # Enamines
        ('[$([N-;X2][C,N]=C)]','N'),
        # Tetrazoles
        ('[n-]','[nH]'),
        # Sulfoxides
        ('[$([S-]=O)]','S'),
        # Amides
        ('[$([N-]C=O)]','N'),
        )
    return [(Chem.MolFromSmarts(x),Chem.MolFromSmiles(y,False)) for x,y in patts]


_reactions=None


def NeutraliseCharges(smiles, reactions=None):
    global _reactions
    if reactions is None:
        if _reactions is None:
            _reactions=_InitialiseNeutralisationReactions()
        reactions=_reactions
    mol = Chem.MolFromSmiles(smiles)
    replaced = False
    for i,(reactant, product) in enumerate(reactions):
        while mol.HasSubstructMatch(reactant):
            replaced = True
            rms = AllChem.ReplaceSubstructs(mol, reactant, product)
            mol = rms[0]
    if replaced:
        return Chem.MolToSmiles(mol, True), True
    else:
        return smiles, False


"""Source: https://sourceforge.net/p/rdkit/mailman/message/36877847/"""


def MolWithoutIsotopesToSmiles(mol):
   atom_data = [(atom, atom.GetIsotope()) for atom in mol.GetAtoms()]
   for atom, isotope in atom_data:
       if isotope:
           atom.SetIsotope(0)
   return mol


"""Below written by me"""


def mol_from_smiles(smiles):
    try:
        mol = MolFromSmiles(smiles)
    except BaseException as error:
        print('An exception occurred: {}'.format(error))
        mol = None
    if mol is None:
        return None
    return MolWithoutIsotopesToSmiles(mol)


def rdkit_smiles_from_mol(mol):
    # res = SaltRemover.SaltRemover(defnFilename=SALTS_FILE).StripMol(mol)
    smiles = MolToSmiles(mol)
    smiles, _ = NeutraliseCharges(smiles)
    if "." in smiles:  # Check if structure after salt removal is chelate
        split_smiles = smiles.split(".", 5)
        if len(set(split_smiles)) == 1:
            return split_smiles[0]
    return smiles


def rdkit_smiles_from_input_smiles(smiles):
    mol = mol_from_smiles(smiles)
    if mol is not None:
        return rdkit_smiles_from_mol(mol)
    return np.nan