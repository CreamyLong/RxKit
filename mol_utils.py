import rdkit
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdChemReactions
from rdkit.Chem import MolStandardize
from tqdm import tqdm
from rdkit import RDLogger
import os
from rdkit import RDLogger

RDLogger.DisableLog('rdApp.*')


def smi_tokenizer(smi):
    """
    Tokenize a SMILES molecule or reaction
    """
    import re
    # pattern =  "(\[[^\]]+]|Br?|Cl?|N|O|S|P|F|I|b|c|n|o|s|p|\(|\)|\.|=|#|-|\+|\\\\|\/|:|~|@|\?|>|\*|\$|\%[0-9]{2}|[0-9])"
    pattern = "(\[[^\]]+]|Br?|Cl?|N|O|S|P|F|I|b|c|n|o|s|p|\(|\)|\.|=|#|-|\+|\\\\|\/|:|~|@|\?|>>?|\*|\$|\%\([0-9]{3}\)|\%[0-9]{2}|[0-9])"

    regex = re.compile(pattern)
    tokens = [token for token in regex.findall(smi)]
    assert smi == ''.join(tokens)
    return ' '.join(tokens)

def slf_tokenizer(smi):
    encoded_selfies = sf.encoder(smi)
    tokens = list(sf.split_selfies(encoded_selfies))
    assert encoded_selfies == ''.join(tokens)
    return ' '.join(tokens)

def canonicalize_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        [a.ClearProp('molAtomMapNumber') for a in mol.GetAtoms()]
        return Chem.MolToSmiles(mol, isomericSmiles=True)
    else:
        return ''

def augm_smile(smi):
    mol = Chem.MolFromSmiles(smi)
    smi = Chem.MolToSmiles(mol, doRandom=True)
    return smi

def inchify(smi, extended_tautomer_check=False):
    mol = Chem.MolFromSmiles(smi)
    options = ""
    if extended_tautomer_check:
        # https://pubs.acs.org/doi/10.1021/acs.jcim.9b01080
        options = "-KET -15T"
    inchi = Chem.MolToInchi(mol, options=options)
    mol = Chem.MolFromInchi(inchi)
    inchiy_smi = Chem.MolToSmiles(mol)
    if "." in inchiy_smi:
        #avoid break ion bond
#         print(" '.' in inchiy_smi: ",smi,"->",inchiy_smi)
        return smi
    else:
        return inchiy_smi

def getNumHeavyAtoms(smi):
    mol = Chem.MolFromSmiles(smi)
    return mol.GetNumHeavyAtoms()

def get_longest_smiles(smis):
    PPh3 = "c1ccc(cc1)P(c1ccccc1)c1ccccc1"
    smis_list = smis.split(".")
    if len(smis_list) == 1:
        return sorted(smis_list, key=lambda s: len(s), reverse=True)[0]
    elif len(smis_list) > 1:
        if PPh3 in smis_list:
            smis_list.remove(PPh3)  # avoid long agent
        return sorted(smis_list, key=lambda s: len(s), reverse=True)[0]
    else:
        return ""
