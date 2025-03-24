
from rdkit import Chem
import itertools

def cano_smiles(smi,chirality=True):
    try:
        return Chem.MolToSmiles(Chem.MolFromSmiles(smi), canonical=True, kekuleSmiles=True,isomericSmiles=chirality)
    except Exception as e:
        print(smi)
        return smi
def remove_charges(mol):
    # 遍历分子中的每个原子，设置它们的电荷为0
    for atom in mol.GetAtoms():
        charge = atom.GetFormalCharge()
        if charge != 0:
           atom.SetFormalCharge(0)
    return mol

def generate_chiral_isomers(mol):
    # 查找手性不确定的原子编号，并生成对映异构体
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    chiral_atom_indices = [atom_idx for atom_idx, chirality in chiral_centers if chirality == '?']
    isomers = []
    chiral_RS = [Chem.ChiralType.CHI_TETRAHEDRAL_CW, Chem.ChiralType.CHI_TETRAHEDRAL_CCW]
    chirality_combinations = itertools.product(chiral_RS, repeat=len(chiral_atom_indices))
    # 对每个手性中心尝试改变手性配置 (R/S 反转)
    for chirality_comb in chirality_combinations:
        mol_copy = Chem.Mol(mol)  
        for idx, chirality_tag in zip(chiral_atom_indices, chirality_comb):
            atom = mol_copy.GetAtomWithIdx(idx)
            atom.SetChiralTag(chirality_tag)  # 设置手性配置
        isomers.append(mol_copy)
    return isomers

