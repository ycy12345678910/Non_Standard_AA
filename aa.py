from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd
from modify_smi import cano_smiles

def convert_to_L_amino_acid(smi):
    smi = smi.replace('*N','N').replace('C(*)','C(O)').replace('*C','OC').replace('N1*','N1').replace('N(*)','N')
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        raise ValueError("无效的 SMILES 表达式")
    # 赋予手性信息
    mol = Chem.AddHs(mol)
    AllChem.AssignAtomChiralTagsFromStructure(mol)
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    chiral_atom_indices = [atom_idx for atom_idx, chirality in chiral_centers if chirality == '?']
    print(chiral_atom_indices)
    # 遍历所有原子，确保 α-碳 是 L-型
    for atom_idx in chiral_atom_indices:
        atom = mol.GetAtomWithIdx(atom_idx)
        Neighbors = [n.GetSymbol() for n in atom.GetNeighbors()]
        if 'C' in Neighbors and 'N' in Neighbors:
            print([n.GetSymbol() for n in atom.GetNeighbors()])
            # 半胱氨酸 (Cys) 例外: 它的 L-型是 (R) 配置
            if 'S' in Neighbors:
                atom.SetChiralTag(Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW)
            else:
                atom.SetChiralTag(Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW)

    # 移式氢
    mol = Chem.RemoveHs(mol)
    # 生成 L-型 SMILES
    smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
    if chiral_atom_indices:
        print(smiles,chiral_atom_indices,'lll')
        cano_smi = cano_smiles(smiles)
        return cano_smi

df = pd.read_csv('data/symbol-new.csv')
df['D'] = df['SMILES'].apply(convert_to_L_amino_acid)
 
print(df)
df.to_csv('data/symbol-new.csv',index= False)

