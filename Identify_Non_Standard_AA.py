from rdkit import Chem
import os
import pandas as pd

from  peptide2AA import peptide_to_AA,cano_smiles
from read_Char import read_known_symbol,read_defined_symbol
from name_AA import Name_new_AA
from modify_PDB import fragment_to_PDB

              
folder_path = "output"
smi = "NCCCC[C@H](NC([C@@H](NC([C@H](CCC(N)=O)N)=O)CCC(O)=O)=O)C(N[C@@H](CO)C(O)=O)=O" 
mol = Chem.MolFromSmiles(smi)
#切成单个氨基酸
fragments =  peptide_to_AA(mol)
#对分子消旋,查看是否在已知氨基酸里
smis = [cano_smiles(Chem.MolToSmiles(Chem.RemoveHs(frag)).replace('*N','N').replace('*',"O"),False) for frag in fragments]
#读取已经定义的氨基酸及其对应smiles
symbol_smi_dict = read_defined_symbol()
#已经用过的三字母
peptide_sequence = ''
known_char = set(read_known_symbol()+list(symbol_smi_dict.values()))
for smi,frag in zip(smis,fragments):
    if smi not in symbol_smi_dict:
        
        residue_name = Name_new_AA(known_char)
        print(residue_name,smi)
        symbol_smi_dict[smi] = residue_name
        if not os.path.exists(folder_path):
            os.makedirs(folder_path)
        file_path = os.path.join(folder_path, residue_name+  ".pdb")
        fragment_to_PDB(frag,residue_name,folder_path)
        peptide_sequence += residue_name + '-'
    else:
        peptide_sequence += symbol_smi_dict[smi] + '-'

print(peptide_sequence)
# #更新新的symbol,smiles列表
# df = pd.DataFrame(list(symbol_smi_dict.items()), columns=['SMILES', 'Symbol'])
# df.to_csv('data/symbol.csv', index=False)
