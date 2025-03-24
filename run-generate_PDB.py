from rdkit import Chem
from rdkit.Chem import rdmolops

import pandas as pd
import os,argparse
from  peptide2AA import peptide_to_AA
from read_Char import read_known_symbol,read_defined_symbol
from name_AA import Name_new_AA
from modify_PDB import fragment_to_PDB
from modify_smi import remove_charges,cano_smiles
from connect_peptide_PDB import connect_AA_to_peptide

def generate_peptide_PDB(mol,known_char,symbol_smi_dict,output_dir,file_name):
    peptide_sequence = []
    #多肽切成单个氨基酸
    fragments =  peptide_to_AA(mol)
    fragments_neutral = [remove_charges(mol) for mol in fragments ]
    #替换dummy原子，并对分子消旋,查看是否在已知氨基酸里
    smis = [cano_smiles(Chem.MolToSmiles(frag)) for frag in fragments_neutral ]
    AA_index_list,pdb_lines = [1],[]
    for smi,frag in zip(smis,fragments):
        if smi not in symbol_smi_dict:
            residue_name = Name_new_AA(known_char)
            print(f'新的氨基酸：{residue_name}, Smiles是: {smi}')
            symbol_smi_dict[smi] = residue_name
        else:
            residue_name = symbol_smi_dict[smi]
        peptide_sequence.append(residue_name)
        #将片段保存为PDB
        new_pdb_line,AA_index_list = fragment_to_PDB(frag,residue_name,AA_index_list)
        pdb_lines += new_pdb_line
    connect_AA_to_peptide(pdb_lines,file_name,output_dir,AA_index_list[:-1])
    return symbol_smi_dict,peptide_sequence

def smiles_to_PDB(mol,symbol_smi_dict,known_residuetypes,output_dir,name):
    #生成不同的手性
    known_char = set(known_residuetypes + list(symbol_smi_dict.values()) )
    symbol_smi_dict,peptide_sequence = generate_peptide_PDB(mol,known_char,symbol_smi_dict,output_dir,name)
    print(f'氨基酸序列是:{'-'.join(peptide_sequence)}')
    #更新新的symbol,smiles列表
    df = pd.DataFrame(list(symbol_smi_dict.items()), columns=['SMILES', 'Symbol'])
    df.to_csv('data/symbol.csv', index=False)
    return symbol_smi_dict

def main():
    parser = argparse.ArgumentParser(description="多肽sdf转PDB工具")
    parser.add_argument('-i', '--input_sdf', required=True, help="输入SDF文件名")
    parser.add_argument('-o', '--output_path', default='Output', help="输出路径")
    args = parser.parse_args()
    # 读取已定义的氨基酸及其SMILES
    script_directory = os.path.dirname(os.path.abspath(__file__))
    output_dir = os.path.join(os.getcwd(),args.output_path)
    file_name = args.input_sdf
    supplier = Chem.SDMolSupplier(f'{file_name}.sdf')
    if not os.path.exists(args.output_path):
        os.makedirs(args.output_path)
    os.chdir(script_directory)
    symbol_smi_dict = read_defined_symbol()
    # 已经用过的三字母氨基酸
    known_residuetypes = read_known_symbol()
    
    print(len(supplier))
    for mol in supplier:
        #已有3d坐标，所以不需要手性信息
        rdmolops.RemoveStereochemistry(mol)
        mol = Chem.RemoveHs(mol)
        mol_name = mol.GetProp('_Name') 
        smiles_to_PDB(mol, symbol_smi_dict, known_residuetypes, output_dir, mol_name)

if __name__ == "__main__":
    main()