from rdkit import Chem
from rdkit.Chem import rdmolops

import pandas as pd
import os,argparse
from  peptide2AA import peptide_to_AA
from read_Char import read_known_symbol,read_defined_symbol
from name_AA import Name_new_AA
from modify_smi import remove_charges,cano_smiles,generate_chiral_isomers

def Smiles2Sequence(mol,symbol_smi_dict,known_residuetypes,output_dir,file_name):
    peptide_sequence = []
    known_char = known_residuetypes + list(symbol_smi_dict.values())    #多肽切成单个氨基酸
    fragments =  peptide_to_AA(mol)
    fragments_neutral = [remove_charges(mol) for mol in fragments ]
    #替换dummy原子,查看是否在已知氨基酸里
    smis = [Chem.MolToSmiles(frag) for frag in fragments_neutral ]
    smis = [smi.replace('*N','N').replace('C(*)','C(O)').replace('*C','OC').replace('N1*','N1').replace('N(*)','N') for smi in smis]
    smis = [cano_smiles(smi) for smi in smis ]
    for smi in smis:
        if smi not in symbol_smi_dict:
            residue_name = Name_new_AA(known_char)
            print(f'新的氨基酸：{residue_name}, Smiles是: {smi}')
            symbol_smi_dict[smi] = residue_name
        else:
            residue_name = symbol_smi_dict[smi]
        peptide_sequence.append(residue_name)
    print("-".join(peptide_sequence),file_name)
    df = pd.DataFrame(list(symbol_smi_dict.items()), columns=['SMILES', 'Symbol'])
    df.to_csv('data/symbol.csv', index=False)
    return symbol_smi_dict

def main():
    parser = argparse.ArgumentParser(description="多肽sdf转PDB工具")
    parser.add_argument('-i', '--csv_file', required=True, help="输入csv文件名")
    parser.add_argument('-o', '--output_path', default='Output', help="输出路径")
    parser.add_argument('-smi', '--smi', default='smiles', help="输入smiles的title")
    parser.add_argument('-name', '--name', default='Example', help="输入分子name的title")
    args = parser.parse_args()
    file_name = args.csv_file
    name_title = args.name
    smiles_title =args.smi
    # 读取已定义的氨基酸及其SMILES
    script_directory = os.path.dirname(os.path.abspath(__file__))
    output_dir = os.path.join(os.getcwd(),args.output_path)
    df_patent = pd.read_csv(f'{file_name}.csv')
    if not os.path.exists(args.output_path):
        os.makedirs(args.output_path)
    os.chdir(script_directory)
    symbol_smi_dict = read_defined_symbol()
    # 已经用过的三字母氨基酸
    known_residuetypes = read_known_symbol()
    for mol_name,smi in zip(df_patent[name_title],df_patent[smiles_title]):
        mol = Chem.MolFromSmiles(smi)
        mols = generate_chiral_isomers(mol)
        for num,mol in enumerate(mols):
            if num:
                mol_name = f'{mol_name}-{num}'
            symbol_smi_dict = Smiles2Sequence(mol, symbol_smi_dict, known_residuetypes, output_dir, mol_name)

if __name__ == "__main__":
    main()