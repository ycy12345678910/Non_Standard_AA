import pandas as pd
from openbabel import pybel
from rdkit import Chem
from modify_smi import generate_chiral_isomers
import argparse

def main():
    parser = argparse.ArgumentParser(description="通过smiles生成3D构象")
    parser.add_argument('-i', '--csv_file', required=True, help="输入csv文件名")
    parser.add_argument('-smi', '--smi', default='smiles', help="输入smiles的title")
    parser.add_argument('-name', '--name', default='Example', help="输入分子name的title")
    args = parser.parse_args()
    file_name = args.csv_file
    name_title = args.name
    smiles_title =args.smi
    df_patent = pd.read_csv(f'{file_name}.csv')
    output_sdf_file = f'{file_name}.sdf'
    with pybel.Outputfile('sdf',output_sdf_file,overwrite=True) as outf:
        # 遍历CSV中的每一行
        for mol_name,smi in zip(df_patent[name_title],df_patent[smiles_title]):
            mol = Chem.MolFromSmiles(smi)
            #当手性不确定时，生成不同enantiomer
            mols = generate_chiral_isomers(mol)
            for num,mol in enumerate(mols):
                smi = Chem.MolToSmiles(mol)
                mol = pybel.readstring('smiles', smi)
                if num:
                    mol.title = f'{mol_name}-{num}'
                else:
                    mol.title = f'{mol_name}'
                mol.make3D() 
                outf.write(mol)
    print(f"转换完成，所有分子已保存到 {output_sdf_file}")

if __name__ == "__main__":
    main()
    