import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem,rdmolops,Draw
import os,sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from peptide2AA import substructure
from modify_smi import cano_smiles
from modify_PDB import fragment_to_PDB
import os


def read_defined_symbol():
    os.chdir(os.path.dirname(__file__))
    df_symbol = pd.read_csv('symbol.csv') 
    df_symbol['SMILES'] = df_symbol['SMILES'].apply(cano_smiles)
    symbol_smi_dict = dict(zip(df_symbol['SMILES'],df_symbol['Symbol']))
    return symbol_smi_dict

def fragment_N_COO(mol):
    N_smarts = '[N;!H0][H]'
    COO_smarts = 'C(=O)[O;H1]'
    matches_N = substructure(mol,N_smarts)
    matches_COO = substructure(mol,COO_smarts)
    path_min = ['0'] *100
    #找主链酰胺
    for match1 in matches_COO:
        for match2 in matches_N:
            path = AllChem.GetShortestPath(mol,match1[0],match2[0])
            if len(path_min) > len(path):
                path_min = path
    cut_idx = []
    #氨基和羧基不同时存在
    if path_min == ['0'] *100:
        if matches_N:
            cut_idx.append((matches_N[0][0],matches_N[0][1]))
        elif matches_COO:
            cut_idx.append((matches_COO[0][0],matches_COO[0][2]))
        else:
            return []
    for match in matches_COO:
        if path_min[0] in match:
            cut_idx.append((path_min[0],match[-1]))
    for match in matches_N:
        if path_min[-1] in match:
            cut_idx.append((path_min[-1],match[-1]))
            break
    bonds_id = [mol.GetBondBetweenAtoms(x,y).GetIdx() for x,y in cut_idx]
    frags = Chem.FragmentOnBonds(mol,bonds_id,dummyLabels = [(0,0)]*len(bonds_id))
    fragments = Chem.GetMolFrags(frags, asMols=True, sanitizeFrags=True)
    return fragments

def AA_smi_to_pdb(smi,residue_name):    
    if '*N' in smi:
        smi = smi.replace('*N','N').replace('*','O')
    else:
        smi = smi.replace('*C','OC').replace('*','').replace ('()','')
    mol = Chem.MolFromSmiles(smi)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)    
    AllChem.MMFFOptimizeMolecule(mol) 
    fragments = fragment_N_COO(mol)
    # img = Draw.MolsToGridImage(mols)
    # img.show()
    if not os.path.exists('New_AA_PDB'):
        os.makedirs('New_AA_PDB')
    for frag in fragments:
        if frag.GetNumHeavyAtoms() > 1 or residue_name == 'NH2':
            rdmolops.RemoveStereochemistry(frag)
            # frag = Chem.RemoveHs(frag)
            new_pdb_lines,AA_index_list = fragment_to_PDB(frag,residue_name,[1])
            with open(f'New_AA_PDB/{residue_name}.pdb', 'w') as outfile:
                outfile.writelines(new_pdb_lines)
            return Chem.MolToSmiles(frag)

symbol_smi_dict = read_defined_symbol()
frags = {}
for smi,residue_name in  symbol_smi_dict.items():
    frags[residue_name] = AA_smi_to_pdb(smi,residue_name) 
df = pd.DataFrame(list(frags.items()), columns=['Symbol','SMILES'])
df.to_csv('output.csv', index=False)
