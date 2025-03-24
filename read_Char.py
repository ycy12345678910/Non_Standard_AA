import pandas as pd
from modify_smi import cano_smiles

def read_defined_symbol():
    df_symbol = pd.read_csv('data/symbol.csv') 
    # df_symbol['cano_smiles'] = df_symbol['SMILES'].apply(cano_smiles)
    symbol_smi_dict = dict(zip(df_symbol['SMILES'],df_symbol['Symbol']))
    return symbol_smi_dict

def read_known_symbol():
    #读取已经用过的三字符，避免覆盖
    with open('data/residuetypes.dat', 'r') as file:
        known_char = [line.split()[0] for line in file if len(line.split()[0]) == 3]
    return known_char
