from rdkit import Chem
import re

def substructure(mol,smarts_amido):
    pattern = Chem.MolFromSmarts(smarts_amido)
    matches = mol.GetSubstructMatches(pattern)
    return matches

def find_N_terminal(fragments_smi_lable):
    for frag_idx,smi in enumerate(fragments_smi_lable):
        if smi.count('*') == 1 and '*]C' in smi :
            #定位酰胺键N原子的原子序号
            connect_N_idx = re.search(r'\[([^\]]+)\*\]', smi).group(1)
            return frag_idx,connect_N_idx

def find_amido_C(connect_N_idx,dummy_idxs,connect_C_idx,sort_key,count):
    for frag_idx,dummy_idx in enumerate(dummy_idxs):
        if connect_C_idx in dummy_idx:
            sort_key[frag_idx] = count
            if count < len(dummy_idxs)-1:
                connect_N_idx  = [idx for idx in dummy_idx if idx != connect_C_idx][0]
    return connect_N_idx,sort_key
def sorted_AA(fragments_smi_lable,fragments,COO_N_idx_list):
    #找N端氨基酸
    sort_key = ['']*len(fragments_smi_lable)
    frag_idx,connect_N_idx = find_N_terminal(fragments_smi_lable)
    sort_key[frag_idx] = 0
    COO_N_idx_list = [(str(x),str(y)) for  x,y in COO_N_idx_list]
    count = 0
    dummy_idxs = [re.findall(r'\[([^\]]+)\*\]', smi) for smi in fragments_smi_lable ]
    #从N端氨基酸开始，对氨基酸进行编号
    while count < len(COO_N_idx_list):
        count += 1
        #定位氨基酸羧基C原子的原子序号
        connect_idx = [idx for idx in COO_N_idx_list if connect_N_idx in idx][0]
        connect_C_idx = [idx for idx in connect_idx if idx != COO_N_idx_list][0]
        connect_N_idx,sort_key = find_amido_C(COO_N_idx_list,dummy_idxs,connect_C_idx,sort_key,count)
            
    #从N端氨基酸开始，对氨基酸进行排序
    # print(f"氨基酸顺序:{sort_key}")
    sorted_fragments = [x for  x,_ in sorted(zip(fragments, sort_key), key=lambda pair: pair[1])]
    return sorted_fragments
def cut_peptide(COO_N_idx,mol):
    bonds_id = [mol.GetBondBetweenAtoms(x,y).GetIdx() for x,y in COO_N_idx]
    frags = Chem.FragmentOnBonds(mol,bonds_id,True,dummyLabels = [(0,0)]*len(bonds_id)) # 切割得到碎片
    fragments = Chem.GetMolFrags(frags, asMols=True, sanitizeFrags=True)
    #判断N端位置，并进行调整为N端开始
    frags_lable = Chem.FragmentOnBonds(mol,bonds_id,True) # 切割得到碎片
    fragments_lable = Chem.GetMolFrags(frags_lable, asMols=True, sanitizeFrags=True)
    fragments_smi_lable = [Chem.MolToSmiles(frag) for frag in fragments_lable ]
    # 对编号为0的dummy原子设置标签为[0*]
    fragments_smi_lable = [re.sub(r'\*(?!])', '[0*]', smi) for smi in fragments_smi_lable ]
    return fragments,fragments_smi_lable
def peptide_to_AA(mol):
    smarts_amido = ['[#6]C(!@-N[#6])=O','O=C(!@-[N;H2])[#6]N']
    COO_N_idx_list = []
    for amido in smarts_amido:
        matches = substructure(mol,amido)
        COO_N_idx_list += [(x,y) for  _,x,y,_,_ in matches]
    COO_N_idx_list = set(COO_N_idx_list)
    fragments,fragments_smi_lable = cut_peptide(COO_N_idx_list,mol)
    dummy_idxs = [re.findall(r'\[([^\]]+)\*\]', smi) for smi in fragments_smi_lable ]
    #检查是否有侧链剪切
    has_side_chain_list = [dummy_idx for dummy_idx in dummy_idxs if len(dummy_idx) >2 ]    
    for dummy_idxs in has_side_chain_list:
        for dummy_idx in dummy_idxs:
            new_COO_N_idx_list = [idxs for idxs in COO_N_idx_list if int(dummy_idx) not in idxs]
            fragments,fragments_smi_lable = cut_peptide(new_COO_N_idx_list,mol)
            dummy_idxs = [re.findall(r'\[([^\]]+)\*\]', smi) for smi in fragments_smi_lable ]
            dummy_idxs_terminal = [dummy_idx[0] for dummy_idx in dummy_idxs if len(dummy_idx) ==1 ]
            if len(dummy_idxs_terminal) == 2:
                COO_N_idx_list = new_COO_N_idx_list
                break
    #对多肽整理为N端到C端
    sorted_fragments = sorted_AA(fragments_smi_lable,fragments,COO_N_idx_list)
    return sorted_fragments

