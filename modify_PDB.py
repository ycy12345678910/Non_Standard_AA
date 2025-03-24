from rdkit.Chem import AllChem,rdmolfiles
from peptide2AA import substructure

def find_dummy(pdb_lines):
    #避免端基氨基酸只有一个剪切位点
    dummy_idx_2,connect_idx2 = -1,-1
    new_pdb_lines= []
    for line in pdb_lines:
        if '*1' in line:
            dummy_idx_1 = line.strip().split()[1]
            continue
        elif '*2' in line:
            dummy_idx_2 = line.strip().split()[1]
            continue
        elif 'CONECT' in line and str(dummy_idx_1) in line:
            #第一个是剪切原子，最后一个是dummy原子
            connect_idx1 = int(line.strip().split()[1]) -1
            line = line[:-5]
        elif 'CONECT' in line and str(dummy_idx_2) in line:
            connect_idx2 = int(line.strip().split()[1]) -1
            line = line[:-5]
        new_pdb_lines.append(line)
    return connect_idx1,connect_idx2,new_pdb_lines

def find_another_connect(match_smarts,pdb_lines,atom_type_dic,mol,connect_idx1):
    matches = substructure(mol,match_smarts)
    if len(matches) == 1:
        connect_idx2 = matches[0][0]
        atom_type_dic.setdefault(connect_idx2,pdb_lines[connect_idx2][-3])
    else:
        #找最短路径
        path_min = ['0'] *100
        for match in matches:
            path = AllChem.GetShortestPath(mol,match[0],connect_idx1)
            if len(path) < len(path_min):
                path_min = path
        if len(path_min) < 100:
            connect_idx2 = path_min[0]
            atom_type_dic.setdefault(connect_idx2,pdb_lines[connect_idx2][-3])
    return atom_type_dic

def find_N_COO_idx(pdb_lines,mol):
    #对氨基N和羧基C定义原子类型
    connect_idx1,connect_idx2,new_pdb_lines = find_dummy(pdb_lines)
    atom_type_dic = {}
    #对氨基N和羧基C从PDB中的元素符号获得
    atom_type_dic.setdefault(connect_idx1,pdb_lines[connect_idx1][-3])
    if connect_idx2 > -1 :
        atom_type_dic.setdefault(connect_idx2,pdb_lines[connect_idx2][-3])
    #如果只有一个剪切位点，需要找到另一个氨基或者羧基
    elif atom_type_dic[connect_idx1] == 'C':
        match_smarts ='[N;!H0]'
        atom_type_dic = find_another_connect(match_smarts,pdb_lines,atom_type_dic,mol,connect_idx1)
    else:
        match_smarts = 'C(=O)[O;!H0,-1]'
        atom_type_dic = find_another_connect(match_smarts,pdb_lines,atom_type_dic,mol,connect_idx1)
    return atom_type_dic,new_pdb_lines

def get_first_level_neighbors(current_atom,visited,atom_type_dic,depth):
    neighbor_idx = []
    for neighbor in current_atom.GetNeighbors():
        #定义羧基氧
        if depth ==1 and neighbor.GetSymbol() =="O":
            symbol = [neighbor.GetSymbol() for neighbor in current_atom.GetNeighbors()]
            visited.append(neighbor.GetIdx())
            #判断是否是C端羧基
            if symbol.count("O") ==1:
                atom_type_dic.setdefault(neighbor.GetIdx(),'O')
            elif neighbor.GetFormalCharge():
                atom_type_dic.setdefault(neighbor.GetIdx(),'OC2' )
            else:
                atom_type_dic.setdefault(neighbor.GetIdx(),'OC1' )
        elif neighbor.GetAtomicNum() > 1 and neighbor.GetIdx() not in visited:
            neighbor_idx.append(neighbor.GetIdx())
            visited.append(neighbor.GetIdx())
    return visited,neighbor_idx,atom_type_dic

def sorted_neighor(neighbor_idx,mol):
    neighbor_Atomic_num = {}
    for idx in neighbor_idx:
        neighbor_Atomic_num[idx] = mol.GetAtomWithIdx(idx).GetAtomicNum() 
    #如果当层序号相同，比较下一层
    if len(set(neighbor_Atomic_num.values())) ==1:
        neighbor_new = {}
        for idx in neighbor_idx:
            Atomic_num = sorted([atom.GetAtomicNum() for atom in mol.GetAtomWithIdx(idx).GetNeighbors() ],reverse=True)
            neighbor_new[idx] = Atomic_num
        neighbor_new = dict(sorted(neighbor_new.items(), key=lambda x: x[1], reverse=True))
        neighbor_idx = list(neighbor_new.keys())
    else:
        sorted_neighbor_Atomic_num = dict(sorted(neighbor_Atomic_num.items(),key=lambda x: x[1], reverse=True))
        neighbor_idx = list(sorted_neighbor_Atomic_num.keys())
    return neighbor_idx

def get_neighbors_without_h(atom_type_dic, mol):
    #如果有羧基C，则从羧基C开始命名，否则从氨基N开始命名
    if 'C' in atom_type_dic.values():
        start_idx  = [key for key, value in atom_type_dic.items() if value == 'C']
    else:
        start_idx =[key for key, value in atom_type_dic.items() if value == 'N']
    atom = mol.GetAtomWithIdx(start_idx[0])
    neighbors = {0: [atom.GetIdx()]} 
    #定义已遍历原子编号
    visited = list(atom_type_dic.keys())
    num_heavy_atoms = mol.GetNumHeavyAtoms()
    depth = 0
    # 循环遍历多度邻居
    while len(visited) < num_heavy_atoms:
        depth += 1
        neighbors[depth] = []
        for idx in neighbors[depth - 1]:
            current_atom = mol.GetAtomWithIdx(idx)
            visited,neighbor_idx,atom_type_dic = get_first_level_neighbors(current_atom,visited,atom_type_dic,depth)
            #不止一个邻居时，需要排序
            if len(neighbor_idx) > 1:   
                neighbor_idx = sorted_neighor(neighbor_idx,mol) 
            neighbors[depth] += neighbor_idx
        #酰胺上的修饰编号
        if neighbors[depth] == []:
            N_idx =[key for key, value in atom_type_dic.items() if value == 'N']
            neighbors[depth].append(N_idx[0])
            # print(f'酰胺上的修饰编号{neighbors}')
    return neighbors,atom_type_dic

def adjust_standard_AA(AA_name,atom_type_dic):
    if AA_name == "TRP":
        for key,value in atom_type_dic.items():
            if value == 'CZ1':
                atom_type_dic[key] = 'CZ2'
            elif value =='CZ2':
                atom_type_dic[key] = 'CZ3'
            elif value =='CH':
                atom_type_dic[key] = 'CH2'
    elif AA_name == "ILE":
        for key,value in atom_type_dic.items():
            if value == 'CD':
                atom_type_dic[key] = 'CD1'
    return atom_type_dic
    
def define_heavy_atom_name(neighbors,mol,atom_type_dic,AA_name):
    alphabet = ['A', 'B', 'G', 'D', 'E', 'Z', 'H', 'Q', 'I', 'K', 'L', 'M', 'N', 'X', 'O', 'P', 'R', 'S', 'T', 'U', 'F', 'C', 'Y', 'W']
    for depth, neighbor in neighbors.items():
        for idx,atom_idx in enumerate(neighbor):
            if len(neighbor) >1:
                atom_type_dic.setdefault(atom_idx,mol.GetAtomWithIdx(atom_idx).GetSymbol() +alphabet[depth-1] +str(idx+1) )
            else:
                atom_type_dic.setdefault(atom_idx,mol.GetAtomWithIdx(atom_idx).GetSymbol() +alphabet[depth-1] )
    atom_type_dic = adjust_standard_AA(AA_name,atom_type_dic)
    return atom_type_dic

def define_H_atom_name(atom_type_dic,mol):
    heavy_atom_idx = list(atom_type_dic.keys())
    for atom_idx in heavy_atom_idx:
        neighbors_H = []
        atom = mol.GetAtomWithIdx(atom_idx)
        neighbors_H = [neighbor_atom.GetIdx() for neighbor_atom in atom.GetNeighbors() if neighbor_atom.GetAtomicNum() == 1 ] 
        for idx,neighbor_idx in enumerate(neighbors_H):
            #定义氨基上的H
            if  atom_type_dic[atom_idx] =='N':
                if len(neighbors_H) == 1:
                    atom_type_dic.setdefault(neighbor_idx,'H')
                else:
                    atom_type_dic.setdefault(neighbor_idx,'H' + str(idx+1))
            elif atom_type_dic[atom_idx] =='O':
                atom_type_dic.setdefault(neighbor_idx,'HO')
            elif len(neighbors_H) == 1:
                atom_type_dic.setdefault(neighbor_idx,'H'+ atom_type_dic[atom_idx][1:])
            else:
                atom_type_dic.setdefault(neighbor_idx, 'H'+ atom_type_dic[atom_idx][1:]+str(idx+1) )
    return atom_type_dic

def define_atom_name(neighbors,mol,atom_type_dic,AA_name):
        atom_type_dic = define_heavy_atom_name(neighbors,mol,atom_type_dic,AA_name)
        atom_type_dic = define_H_atom_name(atom_type_dic,mol)
        return atom_type_dic

def modify_atom_type(atom_type_dic,lines,AA_name,AA_index):
    new_pdb_lines= []
    atom_num = 0
    for idx,line in enumerate(lines):
        if line.startswith('ATOM') or line.startswith('HETATM'):
            idx = idx -1
            line  = line[:12] + atom_type_dic[atom_num].center(4) + line[16:20] + str(AA_index).rjust(6) + line[26:]
            line = line.replace('UNL',AA_name)
            atom_num += 1
        line = line +' \n'
        
        new_pdb_lines.append(line)
    return new_pdb_lines

def fragment_to_PDB(frag,AA_name,AA_index_list):
    if frag.GetNumHeavyAtoms():
        pdb_block = rdmolfiles.MolToPDBBlock(frag)
        pdb_lines = pdb_block.splitlines()
        pdb_lines = [line.replace('HETATM','ATOM  ') for line in pdb_lines if line.startswith('CONECT') or line.startswith('HETATM') ]
        atom_type_dic,new_pdb_lines = find_N_COO_idx(pdb_lines,frag)
        neighbors,atom_type_dic = get_neighbors_without_h(atom_type_dic, frag)
        #删除N的neighbors
        atom_type_dic = define_atom_name(neighbors,frag,atom_type_dic,AA_name)
        #判断是否是封端氨基酸，若是标记为*R
        AA_num = AA_index_list[-1]
        # print(atom_type_dic)
        if set(['C','N']).issubset(set(atom_type_dic.values() )):
            new_pdb_lines = modify_atom_type(atom_type_dic,new_pdb_lines,AA_name,AA_num)
            AA_index_list.append(AA_num+1)
        else:
            if AA_num ==1:
                AA_num = f'{AA_num}R' 
            else:
                 AA_num = f'{AA_num-1}R' 
            new_pdb_lines = modify_atom_type(atom_type_dic,new_pdb_lines,AA_name,AA_num)
            AA_index_list.insert(-1, AA_num)
        return new_pdb_lines,AA_index_list