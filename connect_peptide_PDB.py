import os

def connect_AA_to_peptide(pdb_lines,file_name,output_dir,AA_index_list):
    new_lines_atom,new_lines_conect,idx_atom_type,index_old_new_dic = [],[],{},{}
    new_lines_num = 0
    for line in pdb_lines:
        if line.startswith("ATOM"):
            new_lines_num +=1
            old_line_num = line.split()[1]
            #对原子序号重新编号
            line  = line[:7] + str(new_lines_num).rjust(4) + line[11:]
            #新，旧原子序号映射，便于修正connect信息
            index_old_new_dic[old_line_num] = str(new_lines_num)
            AA_num = line.split()[4]
            AA_type =  f'{AA_num}-{line.split()[2]}'
            #原子序号和元素种类映射，氨基酸之间的酰胺连键
            idx_atom_type.setdefault(AA_type,str(new_lines_num))
            new_lines_atom.append(line)
        else:
            #更新connect 编号
            new_lines = "CONECT"
            for idx in line.split()[1:]:
                new_lines += index_old_new_dic[idx].rjust(5)
            new_lines_conect.append(new_lines+'\n')
    #进行氨基酸之间酰胺连键
    for num in range(len(AA_index_list)-1):
        line = "CONECT" + idx_atom_type[str(AA_index_list[num])+'-C'].rjust(5) + idx_atom_type[str(AA_index_list[num+1])+"-N"].rjust(5) + '\n'
        new_lines_conect.append(line)
    file_path = os.path.join(output_dir, f'{file_name}.pdb')
    with open(file_path, 'w') as outfile:
        outfile.writelines(new_lines_atom)
        outfile.writelines(new_lines_conect)