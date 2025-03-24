import itertools

def get_current_position(position_file):
    try:
        with open(position_file, "r") as f:
            return int(f.read().strip())  
    except FileNotFoundError:
        return 000 

def save_position(position_file,position):
    with open(position_file, "w") as f:
        f.write(str(position))  

def Generate_Char():
    characters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789'
    # 使用itertools.product生成所有3个字符的排列组合
    combinations = itertools.product(characters, repeat=3)
    position_file = "data/position.txt"
    current_position = get_current_position(position_file)
    Three_Char = next(itertools.islice(combinations, current_position, current_position+1))
    save_position(position_file,current_position+1)
    # 将元组转化为字符串并打印
    return  ''.join(Three_Char)

def Name_new_AA(Known_Char):
    while True:
        residue_name = Generate_Char()
        if residue_name not in Known_Char:
            return residue_name

