import csv

# 假设CIF文件路径
cif_file  = "../../../components.cif"
csv_file = 'output.csv'

# 用于存储提取的数据
data = []

# 读取CIF文件
with open(cif_file, 'r') as f:
    for line in f:
        # 查找包含 "SMILES OpenEye OEToolkits" 的行
        # print(line)
        if 'SMILES_CANONICAL' in line and "OpenEye OEToolkits" in line:
            # 按空格拆分该行
            # print(line,'*'*10)
            parts = line.split()
            # 提取第一个和最后一个数据
            first_data = parts[0]
            last_data = parts[-1]

            # 将数据保存到列表
            data.append([first_data, last_data])
# print(data)
# 将数据写入CSV文件
with open(csv_file, mode='w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['First Data', 'Last Data'])  # 写入CSV表头
    writer.writerows(data)  # 写入数据行

