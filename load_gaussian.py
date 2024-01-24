import math
import re

# 原子序号和原子名称的映射
atom_map = {17: "Cl", 9: "F ", 35: "Br", 5: "B ", 46: "Pd", 6: "C ", 1: "H ", 7: "N ", 8: "O ", 16: "S "}


# 解析文件，提取坐标
def extract_coordinates(file):
    with open(file, "r") as f:
        rline = f.readlines()

    start = next(
        (i for i in range(len(rline)) if "Standard orientation:" in rline[i]),
        next((i for i in range(len(rline)) if "Input orientation:" in rline[i]), 0)
    )

    end = next((m for m in range(start + 5, len(rline)) if "---" in rline[m]), len(rline))

    return [(atom_map.get(int(words[1]), words[1]), list(map(float, words[3:6])))
            for line in rline[start + 5: end]
            for words in (line.split(),)]

#提取电荷与自旋多重度
def extract_charge_and_multiplicity(filename):
    with open(filename, 'r') as file:
        content = file.read()
        charge_pattern = r"Charge =  (\d+)"
        multiplicity_pattern = r"Multiplicity = (\d+)"
        charge = re.search(charge_pattern, content)
        charge_value = int(charge.group(1)) if charge else None
        multiplicity = re.search(multiplicity_pattern, content)
        multiplicity_value = int(multiplicity.group(1)) if multiplicity else None
        return f"{multiplicity_value}  {charge_value}"


# 计算两个原子之间的距离
def calculate_distance(atom1, atom2):
    coords1 = atom1[1]
    coords2 = atom2[1]
    return math.sqrt((coords1[0] - coords2[0]) ** 2 + (coords1[1] - coords2[1]) ** 2 + (coords1[2] - coords2[2]) ** 2)


# 寻找过渡原子
def transition_atom(atoms):
    pairs = [('C', 1.15, 1.6), ('H', 0.9, 1.6), ('O', 1.0, 1.6), ('N', 1.15, 1.6)]
    for i, atom1 in enumerate(atoms):
        if atom1[0].strip() != 'H':
            continue
        if sum(any(atom2[0].strip() == atom and d_min < calculate_distance(atom1, atom2) < d_max for atom2 in atoms if
                   atom1 != atom2)
               for atom, d_min, d_max in pairs) >= 2:
            return i + 1
    return 0

#提取INDEF信息，需要Z矩阵坐标
def obtain_index(file):
    input_file = open(file, 'r')
    rline = input_file.readlines()
    # 找前十行最后的空格定位起始位置
    for i in range(10):
        if rline[i] in ['\n', '\r\n']:
            start = i + 3
    # 找输入文件结束的地方
    for j in range(start, len(rline) - 10):
        if rline[j] in ['\n', '\r\n']:
            end = j
    a = 1
    b = 2
    c = 3
    bond = []
    angle = []
    diher = []
    for data in range(start, end):
        rdata = rline[data].split()  # 变量替换
        bond.append(str(a + 1) + '-' + rdata[1])
        a = a + 1
    for data in range(start + 1, end):
        rdata = rline[data].split()
        angle.append(str(b + 1) + '-' + rdata[1] + '-' + rdata[3])
        b = b + 1
    for data in range(start + 2, end):
        rdata = rline[data].split()
        diher.append(str(c + 1) + '-' + rdata[1] + '-' + rdata[3] + '-' + rdata[5])
        c = c + 1
    return bond, angle, diher
# 从文件中读取过渡原子信息
def obtain_transition_atom(file):
    coord = extract_coordinates(file)
    index = transition_atom(coord)
    return index