import sys

def read_xyz_and_format(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
    
    # 跳过第一行（原子数）和第二行（注释行）
    atoms = []
    for line in lines[2:]:
        parts = line.strip().split()
        if len(parts) >= 4:
            symbol = parts[0]
            x, y, z = map(float, parts[1:4])
            atoms.append((symbol, x, y, z))
    
    # 格式化输出
    for i, (symbol, x, y, z) in enumerate(atoms, start=1):
        print(f"A{i} = {{{x:.6f}, {y:.6f}, {z:.6f}}};")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <xyz_file>")
        sys.exit(1)
    
    xyz_file = sys.argv[1]
    read_xyz_and_format(xyz_file)
