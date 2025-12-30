import numpy as np
import matplotlib.pyplot as plt

def get_neighbors(atoms, lattice_vectors, cutoff):
    """
    找出所有在 cutoff 距离内的邻居，并返回 (atom_i, atom_j, R_vector, distance)
    """
    n_atoms = len(atoms)
    neighbors = []
    # 搜索邻居单元格 (简单立方扩展)
    for dx in [-1, 0, 1]:
        for dy in [-1, 0, 1]:
            for dz in [-1, 0, 1]:
                R = dx * lattice_vectors[0] + dy * lattice_vectors[1] + dz * lattice_vectors[2]
                for i in range(n_atoms):
                    for j in range(n_atoms):
                        # 避免同一原子在 R=0 时自比较
                        if dx == 0 and dy == 0 and dz == 0 and i == j:
                            continue
                        
                        pos_i = atoms[i]
                        pos_j = atoms[j] + R
                        dist = np.linalg.norm(pos_i - pos_j)
                        
                        if dist < cutoff:
                            neighbors.append((i, j, np.array([dx, dy, dz]), dist))
    return neighbors

def construct_hamiltonian(k, n_atoms, neighbors, t1, e0=0):
    """
    构造紧束缚哈密顿矩阵 H(k)
    """
    H = np.zeros((n_atoms, n_atoms), dtype=complex)
    np.fill_diagonal(H, e0)
    
    for i, j, R_idx, dist in neighbors:
        # 相位因子 e^(i k . R_real)
        # 注意：这里的 R_idx 是单元格索引 [nx, ny, nz]
        # 实际位移是 R_real = nx*a1 + ny*a2 + nz*a3
        # 我们假设 lattice 是简单的立方 a * I
        phase = np.exp(1j * 2 * np.pi * np.dot(k, R_idx))
        H[i, j] += t1 * phase
        
    return H

def main():
    # 晶格常数
    a = 3.501
    lattice_vectors = np.eye(3) * a
    
    # 12d 轨道在 conventional cell 中的分数坐标 (基于 VASP 文件识别出的 12 个 H 原子)
    # 我们只取 lines 3-14
    h_atoms_frac = np.array([
        [0.25, 0.00, 0.50],
        [0.25, 0.50, 0.00],
        [0.50, 0.25, 0.00],
        [0.00, 0.25, 0.50],
        [0.50, 0.00, 0.75],
        [0.00, 0.50, 0.75],
        [0.75, 0.50, 0.00],
        [0.75, 0.00, 0.50],
        [0.00, 0.75, 0.50],
        [0.50, 0.75, 0.00],
        [0.00, 0.50, 0.25],
        [0.50, 0.00, 0.25]
    ])
    
    h_atoms_cart = h_atoms_frac * a
    n_atoms = len(h_atoms_cart)
    
    # 寻找最近邻
    # d_nn = a * sqrt(2)/4 = 1.2377
    cutoff = 1.3
    neighbors = get_neighbors(h_atoms_cart, lattice_vectors, cutoff)
    
    print(f"找到 {len(neighbors)} 条跳跃路径 (含跨胞)。")
    # 检查配位数
    coordination = [0] * n_atoms
    for i, j, R, d in neighbors:
        coordination[i] += 1
    print(f"每个原子的配位数: {coordination}")

    # 紧束缚参数 (根据 H-H 距离估算，或取典型值 -2.0 eV)
    t1 = -2.0 
    e0 = 0.0
    
    # 定义高对称路径 G-X-M-G-R
    # 坐标以 2pi/a 为单位
    points = {
        'G': [0, 0, 0],
        'X': [0.5, 0, 0],
        'M': [0.5, 0.5, 0],
        'R': [0.5, 0.5, 0.5]
    }
    path_segments = [('G', 'X'), ('X', 'M'), ('M', 'G'), ('G', 'R')]
    n_k = 100
    
    k_list = []
    x_ticks = [0]
    current_x = 0
    
    for start_node, end_node in path_segments:
        start = np.array(points[start_node])
        end = np.array(points[end_node])
        for i in range(n_k):
            frac = i / n_k
            k = start + (end - start) * frac
            k_list.append(k)
        current_x += np.linalg.norm(end - start)
        x_ticks.append(current_x)
    k_list.append(np.array(points[path_segments[-1][1]]))
    
    all_energies = []
    for k in k_list:
        H = construct_hamiltonian(k, n_atoms, neighbors, t1, e0)
        eigenvalues = np.linalg.eigvalsh(H)
        all_energies.append(eigenvalues)
        
    all_energies = np.array(all_energies)
    
    # 绘图
    plt.figure(figsize=(8, 6))
    for i in range(n_atoms):
        plt.plot(all_energies[:, i], color='b', lw=1.5)
    
    # 设置刻度
    k_dist = [0]
    curr = 0
    for i in range(len(k_list)-1):
        curr += np.linalg.norm(k_list[i+1] - k_list[i])
        k_dist.append(curr)
    
    plt.clf() # 重新画以使用正确的 x 轴距离
    for i in range(n_atoms):
        plt.plot(k_dist, all_energies[:, i], color='tab:blue', alpha=0.8)
    
    tick_pos = [0]
    curr = 0
    for start_node, end_node in path_segments:
        curr += np.linalg.norm(np.array(points[end_node]) - np.array(points[start_node]))
        tick_pos.append(curr)
        
    plt.xticks(tick_pos, ['$\Gamma$', 'X', 'M', '$\Gamma$', 'R'])
    for pos in tick_pos:
        plt.axvline(pos, color='k', ls='--', alpha=0.3)
    
    plt.ylabel('Energy (eV)')
    plt.title('Tight-Binding Band Structure of CaH6 Hydrogen Cage (12d sites)')
    plt.grid(True, alpha=0.2)
    # 能带图已保存
    plt.savefig('CaH6_H_cage_TB_bands.png')
    print("能带图已保存至 CaH6_H_cage_TB_bands.png")

    # Gamma 点成键分析
    H_gamma = construct_hamiltonian(np.array([0,0,0]), n_atoms, neighbors, t1, e0)
    evals, evecs = np.linalg.eigh(H_gamma)
    print("\nGamma 点能级与成键性质 (t1 = -2.0 eV):")
    for i in range(n_atoms):
        char = "成键 (Bonding)" if evals[i] < e0 else "反键 (Anti-bonding)"
        print(f"  能级 {i+1}: {evals[i]:.4f} eV - {char}")

if __name__ == "__main__":
    main()
