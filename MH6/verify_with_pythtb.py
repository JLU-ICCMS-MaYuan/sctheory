import sys
import os
import numpy as np
import matplotlib.pyplot as plt

# 将 pythtb 路径加入 sys.path 并放在最前面
# pythtb 位于项目根目录，而本脚本位于 MH6
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'pythtb')))
from pythtb import Lattice, TBModel

def main():
    # 1. 定义晶格
    a = 3.501
    lat = [[a, 0, 0], [0, a, 0], [0, 0, a]]
    
    # 2. 12d 轨道的分数坐标
    orb = [
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
    ]
    
    # 3. 创建 PythTB 模型
    lat_obj = Lattice(lat, orb, periodic_dirs="all")
    my_model = TBModel(lat_obj)
    
    # 4. 设置跳跃 (Hopping)
    # 我们使用 generate_tb_model.py 里的逻辑找到邻居
    t1 = -2.0
    cutoff = 1.3
    
    # 转换为笛卡尔坐标以便计算距离
    h_atoms_cart = np.array(orb) * a
    n_atoms = len(orb)
    lattice_vectors = np.array(lat)
    
    count_hop = 0
    for dx in [-1, 0, 1]:
        for dy in [-1, 0, 1]:
            for dz in [-1, 0, 1]:
                R = dx * lattice_vectors[0] + dy * lattice_vectors[1] + dz * lattice_vectors[2]
                for i in range(n_atoms):
                    for j in range(n_atoms):
                        # 为了避免 pythtb 重复计算共轭对，我们只添加特定顺序的跳跃
                        # 例如只添加 i < j 或者在 i==j 时限制 R 的方向
                        if dx == 0 and dy == 0 and dz == 0 and i == j:
                            continue
                        
                        # 唯一性过滤：只添加 (i,j,R) 使得其共轭对 (j,i,-R) 不会被重复添加
                        if i < j or (i == j and (dx > 0 or (dx == 0 and dy > 0) or (dx == 0 and dy == 0 and dz > 0))):
                            pos_i = h_atoms_cart[i]
                            pos_j = h_atoms_cart[j] + R
                            dist = np.linalg.norm(pos_i - pos_j)
                            
                            if dist < cutoff:
                                # pythtb 的 set_hop 参数: 能量, 轨道i, 轨道j, 晶格位移R
                                my_model.set_hop(t1, i, j, [dx, dy, dz])
                                count_hop += 1
    
    print(f"PythTB 设置了 {count_hop} 条跳跃路径ảng。\n")
    
    # 5. 定义高对称路径
    path = [
        [0, 0, 0],     # G
        [0.5, 0, 0],   # X
        [0.5, 0.5, 0], # M
        [0, 0, 0],     # G
        [0.5, 0.5, 0.5]# R
    ]
    labels = [r'$\Gamma$', 'X', 'M', r'$\Gamma$', 'R']
    
    # 插值点
    (k_vec, k_dist, k_node) = my_model.k_path(path, 400)
    
    # 6. 求解能带
    evals = my_model.solve_ham(k_vec)
    print(f"evals shape: {evals.shape}")
    
    # 7. 绘图
    fig, ax = plt.subplots(figsize=(8, 6))
    for i in range(n_atoms):
        ax.plot(k_dist, evals[:, i], color='r', lw=1.5, alpha=0.7, label='PythTB' if i == 0 else "")
    
    ax.set_xticks(k_node)
    ax.set_xticklabels(labels)
    for node in k_node:
        ax.axvline(node, color='k', ls='--', alpha=0.3)
    
    ax.set_ylabel('Energy (eV)')
    ax.set_title('CaH6 Hydrogen Cage (PythTB vs Manual)')
    ax.grid(True, alpha=0.2)
    
    # 8. 运行原始模型进行对比
    # 我们直接在这里实现简单的能带计算，或者读取 generate_tb_model.py 的逻辑
    # 为了验证，我们在这里手动计算几个点
    print("\n--- 数值对比 (Gamma点) ---")
    k_gamma = np.array([0, 0, 0])
    H_gamma_pythtb = my_model.hamiltonian(k_gamma)
    evals_gamma_pythtb = np.linalg.eigvalsh(H_gamma_pythtb)
    print("PythTB Gamma evals:", evals_gamma_pythtb)
    
    # 手动计算结果 (来自 generate_tb_model.py)
    evals_gamma_manual = np.array([-8.0, -4.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 4.0, 8.0])
    print("Manual Gamma evals:", evals_gamma_manual)
    
    diff = np.max(np.abs(evals_gamma_pythtb - evals_gamma_manual))
    print(f"\nGamma点最大差异: {diff:.2e}")
    if diff < 1e-10:
        print("验证通过：手动模型与 PythTB 结果一致！")
    else:
        print("验证失败：存在数值差异。")

    # 绘图保存
    plt.savefig('MH6/compare_pythtb_manual.png')
    print("\n对比图已保存至 MH6/compare_pythtb_manual.png")

if __name__ == "__main__":
    main()
