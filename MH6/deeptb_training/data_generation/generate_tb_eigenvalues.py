"""
紧束缚模型能带计算核心函数
用于生成DeePTB训练数据
"""
import numpy as np
from typing import List, Tuple, Dict

def get_neighbors(atoms: np.ndarray, lattice_vectors: np.ndarray, cutoff: float) -> List[Tuple]:
    """
    找出所有在cutoff距离内的邻居

    Parameters:
    -----------
    atoms : np.ndarray, shape (n_atoms, 3)
        原子笛卡尔坐标
    lattice_vectors : np.ndarray, shape (3, 3)
        晶格向量
    cutoff : float
        截断距离(Å)

    Returns:
    --------
    neighbors : List[Tuple[int, int, np.ndarray, float]]
        (atom_i, atom_j, R_vector, distance)
    """
    n_atoms = len(atoms)
    neighbors = []

    # 搜索邻居单元格
    for dx in [-1, 0, 1]:
        for dy in [-1, 0, 1]:
            for dz in [-1, 0, 1]:
                R = dx * lattice_vectors[0] + dy * lattice_vectors[1] + dz * lattice_vectors[2]
                for i in range(n_atoms):
                    for j in range(n_atoms):
                        # 避免同一原子在R=0时自比较
                        if dx == 0 and dy == 0 and dz == 0 and i == j:
                            continue

                        pos_i = atoms[i]
                        pos_j = atoms[j] + R
                        dist = np.linalg.norm(pos_i - pos_j)

                        if dist < cutoff:
                            neighbors.append((i, j, np.array([dx, dy, dz]), dist))
    return neighbors


def construct_hamiltonian(k: np.ndarray, n_atoms: int, neighbors: List[Tuple],
                         hopping_params: Dict, onsite_energies: np.ndarray) -> np.ndarray:
    """
    构造紧束缚哈密顿矩阵H(k)

    Parameters:
    -----------
    k : np.ndarray, shape (3,)
        k点坐标(分数坐标)
    n_atoms : int
        原子数量
    neighbors : List[Tuple]
        邻居列表
    hopping_params : Dict
        跳跃积分参数，例如 {('H', 'H'): -2.0}
    onsite_energies : np.ndarray, shape (n_atoms,)
        onsite能级

    Returns:
    --------
    H : np.ndarray, shape (n_atoms, n_atoms), dtype=complex
        k空间哈密顿矩阵
    """
    H = np.zeros((n_atoms, n_atoms), dtype=complex)
    np.fill_diagonal(H, onsite_energies)

    for i, j, R_idx, dist in neighbors:
        # 相位因子 e^(i k . R)
        phase = np.exp(1j * 2 * np.pi * np.dot(k, R_idx))

        # 根据原子类型获取跳跃参数（这里简化为统一参数）
        # 如果需要区分Ca-H, H-H等，可以扩展
        for bond_type, t_value in hopping_params.items():
            H[i, j] += t_value * phase

    return H


def generate_kpath(points_dict: Dict[str, List[float]],
                  path_segments: List[Tuple[str, str]],
                  n_k: int = 100) -> Tuple[np.ndarray, List[float], List[str]]:
    """
    生成高对称路径上的k点

    Parameters:
    -----------
    points_dict : Dict
        高对称点字典，例如 {'G': [0,0,0], 'X': [0.5,0,0]}
    path_segments : List[Tuple[str, str]]
        路径段，例如 [('G', 'X'), ('X', 'M')]
    n_k : int
        每段的k点数量

    Returns:
    --------
    k_list : np.ndarray, shape (total_k_points, 3)
        k点列表
    k_distances : List[float]
        k点距离（用于绘图x轴）
    tick_positions : List[float]
        高对称点位置
    """
    k_list = []
    k_distances = [0]
    tick_positions = [0]
    current_dist = 0

    for start_node, end_node in path_segments:
        start = np.array(points_dict[start_node])
        end = np.array(points_dict[end_node])

        # 生成该段的k点
        for i in range(n_k):
            frac = i / n_k
            k = start + (end - start) * frac
            k_list.append(k)

            if len(k_list) > 1:
                dk = np.linalg.norm(k_list[-1] - k_list[-2])
                current_dist += dk
                k_distances.append(current_dist)

        # 记录段终点位置
        segment_length = np.linalg.norm(end - start)
        tick_positions.append(tick_positions[-1] + segment_length)

    # 添加最后一个点
    k_list.append(np.array(points_dict[path_segments[-1][1]]))
    if len(k_list) > 1:
        dk = np.linalg.norm(k_list[-1] - k_list[-2])
        current_dist += dk
        k_distances.append(current_dist)

    return np.array(k_list), k_distances, tick_positions


def calculate_bands(atoms: np.ndarray, lattice_vectors: np.ndarray,
                   k_list: np.ndarray, hopping_params: Dict,
                   onsite_energies: np.ndarray, cutoff: float) -> np.ndarray:
    """
    计算能带结构

    Parameters:
    -----------
    atoms : np.ndarray, shape (n_atoms, 3)
        原子位置(笛卡尔坐标)
    lattice_vectors : np.ndarray, shape (3, 3)
        晶格向量
    k_list : np.ndarray, shape (n_k, 3)
        k点列表
    hopping_params : Dict
        跳跃参数
    onsite_energies : np.ndarray, shape (n_atoms,)
        onsite能级
    cutoff : float
        跳跃截断距离

    Returns:
    --------
    eigenvalues : np.ndarray, shape (n_k, n_atoms)
        所有k点的本征值
    """
    n_atoms = len(atoms)
    neighbors = get_neighbors(atoms, lattice_vectors, cutoff)

    all_eigenvalues = []
    for k in k_list:
        H = construct_hamiltonian(k, n_atoms, neighbors, hopping_params, onsite_energies)
        eigenvalues = np.linalg.eigvalsh(H)  # 厄米矩阵的本征值
        all_eigenvalues.append(eigenvalues)

    return np.array(all_eigenvalues)


def generate_deeptb_data(structure_file: str, element_list: List[str],
                        hopping_params: Dict, onsite_energies: Dict,
                        cutoff: float, n_k_per_segment: int = 100) -> Tuple[np.ndarray, np.ndarray]:
    """
    从结构文件生成DeePTB训练数据

    Parameters:
    -----------
    structure_file : str
        VASP POSCAR格式结构文件
    element_list : List[str]
        元素列表，例如 ['H'] 或 ['Ca', 'H']
    hopping_params : Dict
        跳跃参数，例如 {('H', 'H'): -2.0}
    onsite_energies : Dict
        onsite能级，例如 {'H': 0.0, 'Ca': -2.0}
    cutoff : float
        跳跃截断距离
    n_k_per_segment : int
        每段k路径的点数

    Returns:
    --------
    eigenvalues : np.ndarray, shape (1, total_k_points, n_bands)
        本征值数组（DeePTB格式：[nframes, nkpoints, nbands]）
    kpoints : np.ndarray, shape (total_k_points, 3)
        k点坐标
    """
    from ase.io import read

    # 读取结构
    structure = read(structure_file)
    lattice_vectors = structure.cell.array
    positions = structure.positions

    # 构造onsite能级数组
    onsite_array = np.array([onsite_energies[elem] for elem in structure.get_chemical_symbols()])

    # 定义高对称路径（立方系统）
    points = {
        'G': [0, 0, 0],
        'X': [0.5, 0, 0],
        'M': [0.5, 0.5, 0],
        'R': [0.5, 0.5, 0.5]
    }
    path_segments = [('G', 'X'), ('X', 'M'), ('M', 'G'), ('G', 'R')]

    # 生成k点路径
    k_list, k_distances, tick_positions = generate_kpath(points, path_segments, n_k_per_segment)

    # 计算能带
    eigenvalues = calculate_bands(positions, lattice_vectors, k_list,
                                  hopping_params, onsite_array, cutoff)

    # 格式化为DeePTB要求的shape: [nframes, nkpoints, nbands]
    eigenvalues_formatted = eigenvalues[np.newaxis, :, :]  # shape: (1, n_k, n_bands)

    return eigenvalues_formatted, k_list


if __name__ == "__main__":
    # 测试代码
    print("紧束缚能带计算模块加载成功")
    print("主要函数：")
    print("  - generate_deeptb_data(): 从结构文件生成训练数据")
    print("  - calculate_bands(): 计算能带结构")
    print("  - generate_kpath(): 生成k点路径")
