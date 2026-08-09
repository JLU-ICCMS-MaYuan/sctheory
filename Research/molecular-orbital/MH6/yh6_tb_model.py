import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh

class YH6TightBinding:
    def __init__(self, a0, y_pos, h_pos, params):
        """
        a0: lattice constant
        y_pos: Nx3 array of Y fractional positions
        h_pos: Mx3 array of H fractional positions
        params: dict containing onsite and SK hopping parameters
        """
        self.a0 = a0
        self.lattice = np.eye(3) * a0
        self.y_pos = np.array(y_pos)
        self.h_pos = np.array(h_pos)
        self.all_pos = np.vstack([self.y_pos, self.h_pos])
        self.params = params
        
        # Orbitals: Y -> [s, dz2, dx2-y2, dxy, dyz, dxz] (6), H -> [s] (1)
        self.n_y = len(y_pos)
        self.n_h = len(h_pos)
        self.dim = self.n_y * 6 + self.n_h * 1
        
        # Build orbital map for indexing
        # 0..n_y*6-1 are Y orbitals, n_y*6..dim-1 are H orbitals
        self.H_onsite = np.zeros((self.dim, self.dim))
        for i in range(self.n_y):
            start = i * 6
            self.H_onsite[start, start] = params['e_y_s']
            for j in range(1, 6):
                self.H_onsite[start + j, start + j] = params['e_y_d']
        for i in range(self.n_h):
            idx = self.n_y * 6 + i
            self.H_onsite[idx, idx] = params['e_h_s']

    def get_sk_block(self, l, m, n, type_pair, p):
        if type_pair == 'HH':
            return np.array([[p['ss_sigma']]])
        
        if type_pair == 'YH':
            v_ss = p['ss_sigma']
            v_sd = p['sd_sigma']
            block = np.zeros((6, 1))
            block[0, 0] = v_ss
            block[1, 0] = v_sd * (n**2 - 0.5*(l**2 + m**2))
            block[2, 0] = v_sd * (np.sqrt(3)/2) * (l**2 - m**2)
            block[3, 0] = v_sd * np.sqrt(3) * l * m
            block[4, 0] = v_sd * np.sqrt(3) * m * n
            block[5, 0] = v_sd * np.sqrt(3) * l * n
            return block

        if type_pair == 'YY':
            v_ss = p['ss_sigma']
            v_dd = p['dd_sigma']
            block = np.zeros((6, 6))
            block[0, 0] = v_ss
            block[1, 1] = v_dd * (n**2 - 0.5*(l**2 + m**2))**2
            return block
        
        return np.zeros((1, 1))

    def find_neighbors(self, cutoff=4.0):
        neighbors = []
        for dx in [-1, 0, 1]:
            for dy in [-1, 0, 1]:
                for dz in [-1, 0, 1]:
                    R_offset = np.array([dx, dy, dz])
                    for i in range(len(self.all_pos)):
                        for j in range(len(self.all_pos)):
                            if dx == 0 and dy == 0 and dz == 0 and i == j:
                                continue
                            
                            diff = self.all_pos[j] + R_offset - self.all_pos[i]
                            dist_vec = np.dot(diff, self.lattice)
                            dist = np.linalg.norm(dist_vec)
                            
                            if 0.1 < dist < cutoff:
                                l, m, n = dist_vec / dist
                                if i < self.n_y and j < self.n_y:
                                    pair = 'YY'; p = self.params['YY']
                                elif i < self.n_y and j >= self.n_y:
                                    pair = 'YH'; p = self.params['YH']
                                elif i >= self.n_y and j < self.n_y:
                                    pair = 'HY'; p = self.params['YH']
                                else:
                                    pair = 'HH'; p = self.params['HH']
                                
                                neighbors.append({
                                    'i': i, 'j': j, 'R': R_offset, 
                                    'lmn': (l, m, n), 'pair': pair, 'p': p
                                })
        return neighbors

    def build_hamiltonian(self, k, neighbors):
        H = self.H_onsite.astype(complex)
        for nb in neighbors:
            i, j = nb['i'], nb['j']
            l, m, n = nb['lmn']
            pair = nb['pair']
            p = nb['p']
            phase = np.exp(1j * 2 * np.pi * np.dot(k, nb['R']))
            
            if pair == 'YY':
                block = self.get_sk_block(l, m, n, 'YY', p)
                H[i*6:(i+1)*6, j*6:(j+1)*6] += block * phase
            elif pair == 'YH':
                block = self.get_sk_block(l, m, n, 'YH', p)
                H[i*6:(i+1)*6, self.n_y*6 + (j-self.n_y)] += block.flatten() * phase
            elif pair == 'HY':
                block_yh = self.get_sk_block(-l, -m, -n, 'YH', p)
                H[self.n_y*6 + (i-self.n_y), j*6:(j+1)*6] += block_yh.flatten().conj() * phase
            elif pair == 'HH':
                block = self.get_sk_block(l, m, n, 'HH', p)
                H[self.n_y*6 + (i-self.n_y), self.n_y*6 + (j-self.n_y)] += block[0, 0] * phase
                
        return H

def calculate_spectral_weights(k, eigvecs):
    """
    Calculates the unfolding weights for BCC supercell to primitive cell.
    Translation vector T = [0.5, 0.5, 0.5]
    """
    T = np.array([0.5, 0.5, 0.5])
    phase = np.exp(1j * 2 * np.pi * np.dot(k, T))
    
    # Mapping pairs of orbitals (Unit 0, Unit 1)
    # Y orbitals: (0-5) and (6-11)
    # H orbitals: (12-17) and (18-23)
    unit0_indices = list(range(0, 6)) + list(range(12, 18))
    unit1_indices = list(range(6, 12)) + list(range(18, 24))
    
    weights = []
    for n in range(eigvecs.shape[1]):
        v = eigvecs[:, n]
        # Sum over primitive cell orbital basis
        w_sum = 0.0
        for i0, i1 in zip(unit0_indices, unit1_indices):
            c0 = v[i0]
            c1 = v[i1]
            w_sum += np.abs(c0 + c1 * phase)**2
        weights.append(w_sum / 2.0)
    return np.array(weights)

def main():
    a0 = 3.3698
    y_pos = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]
    h_pos = [
        [0.25, 0.00, 0.50], [0.25, 0.50, 0.00], [0.50, 0.25, 0.00],
        [0.00, 0.25, 0.50], [0.00, 0.50, 0.25], [0.50, 0.00, 0.25],
        [0.75, 0.50, 0.00], [0.75, 0.00, 0.50], [0.00, 0.75, 0.50],
        [0.50, 0.75, 0.00], [0.50, 0.00, 0.75], [0.00, 0.50, 0.75]
    ]

    params = {
        'e_y_s': -1.5, 'e_y_d': 2.0, 'e_h_s': -4.5,
        'YY': {'ss_sigma': -0.5, 'dd_sigma': -0.3},
        'YH': {'ss_sigma': -1.8, 'sd_sigma': -1.2},
        'HH': {'ss_sigma': -0.8}
    }

    model = YH6TightBinding(a0, y_pos, h_pos, params)
    neighbors = model.find_neighbors(cutoff=3.5)

    # K-path: G-X-M-G-R-X
    k_points = [
        [0, 0, 0], [0, 0.5, 0], [0.5, 0.5, 0], [0, 0, 0], [0.5, 0.5, 0.5], [0, 0.5, 0]
    ]
    labels = [r'$\Gamma$', 'X', 'M', r'$\Gamma$', 'R', 'X']
    
    n_steps = 100
    full_path = []
    tick_indices = [0]
    for i in range(len(k_points)-1):
        for t in np.linspace(0, 1, n_steps, endpoint=False):
            full_path.append(np.array(k_points[i])*(1-t) + np.array(k_points[i+1])*t)
        tick_indices.append(len(full_path))
    full_path.append(k_points[-1])

    all_energies = []
    all_weights = []
    for k in full_path:
        H = model.build_hamiltonian(k, neighbors)
        H = (H + H.conj().T) / 2.0
        evals, evecs = eigh(H)
        all_energies.append(evals)
        weights = calculate_spectral_weights(k, evecs)
        all_weights.append(weights)
    
    all_energies = np.array(all_energies)
    all_weights = np.array(all_weights)

    # Plot
    plt.figure(figsize=(10, 8))
    
    # Flatten data for scatter plot
    x_coords = np.repeat(np.arange(len(full_path)), model.dim)
    y_coords = all_energies.flatten()
    weights = all_weights.flatten()
    
    # Use scatter plot with weights for color intensity and size
    plt.scatter(x_coords, y_coords, s=weights*15, c='tab:blue', alpha=weights, edgecolors='none')
    
    plt.ylim(-10, 10)
    for tick in tick_indices:
        plt.axvline(tick, color='gray', linestyle='--', alpha=0.5)
    
    plt.xticks(tick_indices, labels, fontsize=12)
    plt.ylabel("Energy (eV)", fontsize=12)
    plt.title("YH6 Unfolded Band Structure (Primitive Cell View)", fontsize=14)
    plt.grid(axis='y', alpha=0.2)
    plt.tight_layout()
    plt.savefig('YH6_Unfolded_Bands.png')
    print("Unfolded band plot saved to YH6_Unfolded_Bands.png")

if __name__ == "__main__":
    main()