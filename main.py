import os

import numpy as np
import pandas as pd


def euclidean_distance(protein_coor, ligand_coor):
	return np.linalg.norm(protein_coor - ligand_coor)


def get_v_lj(c12, c6, r):
	return (c12 / r ** 12) - (c6 / r ** 6)


def get_eps_r(r):
	A = 6.02944
	B = 72.37056
	lam = 0.018733345
	kappa = 213.5782
	if (r < 1.32):
		return 8
	return (A + (B / (1 + kappa * np.exp(-1 * lam * B * r))))


def get_v_el(q_i, q_j, r_ij):
	k = 138.935485
	return (k * q_i * q_j) / (get_eps_r(r_ij) * r_ij)


def calc_v(ligand):
	v_all = 0
	for l_idx, lig_atom in enumerate(ligand.iterrows()):
		lig_row = ligand.iloc[l_idx]
		l_nm = lig_row[2]
		cx, cy, cz = lig_row[6], lig_row[7], lig_row[8]
		lig_coord = np.matrix([cx, cy, cz])
		lig_q = lig_row[9]
		# print(lig_row)
		for p_idx, pro_atom in enumerate(protein.iterrows()):
			p_row = protein.iloc[p_idx]
			p_nm = p_row[2]
			px, py, pz = p_row[6], p_row[7], p_row[8]
			p_coord = np.matrix([px, py, pz])
			p_q = p_row[9]
			r = euclidean_distance(p_coord, lig_coord)
			try:
				# param_row =  params.query(f'(@params[0]=="{l_nm}" & @params[1]=="{p_nm}") | (@params[0]=="{p_nm}" & @params[1]=="{l_nm}")')
				param_row = params.loc[
					((params[0] == l_nm) & (params[1] == p_nm) | (params[0] == p_nm) & (params[1] == l_nm))]
			except KeyError:
				# print("nope")
				pass
			c6 = float(param_row[2].values)
			c12 = float(param_row[3].values)
			# print(l_nm, p_nm, r, c6, c12)
			v_lj = get_v_lj(c12, c6, r / 10)
			v_el = get_v_el(lig_q, p_q, r)
			v_tot = v_lj + v_el
			v_all += v_tot
	# print(l_nm, p_nm, v_lj, v_el)
	return v_all


def rotation(beta, ligand_mat):
	lig_x = ligand_mat[6]
	lig_y = ligand_mat[7]
	lig_z = ligand_mat[8]
	ligand_coord = np.matrix([lig_x, lig_y, lig_z])
	"""
	Rotation along y-axis matrix
	[[cos(beta),0,sin(beta)
	    0,      1,      0
	 -sin(beta),0,cos(beta)]]
	"""
	rot_mat = np.matrix([[np.cos(beta), 0, np.sin(beta)], [0, 1, 0], [-np.sin(beta), 0, np.cos(beta)]]) # Rotation along y-axis matrix
	# print(rot_mat)
	rotated = np.matmul(rot_mat, ligand_coord)
	# convert back to matrix
	rot_tran = pd.DataFrame(rotated.transpose())
	ligand_mat[6], ligand_mat[7], ligand_mat[8] = rot_tran[0], rot_tran[1], rot_tran[2]

	return rot_mat, ligand_mat


if __name__ == '__main__':
	ligand = pd.read_csv('./ligand_starting.pdbm', delim_whitespace=True, header=None)
	protein = pd.read_csv('./protein.pdbm', delim_whitespace=True, header=None)
	params = pd.read_csv('./ffG43b1nb.params', delim_whitespace=True, header=None)
	increment = 100
	angles = [float(x) for x in np.linspace(0, 2 * np.pi, increment)]
	V = []
	rot_mats = []
	lig_coords = []
	betas = []
	for idx, angle in enumerate(angles):
		rot_mat, rotated_lig_mat = rotation(angle, ligand)
		v_sum = calc_v(rotated_lig_mat)
		print(f'{idx}---->{v_sum}')
		V.append(v_sum)
		rot_mats.append(rot_mat)
		lig_coords.append(rotated_lig_mat)
		betas.append(angle)
		if not os.path.exists('./output'):
			os.makedirs('output')
		ligand_file = open(f'./output/ligand_{idx}.pdb', 'a')
		ligand_file.write(rotated_lig_mat.to_csv(header=False, index=False, sep='\t'))
		energy_file = open('./output/energy.pdb', 'a')
		energy_file.write(f'{idx}\t{angle}\t{v_sum}\n')

	# print(list(zip(rot_mats, lig_coords, V)))
	result = zip(betas, V)
	print(min(result, key=lambda x: x[1]))
