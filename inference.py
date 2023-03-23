import pandas as pd
import  numpy as np

def get_favorable_pos():
	result = energy
	min_docking_energy_val = result.min()[2]
	best_result = result[result[2] == min_docking_energy_val]
	best_angle = best_result[1].values[0]
	best_idx  = best_result[0].values[0]
	return best_idx, best_angle, min_docking_energy_val

def calc_rmsd(actual,pred):
  mse = (actual - pred)**2
  mse = np.sum(mse, axis=1)
  mse = np.mean(mse)
  rmsd = np.sqrt(mse)
  return rmsd


if __name__ == '__main__':
	energy = pd.read_csv('./output/energy.pdb', delim_whitespace=True, header=None)
	lig_actual = pd.read_csv('./ligand_actual.pdb', delim_whitespace=True, header=None)
	lig_start = pd.read_csv('./ligand_starting.pdbm', delim_whitespace=True, header=None)
	lig_actual = lig_actual.head(12)
	lig_best = pd.read_csv(f'./output/ligand_{get_favorable_pos()[0]}.pdb', delim_whitespace=True, header=None)
	rmsd1 = calc_rmsd(lig_actual, lig_start.iloc[:12, 6:9] )
	rmsd2 = calc_rmsd(lig_actual, lig_best.iloc[:12, 6:9])
	print(f'RMSD for start and actual = {rmsd1}')
	print(f'RMSD for start and actual = {rmsd2}')