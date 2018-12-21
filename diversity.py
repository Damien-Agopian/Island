from rdkit import Chem, DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
import pandas as pd


def main(N_island):
	ms = []

	for i in N_island:

		tab = pd.read_csv(str(i)+'_final_pop.csv')
		ms.append(list(tab.iloc[:,2]))

	fps = [FingerprintMols.FingerprintMol(x) for x in ms]

	t = 0.0
	for x in fps:
	    for y in fps:
	        t += 1-DataStructs.FingerprintSimilarity(x, y, metric=DataStructs.TanimotoSimilarity)
	print(t/(len(fps)*len(fps)))

if __name__ == '__main__':
	main()