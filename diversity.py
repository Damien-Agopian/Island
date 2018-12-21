from rdkit import Chem, DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
import pandas as pd
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type=int, default=0, help='Number of islands')
    parser.add_argument('-a', type=str, default='island',
                        help='"island" or "chem" : Comput diversity for islands or chemGE')
    args = parser.parse_args()

    N_island = args.i
    alg = args.a
    
    ms = []
    if alg == 'island':
        for i in range(N_island):

            tab = pd.read_csv(str(i)+'_final_pop.csv')
            ms += list(tab.iloc[:,2])

        fps = [FingerprintMols.FingerprintMol(Chem.MolFromSmiles(x)) for x in ms]

        t = 0.0
        for x in fps:
            for y in fps:
                t += 1-DataStructs.FingerprintSimilarity(x, y, metric=DataStructs.TanimotoSimilarity)
        print(t/(len(fps)*len(fps)))
    
    if alg == 'chem':
        tab = pd.read_csv('chemGE.csv')
        ms += list(tab.iloc[:,2])

        fps = [FingerprintMols.FingerprintMol(Chem.MolFromSmiles(x)) for x in ms]

        t = 0.0
        for x in fps:
            for y in fps:
                t += 1-DataStructs.FingerprintSimilarity(x, y, metric=DataStructs.TanimotoSimilarity)
        print(t/(len(fps)*len(fps)))


if __name__ == '__main__':
    main()
