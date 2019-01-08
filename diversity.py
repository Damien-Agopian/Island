from rdkit import Chem, DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
import pandas as pd
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type=int, default=0, help='Number of islands')
    parser.add_argument('-m', type=int, help='Migration Interval')
    parser.add_argument('-a', type=str, default='island',
                        help='"island" or "chem" : Comput diversity for islands or chemGE')
    args = parser.parse_args()

    N_island = args.i
    alg = args.a
    m = args.m
    ms = []
    score = []

    if alg == 'island':
        for i in range(N_island):

            tab = pd.read_csv(str(i)+'_final_pop_'+str(N_island)+'_'+str(m)+'.csv')
            ms += list(tab.iloc[:,2])
            score += list(tab.iloc[:,1])
        fps = [FingerprintMols.FingerprintMol(Chem.MolFromSmiles(x)) for x in ms]

        t = 0.0
        for x in fps:
            for y in fps:
                t += 1-DataStructs.FingerprintSimilarity(x, y, metric=DataStructs.TanimotoSimilarity)
        print('diversity: ',t/(len(fps)*len(fps)))
        print('best score: ',max(score))    
        #print(score)
    if alg == 'chem':
        
        tab = pd.read_csv(str(N_island)+'_chemGE.csv')
        ms += list(tab.iloc[:,2])
        score += list(tab.iloc[:,1])
        fps = [FingerprintMols.FingerprintMol(Chem.MolFromSmiles(x)) for x in ms]

        t = 0.0
        for x in fps:
            for y in fps:
                t += 1-DataStructs.FingerprintSimilarity(x, y, metric=DataStructs.TanimotoSimilarity)
        print('diversity: ',t/(len(fps)*len(fps)))
        print('best score: ',max(score))

if __name__ == '__main__':
    main()
