from __future__ import print_function
import argparse
import copy
import nltk
import threading
import time as t
import numpy as np
from rdkit import Chem
from rdkit import rdBase

import cfg_util
import score_util
import zinc_grammar

rdBase.DisableLog('rdApp.error')
GCFG = zinc_grammar.GCFG


def CFGtoGene(prod_rules, max_len=-1):
    gene = []
    for r in prod_rules:
        lhs = GCFG.productions()[r].lhs()
        possible_rules = [idx for idx, rule in enumerate(GCFG.productions())
                          if rule.lhs() == lhs]
        gene.append(possible_rules.index(r))
    if max_len > 0:
        if len(gene) > max_len:
            gene = gene[:max_len]
        else:
            gene = gene + [np.random.randint(0, 256)
                           for _ in range(max_len-len(gene))]
    return gene


def GenetoCFG(gene):
    prod_rules = []
    stack = [GCFG.productions()[0].lhs()]
    for g in gene:
        try:
            lhs = stack.pop()
        except Exception:
            break
        possible_rules = [idx for idx, rule in enumerate(GCFG.productions())
                          if rule.lhs() == lhs]
        #print(possible_rules)
        rule = possible_rules[g % len(possible_rules)]
        prod_rules.append(rule)
        rhs = filter(lambda a: (type(a) == nltk.grammar.Nonterminal)
                     and (str(a) != 'None'),
                     zinc_grammar.GCFG.productions()[rule].rhs())
        stack.extend(list(rhs)[::-1])
    return prod_rules


def selectParent(population, tournament_size=3):
    idx = np.random.randint(len(population), size=tournament_size)
    best = population[idx[0]]
    for i in idx[1:]:
        if population[i][0] > best[0]:
            best = population[i]
    return best


def mutation(gene):
    idx = np.random.choice(len(gene))
    gene_mutant = copy.deepcopy(gene)
    gene_mutant[idx] = np.random.randint(0, 256)
    return gene_mutant


def canonicalize(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if smiles != '' and mol is not None and mol.GetNumAtoms() > 1:
        return Chem.MolToSmiles(mol)
    else:
        return smiles


elapsed_min = 0
best_score = 0
mean_score = 0
std_score = 0
min_score = 0
best_smiles = ""
all_smiles = []


def current_best():
    global elapsed_min
    global best_score
    global best_smiles
    global mean_score
    global min_score
    global std_score
    global all_smiles
    elapsed_min += 1
    print("${},{},{},{}"
          .format(elapsed_min, best_score, best_smiles, len(all_smiles)))
    t = threading.Timer(60, current_best, [])
    t.start()

def migration(Pipes, island_id, nb_of_island, population, migration_nb):
    # Pipes is an array containing every pipe between every islands
    # island_id, the number of the island
    # population the current population (score,smile,gene)
    # migration_nb the migration number, e.g. the k-th migration can be 1 to len(island)-1

    k = migration_nb #between 1 and number of island
    island_list = [i for i in range(nb_of_island)]

    top5 = round(0.05*len(population))
    for i in island_list:
        if i==island_id:
            j = i-k
            (Pipes.iloc[i,j])[0].send(population[0:top5])
            print('send from ',island_id,' with i and j: ',i,j)
            # recv
            j += k
            i += k
            if i>=nb_of_island : 
                i = i - nb_of_island
            migrants = (Pipes.iloc[i,j])[1].recv()
            print('recv from ',island_id)
            population = population[:len(population)-top5] + migrants 
        else :
            pass

    return population


def main(Pipes, island_id, nb_of_island):
    #parser = argparse.ArgumentParser()
    #parser.add_argument('--smifile', default='250k_rndm_zinc_drugs_clean.smi')
    #parser.add_argument('--seed', type=int, default=t.time())
    #args = parser.parse_args()

    smifile = '250k_rndm_zinc_drugs_clean.smi'
    np.random.seed(int(t.time() + island_id))
    #np.random.seed(0)
    global best_smiles
    global best_score
    global all_smiles

    gene_length = 300

    N_mu = 100
    N_lambda = 200

    # initialize population
    seed_smiles = []
    with open(smifile) as f:
        for line in f:
            smiles = line.rstrip()
            seed_smiles.append(smiles)

    initial_smiles = np.random.choice(seed_smiles, N_mu+N_lambda)
    initial_smiles = [canonicalize(s) for s in initial_smiles]
    initial_genes = [CFGtoGene(cfg_util.encode(s), max_len=gene_length)
                     for s in initial_smiles]
    initial_scores = [score_util.calc_score(s) for s in initial_smiles]
    #print(initial_scores)
    population = []
    for score, gene, smiles in zip(initial_scores, initial_genes,
                                   initial_smiles):
        population.append((score, smiles, gene))

    population = sorted(population, key=lambda x: x[0], reverse=True)[:N_mu]

    th = threading.Timer(60, current_best, [])
    th.start()
    print("Start!")
    all_smiles = [p[1] for p in population]
    #print([p[0] for p in population])
    mig_interval = 1000 # A migration every 1000 iteration
    x = [ i for i in range(1000,1000000000,mig_interval)] # All the generation in wich a migration should occur
    k = 1 # First migration
    t0=t.time()
    for generation in range(1000000000):
        scores = [p[0] for p in population]
        mean_score = np.mean(scores)
        min_score = np.min(scores)
        std_score = np.std(scores)
        best_score = np.max(scores)
        idx = np.argmax(scores)
        best_smiles = population[idx][1]
        print("%{},{},{},{},{}".format(generation, best_score,
                                       mean_score, min_score, std_score))

        new_population = []
        for _ in range(N_lambda):
            p = population[np.random.randint(len(population))]
            p_gene = p[2]
            c_gene = mutation(p_gene)

            c_smiles = canonicalize(cfg_util.decode(GenetoCFG(c_gene)))
            if c_smiles not in all_smiles:
                c_score = score_util.calc_score(c_smiles)
                c = (c_score, c_smiles, c_gene)
                new_population.append(c)
                all_smiles.append(c_smiles)

        population.extend(new_population)
        population = sorted(population,
                            key=lambda x: x[0], reverse=True)[:N_mu]

        # Every mig_interval generation make
        if generation in x:
            print('Starting Migration')
            if k >= nb_of_island:
                k = 1
            population = migration(Pipes, island_id, nb_of_island, population, k)
            k+=1
        if t.time() - t0 >= 3600*8 :
            break
    f = open(str(island_id)+'_final_pop.csv','w')
    f.write(population)
    f.close()
if __name__ == "__main__":
    main()
