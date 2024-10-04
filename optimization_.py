import itertools
import numpy as np


def optimization_(df_in, columns, num_drugs: int = 2):
    # Define Heatmap matrix

    Heatmap = df_in.to_numpy()
    # Replace negative values with -1, positive values with 1, and leave 0 as is:
    D = np.sign(Heatmap)

    # Obtain D length
    K = D.shape[0]

    # Combinations/Permutations
    # Original vector with elements from 1 to 24 drugs.
    original_vector = np.arange(1, K + 1)
    # Number of drugs to choose N <= K.
    N = num_drugs

    # Calculate all possible combinations of N elements from the original vector of K elements.
    combinations = list(itertools.combinations(original_vector, N))

    K1 = len(combinations)

    # WEIGHTS, COST, and REFERENCE
    alpha_cr = 1
    alpha_wt = 1
    alpha_cs = 1
    Ref_WT = 0
    Ref_CS = 1
    Ref_CR = 0
    Costo = np.zeros(K1)
    Costok = np.zeros(N)

    for i in range(K1):
        d = combinations[i]
        N = len(d)
        T = np.zeros((N, N))
        for k in range(N):
            for l in range(N):
                T[k, l] = D[d[k]-1, d[l]-1]

        for k in range(N):
            WT = np.sum(T[k, :] == 0) / N
            CS = np.sum(T[k, :] == -1) / N
            CR = np.sum(T[k, :] == 1) / N
            Costok[k] = alpha_wt * (Ref_WT - WT) ** 2 + alpha_cs * abs(Ref_CS - CS) ** 2 + alpha_cr * abs(Ref_CR - CR) ** 2

        Costo[i] = np.sum(Costok)

    # Solution
    # Find the minimum value in the vector.
    min_value = np.min(Costo)
    # Find all positions (indices) where the minimum value is found (is possible there are more than one solution)
    indices = np.where(Costo == min_value)[0]
    CostoTotal = Costo[indices[0]] / N

    # Assign names of drugs to the numbers from 1 to 24.
    nombres = columns
    # Vector c for the optimal combinations (there may be multiple coincidences, we choose the first one)
    c = combinations[indices[0]]
    #  Get the drug names associated with them.
    nombres_asociados = [nombres[i] for i in c]

    return nombres_asociados, c
