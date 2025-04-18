import numpy as np
from scipy.integrate import odeint, solve_ivp
from scipy.optimize import differential_evolution
import plotly.graph_objects as go
import treatmentOptimizer


def set_up_diff_evo(final_graphs, drug1, drug2, drug3, drug_4, mutant_growth, susceptible_growth, death_rate,
                    cycle_time, K, mu, initial, num_cycles, persister_pop):
    df = final_graphs['DataFrame']
    columns = list(df.columns)
    df = df.to_numpy()
    df_ = np.sign(df)
    drugs = [drug1, drug2, drug3, drug_4]
    drugs_final = []
    [drugs_final.append(x) for x in drugs if x is not None]
    figure = set_up_ODE(drugs_final, df_, columns, mutant_growth, susceptible_growth, death_rate,
                        cycle_time, K, mu, initial, num_cycles, persister_pop)
    return figure


def matrix_associations(drugs_final, df, all_drugs):
    nombres_de_las_drogas = all_drugs

    # Replace negative values with -1, positive values with 1, and leave 0s as 0
    D = df  # profile matrix
    print(D)

    # Treatment design
    Entrada = drugs_final
    drogas = np.zeros(len(Entrada), dtype=int)
    # Find each drug in the vector of drug names and get their locations
    for i, droga in enumerate(Entrada):
        drogas[i] = nombres_de_las_drogas.index(droga)
    # drogas = [2, 22, 14]  # Choose drugs for treatment.
    nombres_seleccionados = [nombres_de_las_drogas[i] for i in drogas]
    #
    N = len(drogas)  # Number of drugs chosen for treatment (N=15 is the maximum to handle matrices of 2^15)
    T = np.zeros((N, N))
    for k in range(N):
        for l in range(N):
            T[k, l] = D[drogas[k], drogas[l]]

    # Creation of mij, M = Matrices_M = [M1, M2, M3] where Mi is the mutation matrix for drug i
    #######################################################################
    # Calculate the total number of possible states (2^k)
    num_estados = 2 ** N

    # States = 000, 001, ..
    #######################################
    # Initialize a matrix to store the states
    estados = np.zeros((num_estados, N))

    # Generate the states
    for i in range(num_estados):
        binario = np.binary_repr(i, width=N)  # Convert the decimal number to binary
        estados[i, :] = np.array(list(binario), dtype=int)  # Convert the binary string into a vector of numbers

    # Initialize a list to store the states
    estadosX = []

    # Generate the states
    for i in range(num_estados):
        binario = np.binary_repr(i, width=N)  # Convert the decimal number to binary
        estadoX = ''
        for bit in binario:
            if bit == '0':
                estadoX += 'S'
            else:
                estadoX += 'R'
        estadosX.append(estadoX)

    def operacion(D, fila_D, estado):
        k = len(estado)  # Get the length of the state
        resultado = ''  # Initialize the result as an empty string

        # Perform the operation with the D row according to the defined rules
        for j in range(k):
            if estado[j] == 'R':  # If the current element is 'R'
                # Perform the operation defined for 'R'
                if fila_D[j] == 1:
                    resultado += 'R'  # Concatenate 'R' to the result
                elif fila_D[j] == -1:
                    resultado += 'S'  # Concatenate 'S' to the result
                else:
                    resultado += estado[j]  # Keep the same character if D is 0
            elif estado[j] == 'S':  # If the current element is 'S'
                # Perform the operation defined for 'S'
                if fila_D[j] == 1:
                    resultado += 'R'  # Concatenate 'R' to the result
                elif fila_D[j] == -1:
                    resultado += 'S'  # Concatenate 'S' to the result
                else:
                    resultado += estado[j]  # Keep the same character if D is 0
            else:
                resultado += estado[j]  # Keep the same character if it's neither 'R' nor 'S'

        return resultado

    # Using operations like R:CS=S to calculate Mij with the heatmap T:
    ##########################################################################
    # Mutation matrix for all drugs in T
    # Initialize a list to store the adjacency matrices M1, M2, ..., MN
    Matrices_M = []

    # Calculate the adjacency matrices M1, M2, ..., MN
    for idx in range(T.shape[0]):
        # Initialize the adjacency matrix of the graph Mi
        Mi = np.zeros((num_estados, num_estados))

        # Calculate operacion(T, fila_i, estado) for all possible states
        for i in range(num_estados):
            for j in range(num_estados):
                # Calculate operacion(T, fila_i, estado) using row i of T
                xj = operacion(T, T[idx, :], estadosX[j])

                # If xj is equal to state i, then there is an edge
                if xj == estadosX[i]:
                    Mi[i, j] = 1

        # Store the adjacency matrix Mi in the corresponding list
        Matrices_M.append(Mi)

    return np.array(Matrices_M), estadosX, Entrada


def set_up_ODE(drugs_final, df, all_drugs, mutant_growth, susceptible_growth, death_rate, cycle_time, K, mu, initial, num_cycles, persister_pop):
    # select drugs

    state_associations, states, drugs_considered = matrix_associations(drugs_final, df, all_drugs)

    # growth rates
    ro_Growth_sensitive = susceptible_growth
    ro_Growth_Resistant = mutant_growth
    carrying_capacity = K
    all_drug_matrix = []
    drugs_in_states = len(states[0])
    for state_index in range(drugs_in_states):
        matrix_for_growth = []
        for state_position, state in enumerate(states):
            growth_matrix = [0 for i in range(len(states))]
            if state[state_index] == 'S':
                growth_matrix[state_position] = ro_Growth_sensitive
            else:
                growth_matrix[state_position] = ro_Growth_Resistant

            matrix_for_growth.append(growth_matrix)
        matrix_for_growth_array = np.array(matrix_for_growth)
        all_drug_matrix.append(matrix_for_growth)

    growth_array = all_drug_matrix
    optimal_graph, sequence_graph = diff_evo(growth_array, carrying_capacity, states, drugs_considered,
                                             state_associations, death_rate, cycle_time, mu, initial, num_cycles, persister_pop)

    return optimal_graph, sequence_graph


def diff_evo(matrix_of_growth_dynamics, K, states, drugs_considered, state_associations, death_rate,
             cycles, mu, initial, num_cycles, persister_pop):
    x_i = [0 for i in range(len(states))]
    x_i[0] = initial
    number_of_cycles = num_cycles
    dTsim = 0.1
    t_cycle = cycles
    antibiotics = len(drugs_considered)
    # growth matrix
    R = [np.array(x) for x in matrix_of_growth_dynamics]
    # State associations for mutation matrix
    M = [np.array(y) for y in state_associations]
    Tmax = number_of_cycles * cycles * len(drugs_considered)

    Delta_t = cycles
    times = np.linspace(0, Tmax, int(Tmax / dTsim) + 1)

    treatments = Tmax // Delta_t

    lower_bound = 0
    upper_bound = len(drugs_considered)
    bounds = (lower_bound, upper_bound)

    bactTypes = len(states)
    """result = differential_evolution(objective_function_int, bounds=bounds,
                                args=(x_i, times, R, M, Delta_t, bactTypes, antibiotics, K, death_rate, mu),
                                    strategy='best1bin', mutation=(0.5, 1), recombination=0.7,
                                    tol=0.01, popsize=10, maxiter=50)"""

    bestPop, fitness = treatmentOptimizer.optimizedPopulation(system_dynamics_int, maxGen=100, popNumber=30,
                                                              bounds=bounds,
                                                              timeTreatmentGap=cycles, bactTypes=bactTypes, F0=x_i,
                                                              times=times, R=R, M=M, antibiotics=antibiotics, K=K,
                                                              death_rate=death_rate, mu=mu, persister_pop=persister_pop)
    indBest = np.argmin(fitness)
    best_sequence = bestPop[indBest]

    optimal_graph, sequence_graph = run_optimal_(matrix_of_growth_dynamics, K, states, drugs_considered,
                                                 state_associations, death_rate,
                                                 Delta_t, mu, x_i, best_sequence, num_cycles, persister_pop)

    return optimal_graph, sequence_graph


def run_optimal_(matrix_of_growth_dynamics, K, states, drugs_considered, state_associations, death_rate,
                 cycles, mu, initial, best_sequence, num_cycles, persister_pop):
    dict_for_sim = {}
    t_cycle = cycles
    number_of_cycles = num_cycles



    states.append('Persister')
    x_i = initial
    t_eval = np.linspace(0, t_cycle, 10000)
    time_final = []
    cycle_data = []
    for index_, value in enumerate(x_i):
        dict_for_sim[f'Pop_{index_}'] = []
    dict_for_sim['persist'] = []
    for drug_index in best_sequence:  # drugs to iterate over
        update_x_i = []
        print(drug_index, drugs_considered[drug_index])
        # Pull matrix for drug i
        matrix_i = matrix_of_growth_dynamics[drug_index]

        # get mutation matrix for drug i
        mutation_state_i = state_associations[drug_index]

        # solve for given drug
        solution = solve_ivp(ODE_solve, y0=x_i, t_span=[0, t_cycle],
                             args=(matrix_i, K, death_rate, mu, mutation_state_i, persister_pop),
                             t_eval=t_eval)
        # Keep track of time
        time_final.append(list(solution.t))
        [cycle_data.append(drugs_considered[drug_index]) for i in range(len(list(solution.t)))]
        # update initial conditions
        for pop_index, pop_value in enumerate(solution.y):
            pop_values = dict_for_sim[f'Pop_{pop_index}']
            pop_values.append(list(solution.y[pop_index]))
            dict_for_sim[f'Pop_{pop_index}'] = pop_values
            update_x_i.append(list(solution.y[pop_index])[-1])

        x_i = update_x_i
        persist_pop = dict_for_sim['persist']
        persist_pop.append([persister_pop for x in range(len(solution.t))])

    # store all graph data into one list
    graph_data = []
    # concat each population into a single list for graphing
    t_final_len = [j for i in time_final for j in i]
    days = (number_of_cycles * t_cycle * len(drugs_considered))
    t_final = np.linspace(0, days, len(t_final_len))
    print(t_final)
    for index, key in enumerate(dict_for_sim.keys()):
        values = dict_for_sim[key]
        all_values = [j for i in values for j in i]
        if sum(all_values) == 0:
            pass
        else:
            graph_data.append(go.Scatter(x=t_final, y=all_values, name=f'{states[index]}', showlegend=True))

    # graph total population
    total_pop_values = []
    for val_ in dict_for_sim.keys():
        concat_values = [x for y in dict_for_sim[val_] for x in y]
        total_pop_values.append(concat_values)

    total_pop = []
    for i in range(len(t_final_len)):
        total_pop.append(sum([item[i] for item in total_pop_values]))

    sequence_data = [go.Scatter(x=t_final, y=cycle_data)]

    fig_sequence = go.Figure(sequence_data)
    fig_sequence.update_layout(
        yaxis_title='Drug Used',
        plot_bgcolor='white',
        xaxis_title='Time (hours)',
        title=f'Optimal Sequence')
    fig_sequence.update_layout()
    fig_sequence.update_yaxes(mirror=True,
                              ticks='outside',
                              showline=True,
                              linecolor='black',
                              gridcolor='white'
                              )
    fig_sequence.update_xaxes(mirror=True,
                              ticks='outside',
                              showline=True,
                              linecolor='black',
                              gridcolor='white')

    graph_data.append(go.Scatter(x=t_final, y=total_pop, name='Total Population'))
    fig = go.Figure(data=graph_data)

    fig.update_layout(
        yaxis_title='Bacteria',
        plot_bgcolor='white',
        xaxis_title='Time (hours)',
        title=f'Simulation with optimal sequence of antibiotics')
    fig.update_layout(yaxis_type='log')
    fig.update_yaxes(mirror=True,
                     ticks='outside',
                     showline=True,
                     linecolor='black',
                     gridcolor='white',
                     exponentformat='power',
                     range=[0, 10])
    fig.update_xaxes(mirror=True,
                     ticks='outside',
                     showline=True,
                     linecolor='black',
                     gridcolor='white')

    return fig, fig_sequence


def ODE_solve(t, y, final_matrix, K, ro_Death, mutation_rate, mutation_states, persister_pop):
    x_i = y
    x_i_multiply = np.c_[x_i]

    # update total population
    total_pop = sum(x_i)
    logistic = np.matrix([(1 - ((total_pop + persister_pop) / K)) for i in range(len(x_i))]).transpose()
    logistic = np.array(logistic)

    matrix_i = np.dot(final_matrix, logistic)

    # logistic growth vector
    matrix_i = matrix_i * x_i_multiply

    # subtract death vector from logistic growth vector
    death_matrix = np.matrix([ro_Death for i in range(len(x_i))]).transpose()
    death_matrix = np.array(death_matrix)

    final_matrix = matrix_i - (death_matrix * x_i_multiply)

    mutation_matrix = np.dot(mutation_states, mutation_rate)
    mutation_to_states = np.dot(mutation_matrix, x_i_multiply)

    dx_dt = final_matrix + mutation_to_states

    return dx_dt.transpose().tolist()[0]


def system_dynamics_int(F, t, R, M, Delta_t, sequence, antibiotics, K, death_rate, mu, persister_pop):
    idx = int(t / Delta_t)

    if idx >= len(sequence):
        idx = len(sequence) - 1  # Use the last index if t is greater than Tmax
    ind = int(sequence[idx])
    if ind == antibiotics:
        ind -= 1

    A = R[ind] * (1 - ((np.sum(F) + persister_pop) / K)) - death_rate * np.eye(len(F)) + mu * M[ind]
    dFdt = np.dot(A, F)
    return dFdt


def objective_function_int(sequence, F0, times, R, M, Delta_t, bactTypes, antibiotics, K, death_rate, mu, persister_pop):
    # Solving the differential equations
    fSol = odeint(system_dynamics_int, F0, times, args=(R, M, Delta_t, sequence, antibiotics, K, death_rate, mu, persister_pop))
    X = 0
    for i in range(bactTypes):
        X += fSol[-1, i]

    return X
