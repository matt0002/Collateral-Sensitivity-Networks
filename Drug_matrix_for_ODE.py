import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
import plotly.graph_objects as go
from differential_evo import diff_evo


def diff_eq(final_graphs, drug1, drug2, drug3, drug_4, mutant_growth, susceptible_growth, death_rate,
            cycle_time, K, mu, initial):
    df = final_graphs['DataFrame']
    columns = list(df.columns)
    df = df.to_numpy()
    df_ = np.sign(df)
    drugs = [drug1, drug2, drug3, drug_4]
    drugs_final = []
    [drugs_final.append(x) for x in drugs if x is not None]
    figure, fig_sequence = set_up_ODE(drugs_final, df_, columns, mutant_growth, susceptible_growth, death_rate,
                        cycle_time, K, mu, initial)
    return figure, fig_sequence


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


def ODE_for_Switching_systems(matrix_of_growth_dynamics, K, states, drugs_considered, state_associations, death_rate,
                              cycles, mu, initial):
    # set initial conditions
    x_i = [0 for i in range(len(states))]
    cycle_data = []
    x_i[0] = initial
    number_of_cycles = 6
    ro_Death = death_rate
    ro_Mutation = mu
    dict_for_sim = {}
    t_cycle = cycles
    t_eval = np.linspace(0, t_cycle, 10000)
    time_final = []
    for index_, value in enumerate(x_i):
        dict_for_sim[f'Pop_{index_}'] = []
    for cycle_drug in range(number_of_cycles):  # total simulation time
        for drug_index, drug_element in enumerate(drugs_considered):  # drugs to iterate over
            update_x_i = []
            # Pull matrix for drug i
            matrix_i = matrix_of_growth_dynamics[drug_index]

            # get mutation matrix for drug i
            mutation_state_i = state_associations[drug_index]

            # solve for given drug
            solution = solve_ivp(ODE_solve, y0=x_i, t_span=[0, t_cycle],
                                 args=(matrix_i, K, ro_Death, ro_Mutation, mutation_state_i),
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
        title=f'Sequential Treatment')
    fig_sequence.update_layout()
    fig_sequence.update_yaxes(mirror=True,
                              ticks='outside',
                              showline=True,
                              linecolor='black',
                              gridcolor='white')
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
        title=f'Simulation add drugs considered')
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


def ODE_solve(t, y, final_matrix, K, ro_Death, mutation_rate, mutation_states):
    x_i = y
    x_i_multiply = np.c_[x_i]

    # update total population
    total_pop = sum(x_i)
    logistic = np.matrix([(1 - (total_pop / K)) for i in range(len(x_i))]).transpose()
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


def set_up_ODE(drugs_final, df, all_drugs, mutant_growth, susceptible_growth, death_rate, cycle_time, K, mu, initial):
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
    graph_ode, fig_sequence = ODE_for_Switching_systems(growth_array, carrying_capacity, states, drugs_considered,
                                                        state_associations, death_rate, cycle_time, mu, initial)

    return graph_ode, fig_sequence
