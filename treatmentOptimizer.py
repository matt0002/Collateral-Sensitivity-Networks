import numpy as np
from scipy.integrate import odeint


def objective_function(model, F0, times, R, M, timeTreatmentGap, sequence, bactTypes, antibiotics, K, death_rate, mu, persister_pop):
    # Solving the differential equations
    fSol = odeint(model, F0, times, args=(R, M, timeTreatmentGap, sequence, antibiotics, K, death_rate, mu, persister_pop))

    popWeights = np.ones(bactTypes)

    X = np.dot(popWeights, fSol[-1, :])

    return X


def optimizedPopulation(model, maxGen, popNumber, bounds, timeTreatmentGap, bactTypes, F0, times, R, M,
                        antibiotics, K, death_rate, mu, persister_pop, F=0.5, Cr=0.7, minError=0.0):
    Np = popNumber
    D = int(times[-1] // timeTreatmentGap)
    lowBound, uppBound = bounds

    print(times)

    # Initialize the population with random solutions
    population = np.random.randint(lowBound, uppBound, [Np, D])
    fitness = [
        objective_function(model, F0, times, R, M, timeTreatmentGap, i, bactTypes, antibiotics, K, death_rate, mu, persister_pop) for i
        in population]
    oldFitnessMin = fitness[np.argmin(fitness)]

    # Main optimization loop
    for generation in range(maxGen):
        print(generation)
        for i in range(Np):
            # Create the indices population without index i
            indices = np.delete(np.arange(Np), [i])

            # Select three distinct indices
            r_indeces = np.random.choice(indices, size=3, replace=False)

            # Extract the three random vectors from the population
            random_vecs = population[r_indeces]

            # Create a mutant vector
            mutant = random_vecs[0] + F * (random_vecs[1] - random_vecs[2])  # DE/rand/1

            # Ensure the values of mutant are within the search space
            mutant = np.clip(mutant, lowBound, uppBound - 1)

            # Perform crossover to create a trial solution
            crossover_mask = np.random.rand(D) < Cr
            if crossover_mask.sum() == 0:
                crossover_mask[np.random.randint(D)] = True
            trial = np.where(crossover_mask, mutant, population[i])

            # Ensure integer values
            trial = np.round(trial).astype(int)

            # Evaluate the fitness of the trial vector
            trial_fitness = objective_function(model, F0, times, R, M, timeTreatmentGap, trial, bactTypes,
                                               antibiotics, K, death_rate, mu, persister_pop)

            # Update the population if the trial vector is better
            if trial_fitness < fitness[i]:
                population[i] = trial
                fitness[i] = trial_fitness

        if minError > 0.0:
            newFitnessMin = fitness[np.argmin(fitness)]
            error = np.abs(newFitnessMin - oldFitnessMin)
            if error < minError:
                print(f"Stopped in the generation {generation}")
                return population, fitness
            oldFitnessMin = newFitnessMin

    return population, fitness
