
def conformation_space(Conformations, temperature, threshold, target_conf=None, differby=None, max_iter=1000):
    """Randomly search the conformations landscape for two sequences that
    fold with stability below some threshold and differ at all sites.

    Parameters
    ----------
    Conformations : latticeproteins.conformations.Conformations object
        object that holds all the possible conformations to search.
    temperature : float
        Temperature parameter (ratio to kT)
    threshold : float
        Maximum allowed stability for landscape

    Returns
    -------
    sequences : list of two strings
        List of two sequences that differ at all sites and fold.
    """
    length = Conformations.Length()
    # Check
    if differby is None:
        differby = length
    elif differby > length:
        raise Exception("differby cannot be larger than the length of the sequences.")

    # Find a sequence that's below a certain energy.
    energy = threshold
    while energy >= threshold:
        sequence = RandomSequence(length)
        output = Conformations.FoldSequence(sequence, temperature)
        energy = output[0]

    # Set the resulting sequence as the first variable
    sequence1 = ''.join(sequence)
    energy = threshold # Reset energy

    # Search for next function
    counter = -1
    mutants = NMutants(sequence1, differby, max_iter)
    while energy >= threshold:
        counter += 1
        output = Conformations.FoldSequence(mutants[counter], temperature, target_conf=target_conf)
        energy = output[0]

        # Check looping
        if counter >= max_iter:
            raise Exception("Reached max iteration in search.")

    # Set second sequence
    sequence2 = ''.join(mutants[counter])
    return sequence1, sequence2


def search_fitness_landscape(Fitness, threshold, differby=None, max_iter=1000):
    """ Randomly search the Fitness landscape for two sequences that
    have non-zero fitnesses above threshold and differ at all sites.

    Parameters
    ----------
    Conformations: latticeproteins.conformations.Conformations object
        object that holds all the possible conformations to search.
    temperature: float
        Temperature parameter (ratio to kT)
    threshold: float
        Maximum allowed binding energy for landscape

    Returns
    -------
    sequences: list of two strings
        List of two sequences that differ at all sites and fold.
    """

    length = Fitness.Length()
    # Check
    if differby is None:
        differby = length
    elif differby > length:
        raise Exception("differby cannot be larger than the length of the sequences.")

    # Find a sequence that's below a certain fitness.
    fitness = threshold
    while fitness <= threshold:
        sequence = RandomSequence(length)
        fitness = Fitness.Fitness(sequence)

    # Set the resulting sequence as the first variable
    sequence1 = ''.join(sequence)
    fitness = threshold # Reset fitness

    # Search for next function
    counter = -1
    mutants = NMutants(sequence1, differby, max_iter)
    while fitness <= threshold:
        counter += 1
        fitness = Fitness.Fitness(mutants[counter])

        # Check looping
        if counter >= max_iter:
            raise Exception("Reached max iteration in search.")

    # Set second sequence
    sequence2 = ''.join(mutants[counter])
    return sequence1, sequence2
