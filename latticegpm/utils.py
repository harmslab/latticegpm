import itertools as it
# ------------------------------------------------------------
# Jesse Bloom's Lattice Model imports
# ------------------------------------------------------------
from latticeproteins.sequences import HammingDistance, RandomSequence, NMutants

""" ADD NMUTANTS STUFF. """

def compare_sequences(s1, s2):
    """ Return the indice where two strings differ. """
    return [i for i in range(len(s1)) if s1[i] != s2[i]]

def search_conformation_space(Conformations, temperature, threshold, differby=None, max_iter=1000):
    """ Randomly search the conformations landscape for two sequences that 
        fold with energy below some threshold and differ at all sites. 
        
        Args:
        ----
        Conformations: latticeproteins.conformations.Conformations object
            object that holds all the possible conformations to search.
        temperature: float
            Temperature parameter (ratio to kT)
        threshold: float
            Maximum allowed binding energy for landscape
        
        Returns:
        -------
        sequences: list of two strings
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
        output = Conformations.FoldSequence(mutants[counter], temperature)
        energy = output[0]
    
    # Check looping
    if counter >= max_iter:
        raise Exception("Reach max iteration in search.")
    
    # Set second sequence
    sequence2 = ''.join(mutants[counter])
    return sequence1, sequence2


def search_fitness_landscape(Fitness, threshold, max_iter=1000):
    """ Randomly search the Fitness landscape for two sequences that 
        have non-zero fitnesses above threshold and differ at all sites. 
        
        Args:
        ----
        Conformations: latticeproteins.conformations.Conformations object
            object that holds all the possible conformations to search.
        temperature: float
            Temperature parameter (ratio to kT)
        threshold: float
            Maximum allowed binding energy for landscape
        
        Returns:
        -------
        sequences: list of two strings
            List of two sequences that differ at all sites and fold.
    """
    length = Fitness.Length()
    counter = 0
    sequences = list()
    while len(sequences) < 2 and counter < max_iter:
        sequence = RandomSequence(length)
        fitness = Fitness.Fitness(sequence)
        # Check that fitness value is above the threshold
        if fitness > threshold:
            # Check Hamming distance once sequences list contains more than 2 sequences.
            if len(sequences) > 0:
                if HammingDistance(sequences[0], sequence) is length:
                    sequences.append("".join(sequence))
            else:
                sequences.append("".join(sequence))   
        counter +=1
    
    # Raise error if random search reaches maximum iterations
    if counter == max_iter:
        raise Exception("Random search reached max iterations before finding satisfying sequences.")
        
    return sequences

def enumerate_space(wildtype, mutant):
    """ Build a list of sequences that represent all binary combinations of the wildtype 
        and mutant sequence. Systematically flips all sites in the wildtype sequence, 
        mutating towards the mutant.
    
        Args:
        ----
        wildtype: str
            Wildtype sequence
        mutant: str
            Mutant sequence that differs at all sites from wildtype
    """
    
    # Check that wildtype and mutant are the same length
    if len(wildtype) != len(mutant):
        raise IndexError("ancestor_sequence and derived sequence must be the same length.")
    
    # Count mutations and keep indices
    mutations = compare_sequences(wildtype, mutant)
    n_mut = len(mutations)
    mutation_map = dict(zip(range(n_mut), mutations))
    
    # Build a binary representation
    combinations = [list(j) for i in range(1,n_mut+1) for j in it.combinations(mutations, i)]
    sequence_space = [wildtype]
    for c in combinations:
        sequence = list(wildtype)
        for el in c:
            sequence[-el-1] = mutant[-el-1]
        sequence_space.append("".join(sequence))
    return sequence_space