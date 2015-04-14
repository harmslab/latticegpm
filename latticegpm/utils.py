import itertools as it
# ------------------------------------------------------------
# Jesse Bloom's Lattice Model imports
# ------------------------------------------------------------
from latticeproteins.sequences import HammingDistance, RandomSequence

def search_fitness_landscape(Fitness, threshold, max_iter=1000):
    """ Randomly search the Fitness landscape for two sequences that 
        have non-zero fitnesses above threshold and differ at all sites. 
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

def generate_binary_space(wildtype, mutant):
    """ Generate binary genotype space between two sequences (that should differ at all sites) """
    if len(wildtype) != len(mutant):
        raise IndexError("ancestor_sequence and derived sequence must be the same length.")

    binaries = sorted(["".join(list(s)) for s in it.product('01', repeat=len(wildtype))])
    sequence_space = list()
    for b in binaries:
        binary = list(b)
        sequence = list()
        for i in range(len(wildtype)):
            if b[i] == '0':
                sequence.append(wildtype[i])
            else:
                sequence.append(mutant[i])
        sequence_space.append(''.join(sequence))
    return sequence_space
