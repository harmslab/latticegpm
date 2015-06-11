import itertools as it
import numpy as np
# ------------------------------------------------------------
# Jesse Bloom's Lattice Model imports
# ------------------------------------------------------------
from latticeproteins.sequences import HammingDistance, RandomSequence, NMutants
from latticeproteins.interactions import miyazawa_jernigan

class ConformationError(Exception):
    """ If protein doesn't fold, raise error. """

def compare_sequences(s1, s2):
    """ Return the indice where two strings differ. """
    return [i for i in range(len(s1)) if s1[i] != s2[i]]

def fold_energy(sequence, conformation):
    """ Calculate the energy of the sequence with the given conformation. 
    
        Args:
        ----
        sequence: str
            Amino acid sequence to fold.
        conformation: str
            Conformation according to latticemodel's conformations format (e.g. 'UDLLDRU')
        
        Returns:
        -------
        energy: float
            energy of the conformation (sum of all contact energies)
    """
    contacts = lattice_contacts(sequence, conformation)
    energy = sum([miyazawa_jernigan[c] for c in contacts])
    return energy
    
def lattice_contacts(sequence, conformation):
    """ Find all contacts in conformation.
    
        Args:
        ----
        sequence: str
            Amino acid sequence to fold.
        conformation: str
            Conformation according to latticemodel's conformations format (e.g. 'UDLLDRU')
        
        Returns:
        -------
        contacts: list
            list of contact pairs
    """
    sites = list(sequence)
    length = len(sites)
    
    try:
        moves = list(conformation)  
    except TypeError:
        raise ConformationError("""Protein conformation is None; is there a native state? """)
    
    # build a coordinate system, note that odd rotation of intuitive coordinates
    # since we are working in numpy array grid.
    coordinates = {"U": [-1,0], "D":[1,0], "L":[0,-1], "R":[0,1]}
    grid = np.zeros((length, length), dtype=str)
    x = y = round(length/2.0) # initial position on the grid is at the center of the 2d array
    grid[x,y] = sites[0]
    
    # move on grid, populate with amino acid at that site, and store all contacting neighbors. 
    contacts = []
    for i in range(length-1):
        step = coordinates[moves[i]]
        x += step[0]
        y += step[1]
        grid[x,y] = sites[i+1]
        neighbors = [sites[i+1] + grid[x+c[0], y+c[1]] for c in coordinates.values()]
        contacts += [n for n in neighbors if n in miyazawa_jernigan]

    # subtract the contacts that have bonds between them.    
    for i in range(1,length):
        try:
            contacts.remove(sequence[i-1:i+1])
        except ValueError:
            contacts.remove(sequence[i] + sequence[i-1])
    
    return contacts
    

def search_conformation_space(Conformations, temperature, threshold, target_conf=None, differby=None, max_iter=1000):
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
    rev_mutations = [mutations[i] for i in range(n_mut-1, -1, -1)]
    mutation_map = dict(zip(range(n_mut), mutations))
    
    # Build a binary representation
    combinations = [list(j) for i in range(1,n_mut+1) for j in it.combinations(rev_mutations, i)]
    sequence_space = [wildtype]
    for c in combinations:
        sequence = list(wildtype)
        for el in c:
            sequence[el] = mutant[el]
        sequence_space.append("".join(sequence))
     
    return sequence_space