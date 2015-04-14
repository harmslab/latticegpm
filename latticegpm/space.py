import numpy as np
from latticeproteins.conformations import Conformations, BindLigand, PrintConformation
from latticegpm.utils import generate_binary_space
from latticegpm.mapping import LatticeMap

# ------------------------------------------------------
# Build a binary protein lattice model sequence space
# with fitness defined by function in Jesse Blooms'
# latticeproteins python package.
# ------------------------------------------------------

class LatticeSequenceSpace(LatticeMap):
    
    def __init__(self, wildtype, mutant, Fitness):
        """ Build a protein lattice model sequence space from a given fitness function. """
        self.sequences = generate_binary_space(wildtype, mutant)
        self.wildtype = wildtype
        self.temperature = Fitness._temp
        # Check that the fitness landscape describes ligand binding.
        if Fitness._ligand is None:
            raise Exception("The Fitness object must be associated with ligand binding for this analysis. ")
        else:
            self.ligand_tup = Fitness._ligand
        
        # Fold proteins and extract stability parameter and native conformations
        folds = np.array([Fitness._conformations.FoldSequence(s, self.temperature) for s in self.sequences])
        self.stabilities = folds[:,0]
        self.conformations = folds[:,1]
        
        # Calculate binding energies and fitnesses (fitness = exp(-be))
        binding_energies = np.empty(self.n, dtype=float)
        for i in self._indices:
            if self.conformations[i] is None:
                binding_energies[i] = 0
            else:
                binding_energies[i] = BindLigand(self.sequences[i], self.conformations[i], self.ligand, self.ligand_conformation, numsaved=100000)[0]
        self.fitnesses = np.exp(-1*binding_energies)
                
    def print_sequences(self, sequences, with_ligand=False):
        """ Print sequence conformation with/without ligand bound. """
        seq2conf = self.seq2conformation
        for s in sequences:
            be, x,y, conf = BindLigand(s, seq2conf[s], self.ligand, self.ligand_conformation)
            ligand_stuff = (self.ligand, conf, x, y)
            PrintConformation(s, seq2conf[s], ligand_tup = ligand_stuff)