import numpy as np
from latticeproteins.conformations import Conformations, BindLigand, PrintConformation
from latticegpm.utils import generate_binary_space
from latticegpm.mapping import LatticeMap

# ------------------------------------------------------
# Build a binary protein lattice model sequence space
# with fitness defined by function in Jesse Blooms'
# latticeproteins python package.
# ------------------------------------------------------

class LatticeConformationSpace(LatticeMap):
    
    def __init__(self, wildtype, mutant, Conformations, temperature=1.0):
        """ Build a protein lattice model sequence space from a conformation space. """
        self.sequences = generate_binary_space(wildtype, mutant)
        self.wildtype = wildtype
        self.temperature = temperature
        
        # Fold proteins and extract stability parameter and native conformations
        folds = np.array([Conformations.FoldSequence(s, self.temperature) for s in self.sequences])
        self.stabilities = np.array(folds[:,0], dtype=float)
        self.conformations = folds[:,1]
        
     def print_sequences(self, sequences):
         """ Print sequence conformation with/without ligand bound. """
         seq2conf = self.seq2conformation
         for s in sequences:
             PrintConformation(s, seq2conf[s])   
        

class LatticeFitnessSpace(LatticeFitnessMap):
    
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
        self.stabilities = np.array(folds[:,0], dtype=float)
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
        # Include ligand?
        if with_ligand:
            for s in sequences:
                be, x,y, conf = BindLigand(s, seq2conf[s], self.ligand, self.ligand_conformation)
                ligand_stuff = (self.ligand, conf, x, y)
                PrintConformation(s, seq2conf[s], ligand_tup = ligand_stuff)
        else:
            for s in sequences:
                PrintConformation(s, seq2conf[s])