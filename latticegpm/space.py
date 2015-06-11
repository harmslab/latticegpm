import numpy as np
from latticeproteins.conformations import Conformations, BindLigand, PrintConformation
from latticegpm.utils import enumerate_space, fold_energy, ConformationError
from latticegpm.mapping import LatticeMap, LatticeFitnessMap

# ------------------------------------------------------
# Build a binary protein lattice model sequence space
# with fitness defined by function in Jesse Blooms'
# latticeproteins python package.
# ------------------------------------------------------

class LatticeConformationSpace(LatticeMap):
    
    def __init__(self, wildtype, mutant, conformations, target_conf=None, temperature=1.0):
        """ Build a protein lattice model sequence space from a conformation space. """
        self.sequences = enumerate_space(wildtype, mutant)
        self.wildtype = wildtype
        self.temperature = temperature
        self.target_conf = target_conf
        
        # Fold proteins and extract stability parameter and native conformations
        self.stabilities = np.zeros(len(self.sequences),dtype=float)
        self.conformations = np.zeros(len(self.sequences),dtype="<U" +str(self.length))
        
        for i in range(len(self.sequences)):
            fold = conformations.FoldSequence(self.sequences[i], self.temperature, target_conf=self.target_conf)
            if fold[1] is None:
                self.stabilities[i] = 0
                self.conformations[i] = "U" * (self.length-1)
            else:
                self.stabilities[i] = fold[0]
                self.conformations[i] = fold[1]
        unique = np.unique(self.conformations)
        print(unique)
        # Calculate the energies of all folds
        self.energies = np.zeros(len(self.sequences), dtype=float)
        self.delta_e = np.zeros(len(self.sequences), dtype=float)
        for i in range(len(self.sequences)):            
            try:
                self.energies[i] = fold_energy(self.sequences[i], self.conformations[i])
                if len(unique) <= 2:
                    self.delta_e[i] = self.energies[0] - self.energies[i]
                else:
                    raise Exception("""More than two state system.""")
            except ConformationError:
                self.energies[i] = 0
                self.delta_e[i] = self.energies[0] - 0
                
                
    def fit_ddG(self):
        """ Calculate the delta delta Gs for the lattice model genotype-phenotype map."""
        
        
        
        
    def print_sequences(self, sequences):
        """ Print sequence conformation with/without ligand bound. """
        seq2conf = self.seq2conformation
        for s in sequences:
            PrintConformation(s, seq2conf[s])   
        

class LatticeFitnessSpace(LatticeFitnessMap):
    
    def __init__(self, wildtype, mutant, Fitness):
        """ Build a protein lattice model sequence space from a given fitness function. """
        self.sequences = enumerate_space(wildtype, mutant)
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