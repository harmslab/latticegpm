import numpy as np
from latticeproteins.conformations import Conformations, BindLigand, PrintConformation
from latticegpm.utils import enumerate_space, fold_energy, ConformationError
from latticegpm.mapping import LatticeMap, LatticeFitnessMap

# use space enumeration
from seqspace.utils import binary_mutations_map
# ------------------------------------------------------
# Build a binary protein lattice model sequence space
# with fitness defined by function in Jesse Blooms'
# latticeproteins python package.
# ------------------------------------------------------

class LatticeConformationSpace(LatticeMap):
    
    def __init__(self, wildtype, mutant, conformations, target_conf=None, temperature=1.0, n_conformations=2):
        """ Build a protein lattice model sequence space from a conformation space. 
        
            Parameters:
            ----------
            wildtype: str
                Wildtype sequence
            mutant: str
                Mutant sequence
            conformations: Conformations object
                latticeproteins.conformations object for all conformations for
                strings with len(wildtype)
            target_conf: str (optional)
                String that describes the target conformation to fold each sequence to.
            temperature: float
                temperature parameter for calculating folding stability.
            n_conformations: int
                number of conformations that should be in space (raise error if more)
        
        """
        mutations = binary_mutations_map(wildtype, mutant)
        genotypes = enumerate_space(wildtype, mutant)
        phenotypes = np.zeros(len(genotypes),dtype=float)
        length = len(genotypes[0]) # string length
        
        if target_conf != None:
            n_conformations = 1
        
        # Fold proteins and extract stability parameter and native conformations
        confs = np.zeros(len(genotypes),dtype="<U" +str(length))
        
        # Determine lattice protein native fold conformation.
        for i in range(len(genotypes)):
            fold = conformations.FoldSequence(genotypes[i], temperature, target_conf=target_conf)
            # If the lattice protein does not have a stable native state, do not fold protein (i.e. stability = 0)
            if fold[1] is None:
                phenotypes[i] = 0
                confs[i] = "U" * (length-1)
            # Else -- store stabilities and conformations
            else:
                phenotypes[i] = fold[0]
                confs[i] = fold[1]
        
        # Find all unique conformations
        unique_confs = np.unique(confs)
        
        # Calculate the energies of all folds
        energies = np.zeros(len(genotypes), dtype=float)
        delta_e = np.zeros(len(genotypes), dtype=float)
        for i in range(len(genotypes)):            
            try:
                energies[i] = fold_energy(genotypes[i], confs[i])
                # If the space has more than n_conformations, raise error
                if len(unique_confs) == n_conformations:
                    delta_e[i] = energies[0] - energies[i]
                else:
                    raise Exception("More than " + str(n_conformations) + " state system.")
            except ConformationError:
                # Catch the unfolded proteins
                energies[i] = 0
                delta_e[i] = energies[0] - 0
                
        super(LatticeMap, self).__init__(wildtype, genotypes, phenotypes, mutations=mutations)        
        self.conformations = confs
        self.temperature = temperature
        self.target_conf = target_conf
        self.energies = energies
        self.delta_e = delta_e
        
        # Find all unique conformations.
        self.n_conformations = np.unique(confs)
                
    def fit_global_dG(self):
        """ Calculate the delta delta Gs (i.e. global stability) for the lattice model genotype-phenotype map."""
        pass
        
    def print_sequences(self, sequences):
        """ Print sequence conformation with/without ligand bound. """
        seq2conf = self.get_map("genotypes", "conformations")
        for s in sequences:
            PrintConformation(s, seq2conf[s])   
        

class LatticeFitnessSpace(LatticeFitnessMap):
    
    def __init__(self, wildtype, mutant, Fitness):
        """ Build a protein lattice model sequence space from a given fitness function. 
        
            Parameters:
            ----------
            wildtype: str
                Wildtype sequence
            mutant: str
                Mutant sequence
            Fitness: Fitness object
                latticeproteins.fitness object for protein that binds a ligand. 
        """
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