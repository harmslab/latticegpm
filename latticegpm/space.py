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

    def __init__(self, wildtype, mutant, conformations, target_conf=None, temperature=1.0):
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
        # Construct a mutational mapping dictionary
        mutations = binary_mutations_map(wildtype, mutant)

        # Produce a binary sequence space between two sequences.
        genotypes = enumerate_space(wildtype, mutant)

        # Initialize a phenotypes array.
        phenotypes = np.zeros(len(genotypes), dtype=float)

        # Set initial properties of sequence space.
        self.temperature = temperature
        self.target_conf = target_conf

        super(LatticeMap, self).__init__(wildtype, genotypes, phenotypes, mutations=mutations)
        self.conformations = conformations

        # Fold proteins and calculate stabilities
        self.fold_proteins()

    def fold_proteins(self):
        """ Fold all sequences in protein sequence space. """
        # Fold proteins and extract stability parameter and native conformations
        confs = np.zeros(self.n,dtype="<U" +str(self.length))

        # Determine lattice protein native fold conformation.
        for i in range(self.n):
            fold = self.conformations.FoldSequence(self.genotypes[i], self.temperature, target_conf=self.target_conf)
            
            # If the lattice protein does not have a stable native state, do not fold protein (i.e. stability = 0)
            if fold[1] is None:
                self._phenotypes[i] = 0
                confs[i] = "U" * (self.length-1)
            
            # Else store stabilities and conformations
            else:
                self._phenotypes[i] = fold[0]
                confs[i] = fold[1]

        # Find all unique conformations
        self.unique_confs = np.unique(confs)
        self.n_conformations = len(self.unique_confs)

        # Set all conformations in the space.
        self.confs = confs

    def redefine_partition(self, z_confs):
        """ Calculate stabilities with a new set of states in partition function.
            Note that this changes the phenotypes in place.

            __Arguments__:

            `z_confs` [list] : conformations to include in partition func
            for stability calculations.

            __Returns__:

            `phenotypes` [array] : Array of stabilities recalculated.
        """
        # Set the partition function conformations as a attribute of the space.
        self.z_confs = z_confs

        # Calculate the partition functions.
        partition = np.empty(self.n, dtype=float)
        energies = np.empty(self.n, dtype=float)

        for i in range(self.n):
            # Sum over all conformations in z_conf for genotype i
            z = 0

            conf_energies = []
            for conf in self.z_confs:
                # Calculate folding energies of configuration
                fe = fold_energy(self.genotypes[i],conf)
                
                # Add config to partition function
                z += np.exp(-fe/self.temperature)
                
                # Store this energie for latter
                conf_energies.append(fe)

            # Set partition functions
            partition[i] = z

            # Find energy minimum from allowed conformations, this becomes native state.
            energies[i] = min(conf_energies)

        # Calculate the stabilities
        stabilities = energies - self.temperature * np.log(partition - np.exp(-energies/self.temperature))
        
        # Quality control... any NaN stabilities get set to 0 stability
        self._phenotypes = np.nan_to_num(stabilities)

        return self.phenotypes

    def print_sequences(self, sequences):
        """ Print sequence conformation with/without ligand bound. """

        # Get the sequence to conformation mapping from `seqspace` machinery.
        seq2conf = self.get_map("genotypes", "confs")
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
