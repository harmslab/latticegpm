import os
import numpy as np
from latticeproteins.conformations import Conformations, BindLigand, PrintConformation
from latticeproteins.interactions import miyazawa_jernigan
from latticeproteins.sequences import RandomSequence
from latticegpm.utils import fold_energy, ConformationError
from . import search

# use space enumeration
from seqspace.gpm import GenotypePhenotypeMap
from seqspace.utils import binary_mutations_map, enumerate_space

# ------------------------------------------------------
# Build a binary protein lattice model sequence space
# with fitness defined by function in Jesse Blooms'
# latticeproteins python package.
# ------------------------------------------------------

class LatticeGenotypePhenotypeMap(GenotypePhenotypeMap):
    """Construct a genotype-phenotype map from 2d protein lattice models. Genotypes
    are represented by

    Parameters
    ----------
    wildtype : str
        Wildtype sequence
    mutant : str
        Mutant sequence
    target_conf : str (optional)
        String that describes the target conformation to fold each sequence to.
    temperature : float
        temperature parameter for calculating folding stability.
    n_conformations : int
        number of conformations that should be in space (raise error if more)
    conformations : Conformations object
        latticeproteins.conformations object for all conformations for
        strings with len(wildtype)

    Attributes
    ----------
    temperature : float
        temperature of the system
    target_conf : string
        target conformation for all sequences
    interaction_energies : dictionary
        mapping all contact pairs to their corresponding energies
    energies : array of floats
        native energies of all sequences
    stabilities : array of floats
        stabilities of native states for all sequences
    fitnesses : array of floats
        fitnesses of all sequences
    conformations : array of strings
        native conformations of all sequences
    partition_sum : array of float
        boltzmann weighted sum of all conformation/energies
    fold : boolean
        folded or not?

    For more attributes, see GenotypePhenotypeMap in `seqspace` package.
    """
    def __init__(self, wildtype,
            mutations,
            target_conf=None,
            temperature=1.0,
            interaction_energies=miyazawa_jernigan,
            database_dir="database/",
        ):
        # Construct a mutational mapping dictionary
        self.wildtype = wildtype
        self.mutations = mutations
        # Set initial properties of sequence space.
        self.temperature = temperature
        self.target_conf = target_conf
        self.interaction_energies = interaction_energies
        # Initialize a phenotypes array.
        phenotypes = np.empty(len(genotypes), dtype=float)
        # Construct a genotype-phenotype map
        super(LatticeMap, self).__init__(wildtype, genotypes, phenotypes, mutations=mutations)
        # Make a conformations database.
        if not os.path.exist(database_dir):
            os.makedirs(database_dir)
        # Construct conformations database.
        self.Conformations = Conformations(self.length, database_dir, interaction_energies=interaction_energies)
        # Set the fitness.
        self.Fitness = Fitness(self.temperature, self.Conformations,
            dGdependence=None,
            targets=target_conf)
        # Fold proteins and calculate stabilities
        self.build()

    @property
    def phenotypes(self):
        """Get the phenotypes, specified by phenotype_type attribute.
        """
        return getattr(self, self.phenotype_type)

    def build(self):
        """Build the LatticeGenotypePhenotypeMap from `latticeprotein` package.

        Uses the `latticeprotein` package to calculate the following attributes
        and add them to the genotype-phenotype map.

        Attributes
        ----------
        energies : array of floats
            native energies of all sequences
        stabilities : array of floats
            stabilities of native states for all sequences
        fitnesses : array of floats
            fitnesses of all sequences
        conformations : array of strings
            native conformations of all sequences
        partition_sum : array of float
            boltzmann weighted sum of all conformation/energies
        fold : boolean
            folded or not?
        """
        self.energies = np.empty(self.n, dtype=float)
        self.stabilities = np.empty(self.n, dtype=float)
        self.fitnesses = np.empty(self.n, dtype=float)
        self.conformations = np.empty(self, dtype=str)
        self.partition_sum = np.empty(self.n, dtype=float)
        self.fold = np.empty(self.n, dtype=bool)
        for i, g in enumerate(self.genotypes):
            results = self.Fitness._AllMetrics(g)
            self.energies[i] = results[0]
            self.stabilites = results[1]
            self.fitnesses[i] = results[2]
            self.conformations[i] = results[3]
            self.partition_sum = results[4]
            self.fold[i] = results[-1]

    @classmethod
    def with_length(cls, length, temperature=1.0, differby=None, **kwargs):
        """Searches regions of sequences space for a lattice proteins
        with the given length on calculates their fitness.
        """
        seq1, seq2 = search.sequence_space(length, temperature,
            differby=differby,
            **kwargs)

        return cls(seq1, )
        #return cls()

    @classmethod
    def with_genotypes(cls, sequences, **kwargs):
        """Calculates a set of lattice proteins from sequences.
        """
        pass
        #return cls()

    def set_phenotypes(self, attr="stabilities"):
        """Set the phenotype to a different attribute.
        """
        self.phenotype_type = attr

    def recalculate_partition_sum(self, zconfs, target_conf=None):
        """Recalculate stabilities for all sequences with new manually defined
        states in partition function.

        Note that the partition conformations always include the E=0 state (unfolded).

        Parameters
        ----------
        z_confs : lists
            conformations to include in partition func for stability calculations.

        Returns
        -------
        phenotypes : array
            Array of stabilities recalculated.
        """
        # Set the partition function conformations as a attribute of the space.
        self.z_confs = z_confs
        for i in range(self.n):
            # Sum over all conformations in z_conf for genotype i
            z = 0   # Start with the completely unfolded conformation
            # The completely unfolded state is a possible configuration
            conf_energies = [0]
            min_conf = ""
            min_energy = 0
            fold = True
            for conf in self.z_confs:
                # Calculate folding energies of configuration
                fe = fold_energy(self.genotypes[i], conf, self.interaction_energies)
                # Add config to partition function
                z += np.exp(-fe/self.temperature)
                # Store this energie for latter
                conf_energies.append(fe)
                # Set the minimum energy.
                if fe < min_energy:
                    min_energy = fe
                    min_conf = conf
                    fold=True
                # If it's not a single lowest state, its nonnative.
                elif fe == min_energy:
                    fold=False
            # Set partition functions
            self.fold[i] = fold
            self.partition_sum[i] = z
            # If target conformation is given, use it as native state.
            if target_conf is not None:
                self.energies[i] = fold_energy(self.genotypes[i], target_conf, self.interaction_energies)
            # set the non-native
            elif fold:
                # Find energy minimum from allowed conformations, this becomes native state.
                self.energies[i] = min_energy
                self.conformations[i] = min_conf
            else:
                self.conformations[i] = "U" * (self.length-1)
                self.energies[i] = 0
            # Prepare options for calcs below.
            info = [self.conformations[i], self.partition_sum[i], 0, self.fold[i]]
            # Calculate stability
            self.stabilities[i] = self.Fitness._Stabilty(self.energies[i], *info)
            # Calculate the fold
            self.fitnesses[i] = self.Fitness._Fitness(self.energies[i], *info)

    def print_sequences(self, sequences):
        """ Print sequence conformation with/without ligand bound. """

        # Get the sequence to conformation mapping from `seqspace` machinery.
        seq2conf = self.get_map("genotypes", "confs")
        for s in sequences:
            PrintConformation(s, seq2conf[s])
