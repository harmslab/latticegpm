from latticeproteins import LatticeProteins

# use space enumeration
from gpmap.gpm import GenotypePhenotypeMap
from gpmap.utils import mutations_to_genotypes

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

    conformations : latticeproteins.conformations.Conformations object
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

    For more attributes, see GenotypePhenotypeMap in `gpmap` package.
    """
    def __init__(self, wildtype,
        mutations,
        conformations=None,
        target=None,
        temp=1.0,
        phenotype_type="stability",
        **kwargs):

        # Get list of genotypes
        genotypes = mutations_to_genotypes(wildtype, mutations)

        # Calculate lattice proteins.
        self.latticeproteins = LatticeProteins(
            genotypes,
            conformations=conformations,
            target=target
        )

        # Get phentoype of interest.
        self._phenotype_type = phenotype_type
        phenotypes = getattr(self.latticeproteins, phenotype_type)

        # Build genotype-phenotype map.
        super(LatticeGenotypePhenotypeMap, self).__init__(
            wildtype,
            genotypes,
            phenotypes,
            mutations=mutations
        )

    @property
    def phenotype_type(self):
        return self._phenotype_type

    @phenotype_type.setter
    def phenotype_type(self, phenotype_type):
        self._phenotype_type = phenotype_type
        self.data['phenotypes'] = getattr(self.latticeproteins, phenotype_type)

    def print_sequences(self, sequences):
        """ Print sequence conformation with/without ligand bound. """
        # Get the sequence to conformation mapping from `seqspace` machinery.
        seq2conf = self.map("genotypes", "confs")
        for s in sequences:
            PrintConformation(s, seq2conf[s])
