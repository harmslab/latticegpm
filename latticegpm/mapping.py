import numpy as np
from collections import OrderedDict
from latticeproteins.conformations import PrintConformation
from latticeproteins.sequences import HammingDistance
from seqspace.gpm import GenotypePhenotypeMap

class LatticeBaseMap(GenotypePhenotypeMap):
    """
    """
    def __init__(self, wildtype, genotypes, phenotypes, mutations=None):
        super(LatticeMap, self).__init__(self, wildtype, genotypes, phenotypes, mutations=mutations)

    # ----------------------------------------
    # Get Properties of the Lattice Map
    # ----------------------------------------

    @property
    def temperature(self):
        """Get temperature of the system. """
        return self._temperature

    @property
    def stabilities(self):
        """Get the dGs of the each sequence. """
        return self._stabilities

    @property
    def conformations(self):
        """Get the native conformations of the each sequence. """
        return self._conformations

    @property
    def energies(self):
        """Get the native energies of all genotypes"""
        self._energies

    @property
    def fitnesses(self):
        """Get the fitness of genotypes.
        """
        self._fitnesses

    @property
    def unique_conformations(self):
        """return an array of conformations that are unique.
        """
        return np.unique(self._conformations)

    @property
    def n_conformations(self):
        """Return the number of unique conformations
        """
        # Find all unique conformations
        return len(self.unique_conformations)

    # -------------------------------------------
    # Setting methods for properties
    # -------------------------------------------

    @temperature.setter
    def temperature(self, temperature):
        """ Set the temperature of the system. """
        self._temperature = temperature

    @conformations.setter
    def conformations(self, conformations):
        """ Set native conformations of all sequences in space.
        """
        self._conformations = conformations

    @stabilities.setter
    def stabilities(self, stabilities):
        """ Set dGs of all native conformations of sequences in map.
        """
        self._stabilities = stabilities
