import numpy as np
from collections import OrderedDict
from latticeproteins.conformations import PrintConformation
from latticeproteins.sequences import HammingDistance
from seqspace.gpm import GenoPhenoMap

class LatticeMap(GenoPhenoMap):
    
    def __init__(self, wildtype, genotypes, phenotypes, mutations=None):
        """ """
        super(LatticeMap, self).__init__(self, wildtype, sgenotypes, phenotypes, mutations=mutations)
    
    # ----------------------------------------
    # Get Properties of the Lattice Map
    # ----------------------------------------
    
    @property
    def temperature(self):
        """ Get temperature of the system. """
        return self._temperature
        
    @property
    def stabilities(self):
        """ Get the dGs of the each sequence. """
        return self._stabilities
    
    @property
    def conformations(self):
        """ Get the native conformations of the each sequence. """
        return self._conformations
        
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
        if type(conformations) is dict:
            self._conformations = self._if_dict(conformations)
        else:
            if len(conformations) != len(self.genotypes):
                raise("Number of conformations does not equal number of sequences.")
            else:
                self._conformations = conformations
                
    @stabilities.setter
    def stabilities(self, stabilities):
        """ Set dGs of all native conformations of sequences in map. 
        """
        if type(stabilities) is dict:
            self._stabilities = self._if_dict(stabilities)
        else:
            if len(stabilities) != len(self.genotypes):
                raise("Number of stabilities does not equal number of sequences.")
            else:
                self._stabilities = stabilities
        
class LatticeFitnessMap(LatticeMap):
    
    # ----------------------------------------
    # Get Properties of the Lattice Map
    # ----------------------------------------
    
    @property
    def fitnesses(self):
        """ Get fitness. """
        return self._fitnesses
        
    @property
    def ligand_tup(self):
        """ Get ligand tuple (ligand sequence, ligand conformation, stability cutoff). """
        return self.ligand_tup
    
    @property
    def ligand(self):
        """ Get ligand sequence"""
        return self._ligand_tup[0]
        
    @property
    def ligand_conformation(self):
        """ Get ligand conformation. """
        return self._ligand_tup[1]
        
    # -------------------------------------------
    # Setting methods for properties
    # -------------------------------------------
        
    @ligand.setter
    def ligand_tup(self, ligand):
        """ Set the ligand for binding fitnesses. """
        if type(ligand) != tuple:
            raise TypeError("ligand must be type==tuple.")
        self._ligand_tup = ligand
        
    @fitnesses.setter
    def fitnesses(self, fitnesses):
        """ NORMALIZE and set fitnesses from ordered list of fitnesses 
            
            Args:
            -----
            fitnesses: array-like or dict
                if array-like, it musted be ordered by sequences; if dict,
                this method automatically orders the fitnesses into numpy
                array.
        """
        if type(fitnesses) is dict:
            self._fitnesses = self._if_dict(fitnesses)/fitnesses[self.wildtype]
        else:
            if len(fitnesses) != len(self._sequences):
                raise("Number of phenotypes does not equal number of sequences.")
            else:
                wildtype_index = self.seq2index[self.wildtype]
                self._fitnesses = fitnesses/fitnesses[wildtype_index]