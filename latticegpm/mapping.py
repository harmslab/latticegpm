import numpy as np
from collections import OrderedDict
from latticeproteins.conformations import PrintConformation
from latticeproteins.sequences import HammingDistance

class LatticeMap(object):
    
    # ----------------------------------------
    # Get Properties of the Lattice Map
    # ----------------------------------------
    
    @property
    def length(self):
        """ Get length of the sequences. """
        return self._length
    
    @property
    def n(self):
        """ Get size of the sequence space. """
        return self._n
    
    @property
    def wildtype(self):
        """ Get reference sequences for interactions. """
        return self._wildtype
    
    @property
    def mutations(self):
        """ Get possible that occur from reference system. """
        return self._mutations
    
    @property
    def temperature(self):
        """ Get temperature of the system. """
        return self._temperature
    
    @property
    def sequences(self):
        """ Get sequences. """
        return self._sequences
        
    @property
    def stabilities(self):
        """ Get the dGs of the each sequence. """
        return self._stabilities
    
    @property
    def conformations(self):
        """ Get the native conformations of the each sequence. """
        return self._conformations
         
    # -------------------------------------------
    # Helpful maps that are built on the fly
    #--------------------------------------------
    @property
    def seq2fitness(self):
        """ Return dict of sequences mapped to fitnesses. """
        return self._map(self.sequences, self.fitnesses)
    
    @property
    def seq2index(self):
        """ Return an ordered dictionary mapping sequences to indices. """
        return self._map(self.sequences, self._indices)
    
    @property
    def seq2conformation(self):
        """ Return an ordered dictionary of sequences mapped to their native conformations. """
        return self._map(self.sequences, self.conformations)
        
    # -------------------------------------------
    # Setting methods for properties
    # -------------------------------------------
    
    @sequences.setter
    def sequences(self, sequences):
        """ Set sequences from ordered list of sequences. """
        self._n = len(sequences)
        self._length = len(sequences[0])
        self._sequences = np.array(sequences)
        self._indices = np.arange(self._n)
        
    @wildtype.setter
    def wildtype(self, wildtype):
        """ Set the reference sequence among the mutants in the system. """
        self._wildtype = wildtype
        self._mutations = self._differ_all_sites(wildtype)
        
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
            if len(conformations) != len(self._sequences):
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
            if len(stabilities) != len(self._sequences):
                raise("Number of stabilities does not equal number of sequences.")
            else:
                self._stabilities = stabilities
                

    # --------------------------------------------
    # Useful methods for mapping object
    # --------------------------------------------
        
    def _map(self, keys, values):
        """ Return ordered dictionary mapping two properties in self. """
        return OrderedDict([(keys[i], values[i]) for i in range(self.n)])

    def _differ_all_sites(self, reference):
        """ Find the sequence in the system that differs at all sites from reference.""" 
        for sequence in self.sequences:
            if HammingDistance(sequence, reference) == self.length:
                differs = sequence
                break
        return sequence
        
        
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