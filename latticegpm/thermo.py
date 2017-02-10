#
# This module contains random functions that I found useful for probing
# lattice protein sequence space. Note: these were not designed with
# speed/efficiency in mind. They are bit crude in their implementation.
#
import itertools as it
import numpy as np

# ------------------------------------------------------------
# Jesse Bloom's Lattice Model imports
# ------------------------------------------------------------

from latticeproteins.interactions import miyazawa_jernigan

class LatticeThermodynamics(object):
    """Calculate Lattice thermodynamics for a sequence from a list of conformations.

    Currently, doesn't do a lot of quality control
    """
    def __init__(self, sequence, conf_list, temperature, interaction_energies=miyazawa_jernigan):
        self.sequence = sequence
        self.conf_list = conf_list
        self.temperature = temperature
        self.interaction_energies = interaction_energies

    @property
    def energies(self):
        """Get the energies of all conformations
        """
        try:
            return self._energies
        except AttributeError:
            self._energies = energy_list(self.sequence,
                self.conf_list,
                interaction_energies=self.interaction_energies)
            conf_i = np.where(self._energies == self._energies.min())[0]
            self.native_conf = self.conf_list[conf_i][0]
            return self._energies

    @property
    def partition_function(self):
        """Get the partition sum for lattice."""
        try:
            return self._partition_sum
        except AttributeError:
            self._partition_sum = partition_function_from_energies(
                self.energies,
                self.temperature)
            return self._partition_sum

    @property
    def stability(self):
        """Get stability of the lattice protein
        """
        try:
            return self._stability
        except AttributeError:
            self._stability, self._folded = stability_from_energies(
                self.energies,
                self.temperature)
            return self._stability

    @property
    def folded(self):
        """Get folded attribute."""
        try:
            return self._folded
        except AttributeError:
            self._stability, self._folded = stability_from_energies(
                self.energies,
                self.temperature)
            return self._folded

    @property
    def fracfolded(self):
        """Get fraction folded."""
        try:
            return self._fracfolded
        except AttributeError:
            self._fracfolded = fracfolded_from_stability(self.stability,
                self.temperature)
            return self._fracfolded

def energy_list(sequence, conf_list, interaction_energies=miyazawa_jernigan):
    """Calculate a energies from a list of conformations for a given sequence.
    """
    return np.array([fold_energy(sequence, c) for c in conf_list])

def partition_function_from_energies(energies, temperature):
    """Calculate a partition function from a list of energies.
    """
    energies = np.array(energies)
    boltzmann = np.exp(-energies / temperature)
    return sum(boltzmann)

def partition_function(sequence, conf_list, temperature, interaction_energies=miyazawa_jernigan):
    """Calculate a partition sum from a list of conformations.
    """
    energies = energy_list(sequence, conf_list, interaction_energies=interaction_energies)
    return partition_function_from_energies(energies, temperature)

def stability_from_conf_list(sequence, conf_list, temperature, interaction_energies=miyazawa_jernigan):
    """Calculate stabilities from list of conformations.

    Returns
    -------
    stability : float
        Stability of sequence given conf_list
    folded : bool
        True if the protein folded, False if not.
    """
    # energies
    energies = energy_list(sequence, conf_list, interaction_energies=interaction_energies)
    # partition function
    partition = partition_function_from_energies(energies, temperature)
    # native energy
    minE = energies[energies==energies.min()]
    if len(minE) > 1:
        return 0, False
    # Calculate stabilities
    return minE[0] + temperature * np.log(partition - np.exp(-minE[0] / temperature)), True

def stability_from_energies(energies, temperature):
    """Calculate stability from list of energies.

    Returns
    -------
    stability : float
        Stability of sequence given conf_list
    folded : bool
        True if the protein folded, False if not.
    """
    # partition function
    partition = partition_function_from_energies(energies, temperature)
    # native energy
    minE = energies[energies==energies.min()]
    if len(minE) > 1:
        return 0, False
    # Calculate stabilities
    return minE[0] + temperature * np.log(partition - np.exp(-minE[0] / temperature)), True

def fracfolded_from_conf_list(sequence, conf_list, temperature, interaction_energies=miyazawa_jernigan):
    """Calculate staiblity from a list of conformations
    """
    stability, folded = stability_from_conf_list(sequence,
        conf_list,
        temperature,
        interaction_energies=interaction_energies)
    return 1.0 / (1 + np.exp(stability/temperature))

def fracfolded_from_energies(energies, temperature):
    """Calculate a fraction folded from a list of energies.
    """
    stability, folded = stability_from_energies(energies, temperature)
    return 1.0 / (1 + np.exp(stability/temperature))

def fracfolded_from_stability(stability, temperature):
    """Calculate a fraction folded from stability
    """
    return 1.0 / (1 + np.exp(stability/temperature))

def fold_energy(sequence, conformation, interactions=miyazawa_jernigan):
    """Calculate the energy of the sequence with the given conformation.

    Parameters
    ----------
    sequence : str
        Amino acid sequence to fold.
    conformation : str
        Conformation according to latticemodel's conformations format (e.g. 'UDLLDRU')

    Returns
    ------
    energy : float
        energy of the conformation (sum of all contact energies)
    """
    contacts = lattice_contacts(sequence, conformation)
    energy = sum([interactions[c] for c in contacts])
    return energy

def lattice_contacts(sequence, conformation):
    """ Find all contacts in conformation.

    Parameters
    ----------
    sequence : str
        Amino acid sequence to fold.
    conformation : str
        Conformation according to latticemodel's conformations format (e.g. 'UDLLDRU')

    Returns
    -------
    contacts : list
        list of contact pairs
    """
    sites = list(sequence)
    length = len(sites)
    try:
        moves = list(conformation)
    except TypeError:
        raise ConformationError("""Protein conformation is None; is there a native state? """)
    # build a coordinate system, note that odd rotation of intuitive coordinates
    # since we are working in numpy array grid.
    coordinates = {"U": [-1,0], "D":[1,0], "L":[0,-1], "R":[0,1]}
    grid = np.zeros((2*length+1, 2*length+1), dtype=str)
    x = y = round(length/2.0) # initial position on the grid is at the center of the 2d array
    grid[x,y] = sites[0]
    # move on grid, populate with amino acid at that site, and store all contacting neighbors.
    contacts = []
    for i in range(length-1):
        step = coordinates[moves[i]]
        x += step[0]
        y += step[1]
        grid[x,y] = sites[i+1]
        neighbors = [sites[i+1] + grid[x+c[0], y+c[1]] for c in coordinates.values()]
        contacts += [n for n in neighbors if n in miyazawa_jernigan]
    # subtract the contacts that have bonds between them.
    for i in range(1,length):
        try:
            contacts.remove(sequence[i-1:i+1])
        except ValueError:
            contacts.remove(sequence[i] + sequence[i-1])
    return contacts
