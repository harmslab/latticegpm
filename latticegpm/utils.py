

class ConformationError(Exception):
    """ If protein doesn't fold, raise error. """

def compare_sequences(s1, s2):
    """ Return the indice where two strings differ. """
    return [i for i in range(len(s1)) if s1[i] != s2[i]]

def mutations_map(s1, s2):
    """ Construct a mutations dictionary for latticegpm between
    to sequences, s1 and s2.

    Example
    -------
        {site-number: [alphabet]}
    """
    mutations = dict()
    for i in range(len(s1)):
        if s1[i] == s2[i]:
            mutations[i] = None
        else:
            mutations[i] = [s1[i], s2[i]]
    return mutations
