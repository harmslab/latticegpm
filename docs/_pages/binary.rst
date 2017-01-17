Binary lattice model map
========================

Binary lattice model maps provide straight-forward systems for studying evolutionary
dynamics, trajectories, and epistatic interactions between mutations. To construct
one of these maps, ``latticegpm`` includes a set of functions for selecting two
random sequences that fold below a given $\Delta G$.

.. code-block:: python

    # Search function in utils module
    from latticegpm.utils import search_conformation_space

    # Enumerate conformations on lattice
    length = 10
    temperature = 1.0
    stability_max = 0.0
    database = "path/to/database"
    c = conformations(length, database)

    # Get two sequences that differ at all sites
    sequence1, sequence2 = search_conformation_space(c, temperature, stability_max)


The full documentation for searching conformation space:

.. autofunction:: latticegpm.utils.search_conformation_space
