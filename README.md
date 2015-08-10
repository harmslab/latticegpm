# LatticeGPM

LatticeGPM is a python package for building complete binary sequence spaces from protein lattice models. Lattice models provide a nice way to simulate the relationships between sequence, conformation, and function. Sequence Space provides a metaphor for thinking about evolution of sequence and function through this space. All lattice model calculations are done using [Dr. Jesse Bloom's](https://github.com/jbloom) `latticeproteins` package [here](https://github.com/jbloom/latticeproteins).


## Dependencies

The main dependency for this repository is the `seqspace` module found [here](https://github.com/harmslab/seqspace).

To run the notebooks, IPython and Jupyter notebooks must be installed.

Also, for constructing notebooks, NetworkX must be properly installed.

## Examples/Tutorials.

Check out the tutorials [here](https://github.com/harmslab/seqspace/blob/master/examples/Introduction%20to%20Genotype-Phenotype%20Map%20Module.ipynb).

## Installation

Git must be installed to clone and contribute to this project

### Setting up for use.

- Clone this repo locally:
```
git clone https://github.com/<user-name>/latticegpm
```
- Navigate to this directory, and install this python package with 
```
python setup.py develop
```

### Setting up for development.

- Fork this repository on Github
- Clone that repository locally
```
    git clone https://github.com/<user-name>/latticegpm
```
- Navigate to this directory, and install (softly) this python package with 
```
cd latticegpm
python setup.py develop
```
- Add another remote link to the master version, call it `upstream`.
```
git remote add upstream https://github.com/harmslab/latticegpm
```
- Start a branch locally from local master
```
git checkout -B <branch-name>
```
- Make changes and commit to that branch.
```
git commit -a -m "<commit message>"
```
- Push to your fork on github (which you called `upstream`).
```
git push upstream <branch-name>
```
- Pull request the branch on Github into this master repository on Github.


