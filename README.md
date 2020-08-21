# About

This repository contains the code I developed as part of my PhD thesis, <a href="https://orbit.dtu.dk/en/publications/topology-optimization-and-lattice-boltzmann-methods">"Topology optimization and lattice Boltzmann methods"</a>.

# Usage

The code is _not_ in a "plug and play" state to use as a library, as it was developed with myself as the sole user. Similarly, it contains bits of code that were never finished. However, bits and pieces of the code might be useful for others to incorporate into their own research projects.

## Potentially useful code for others

The following bits of code might be of interest to others:

* The <a href="https://github.com/northgaard/LBM-TopOpt/tree/master/src/topOpt">topOpt</a> folder contains a low overhead implementation of a memory efficient checkpointing algorithm for adjoint computations. See <a href="https://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.527.8486">original paper</a> by Wang _et al_ for theoretical details.
* <a href="https://github.com/northgaard/LBM-TopOpt/blob/master/src/core/MPISplit.cc">MPISplit.cc</a> contains the implementation of a helper class used to solve multiple systems in parallel in separate MPI communicators, useful for <a href="https://link.springer.com/article/10.1007/s00158-010-0602-y">robust topology optimization</a>.

These can both be used in another context than lattice Boltzmann, e.g. for topology optimization for other kinds of physics or for other kinds of optimization problems.

# Dependencies

The code depends on <a href="https://www.mcs.anl.gov/petsc/">PETSc</a>, and optionally on <a href="https://www.scicomp.uni-kl.de/software/codi/">CoDiPack</a> for differentiation of arbitrary lattice Boltzmann models using automatic differentiation.

# My papers

* <a href="https://www.sciencedirect.com/science/article/pii/S0021999115008426">Topology optimization for unsteady flow problems using the lattice Boltzmann method</a>
* <a href="https://link.springer.com/article/10.1007%2Fs00158-017-1708-2">Applications of automatic differentiation in topology optimization</a>

# Contact

In case any of the above is potentially useful to you, feel free to open an issue in case you need help.
