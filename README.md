# ALPACA

## About

*ALPACA* is an MPI-parallelized C++ code framework to simulate compressible multiphase flow physics. It allows for advanced high-resolution sharp-interface modeling empowered with efficient multiresolution compression. The modular code structure offers a broad flexibility to select among many most-recent numerical methods covering WENO/T-ENO, Riemann solvers (complete/incomplete), strong-stability preserving Runge-Kutta time integration schemes, level-set methods and many more.

## Getting Started

### Installation

Recursively cloning in the case of a fresh installation

```bash
git clone --recursive https://github.com/tumaer/ALPACA.git
```

or in the case of an existing download

```bash
git fetch && git submodule update --init --recursive
```

After which we first need to install ALPACA's dependencies, ALPACA depends on

* MPI
* HDF5

> Further build instructions to be added.

### Testing

To validate the installation, we recommend running unit-tests after the completed installation. To do so

```bash
make Paco -j 4
```

after which we can run single-, as well as two-core tests to verify the correctness of the installation.

```bash
mpiexec -n 1 ./Paco [1rank]
mpiexec -n 2 ./Paco [2rank]
```

For further instructions, first steps, and API documentation, please consult the ReadTheDocs.

## Academic Usage

If you use *ALPACA* in an academic setting, please cite [our papers](./CITATION.cff).

## Acknowledgments

*ALPACA* has received support from multiple funding bodies over the course of its inception:

* This project has received funding from the European Research Council (ERC) under the European Union’s Horizon 2020 research and innovation program: ERC Advanced Grant No. 667483, Prof. Dr. Nikolaus A. Adams, "NANOSHOCK - Manufacturing Shock Interactions for Innovative Nanoscale Processes"
* This project has received computing time on the GCS Supercomputer SuperMUC at Leibniz Supercomputing Centre (www.lrz.de) from the Gauss Centre for Supercomputing e.V. (www.gauss-centre.eu).
* This project has received funding from German Research Foundation (DFG).
* This project has received funding from the Bavarian State Ministry of Science and the Arts through the Competence Network for Scientific High Performance Computing in Bavaria (KONWIHR).
