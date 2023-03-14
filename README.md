# ALPACA

[![Supported Platforms](https://img.shields.io/badge/platforms-linux%20|%20osx-blue)](https://warpx.readthedocs.io/en/latest/install/users.html)
[![Language: C++17](https://img.shields.io/badge/language-C%2B%2B17-orange.svg)](https://isocpp.org/)
[![Language: Python](https://img.shields.io/badge/language-Python-orange.svg)](https://python.org/)
[![License ALPACA](https://img.shields.io/badge/license-GPL--3-blue)](https://spdx.org/licenses/GPL-3.0-only.html)

## Overview

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

On clusters, the two are likely going to be available as module to load. Outside of such computing environment, we need to make sure that we have them available on our system.

<details>
  <summary>MPI Installation Instructions</summary>
  
  To install and setup MPI, we have the choice of using [OpenMPI](https://www.open-mpi.org), and [MPICH](https://www.mpich.org). This instruction here is for OpenMPI, but applies equally as much for MPICH. Creating the build directory:

  ```bash
  mkdir mpi-build && export MPI_BUILD_DIR=$(PWD)/mpi-build
  ```
  
  To then begin the installation of MPI, we first have to download the source:

  ```bash
  wget https://download.open-mpi.org/release/open-mpi/v4.1/openmpi-4.1.5.tar.gz
  tar -xzf openmpi-4.1.5.tar.gz && cd openmpi-4.1.5
  ```

  We then have to configure our installation, and compile the library:

  ```bash
  ./configure --prefix=$MPI_BUILD_DIR
  make -j && make install
  ```
  
  After which we are left to export the MPI directories:

  ```bash
  export PATH=$MPI_BUILD_DIR/bin:$PATH
  export LD_LIBRARY_PATH=$MPI_BUILD_DIR/lib:$LD_LIBRARY_PATH
  ```

  > If your cluster environment comes with its own MPI library, you should **always** prefer using the system MPI library over doing a source install.

</details>

<details>
  <summary>HDF5 Installation Instructions</summary>
  
  To install HDF5, we roughly follow the same outlines as the ones for the MPI installation. Creating the build directory:

  ```bash
  mkdir hdf5-build && export HDF5_BUILD_DIR=$(pwd)/hdf5-build
  mkdir hdf5-install && export HDF5_INSTALL_DIR=$(pwd)/hdf5-install
  ```

  To then begin the installation of [HDF5](https://www.hdfgroup.org/downloads/hdf5/source-code/), we have to get the source, and then unpack it:

  ```bash
  wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8/hdf5-1.8.23/src/hdf5-1.8.23.tar.gz
  tar -xzf hdf5-1.8.23.tar.gz && cd hdf5-1.8.23
  ```

  Set the compilers to be the MPI-compilers:

  ```bash
  export CXX=mpic++
  export CC=mpicc
  ```

  After which we have to configure our installation, and then compile the library:

  ```bash
  cmake -GNinja -B ../hdf5-build/ -S . \
      -DCMAKE_INSTALL_DIR=$(pwd)/../hdf5-install \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_C_COMPILER=$(pwd)/../mpi-build/bin/mpicc \
      -DCMAKE_CXX_COMPILER=$(pwd)/../mpi-build/bin/mpic++ \
      -DHDF5_ENABLE_PARALLEL=On \
      -DHDF5_BUILD_CPP_LIB=On \
      -DALLOW_UNSUPPORTED=On
  ```
  
  To then build and install from the build directory
  
  ```bash
  cd $HDF5_BUILD_DIR
  ninja && ninja install
  ```

</details>

Having MPI & HDF5, we can then install ALPACA with

```bash
cmake -GNinja -B ../alpaca-build/ -S . \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_C_COMPILER=mpicc \
    -DCMAKE_CXX_COMPILER=mpicxx \
    -DHDF5_DIR=$HDF5_INSTALL_DIR/cmake
```

to build, we then invoke CMake again

```bash
cmake --build ../alpaca-build/
```

> We highly recommend using ``ccache`` together with CMake. To do so, add the following flags to the configuration step of CMake:
>
> ```bash
> -DCMAKE_C_COMPILER_LAUNCHER=ccache
> -DCMAKE_CXX_COMPILER_LAUNCHER=ccache
> ```

### Testing

To validate the installation, we recommend running unit-tests after the completed installation. To do so

```bash
ninja Paco -j 4
```

after which we can run single-, as well as two-core tests to verify the correctness of the installation.

```bash
mpiexec -n 1 ./Paco [1rank]
mpiexec -n 2 ./Paco [2rank]
```

For further instructions, first steps, and API documentation, please consult the ReadTheDocs.

## Academic Usage

If you use *ALPACA* in an academic setting, please cite [our papers](./CITATION.bib).

## Acknowledgments

*ALPACA* has received support from multiple funding bodies over the course of its inception:

* This project has received funding from the European Research Council (ERC) under the European Union’s Horizon 2020 research and innovation program: ERC Advanced Grant No. 667483, Prof. Dr. Nikolaus A. Adams, "NANOSHOCK - Manufacturing Shock Interactions for Innovative Nanoscale Processes"
* This project has received computing time on the GCS Supercomputer SuperMUC at Leibniz Supercomputing Centre (www.lrz.de) from the Gauss Centre for Supercomputing e.V. (www.gauss-centre.eu).
* This project has received funding from German Research Foundation (DFG).
* This project has received funding from the Bavarian State Ministry of Science and the Arts through the Competence Network for Scientific High Performance Computing in Bavaria (KONWIHR).
