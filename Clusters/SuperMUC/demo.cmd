#!/bin/bash
##
#@ energy_policy_tag = NONE
##
#@ island_count = 1 
#@ job_type = MPICH 
#@ class = general
#@ node = 2
#@ tasks_per_node = 28
#@ wall_clock_limit = 24:00:00
#@ job_name = Name
#@ network.MPI = sn_all,not_shared,us
#@ initialdir = dir
#@ output = $(home)/Results/Name.Out.$(jobid)
#@ error =  $(home)/Results/Name.Err.$(jobid)
#@ notification=complete
#@ notify_user=your.name@provider.domain
#@ queue
. /etc/profile
. /etc/profile.d/modules.sh
#setup environment
module unload mpi.ibm
module unload poe
module unload intel
module load gcc/5
module load intel
module load mpi.intel
module load hdf5/mpi

mpiexec -n 56 ./ALPACA ./inputfile.xml
