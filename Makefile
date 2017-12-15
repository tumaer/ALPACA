##########################################################################################
#                                                                                        #
# This file is part of ALPACA                                                            #
#                                                                                        #
##########################################################################################
#  \\\\                                                                                  #
#  l '>                                                                                  #
#  | |                                                                                   #
#  | |                                                                                   #
#  | alpaca~                                                                             #
#  ||    ||                                                                              #
#  ''    ''                                                                              #
#                                                                                        #
# ALPACA                                                                                 #
# Copyright (c) 2017 Nikolaus A. Adams and contributors (see AUTHORS list)               #
# All rights reserved.                                                                   #
#                                                                                        #
# Chair of Aerodynamics and Fluid Mechanics                                              #
# Technical University of Munich                                                         #
#                                                                                        #
# This code is developed by the 'Nanoshock group' at the Chair of Aerodynamics and       #
# Fluid Mechanics, Technical University of Munich.                                       #
#                                                                                        #
# This project has received funding from the European Reseach Council (ERC)              #
# under the European Union's Horizon 2020 research and innovation programme              #
# (grant agreement No 667483).                                                           #
#                                                                                        #
# ERC Advanced Grant No 667483, Prof. Dr. Nikolaus A. Adams:                             #
# "NANOSHOCK - Manufacturing Shock Interactions for Innovative Nanoscale Processes"      #
#                                                                                        #
##########################################################################################
#                                                                                        #
# Redistribution and use in source and binary forms, with or without                     #
# modification, are permitted provided that the following conditions are met:            #
#                                                                                        #
# 1. Redistributions of source code must retain the above copyright notice,              #
#    this list of conditions and the following disclaimer.                               #
#                                                                                        #
# 2. Redistributions in binary form must reproduce the above copyright notice            #
#    this list of conditions and the following disclaimer in the documentation           #
#    and/or other materials provided with the distribution.                              #
#                                                                                        #
# 3. Neither the name of the copyright holder nor the names of its                       #
#    contributors may be used to endorse or promote products derived from this           #
#    software without specific prior written permission.                                 #
#                                                                                        #
# 4. Any redistribution of substantial fractions of the code as a                        #
#    different project should preserve the word ALPACA in the name                       #
#    of the code                                                                         #
#                                                                                        #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"            #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE              #
# IMPLIED WARRANTIES OF  MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE            #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE              #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR                    #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF                   #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS               #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN                #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)                #
# ARISING IN ANY WAY OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE            #
# POSSIBILITY OF SUCH DAMAGE.                                                            #
#                                                                                        #
# Please note, several third-party tools are used within the ALPACA code under           #
# their own license agreement.                                                           #
#                                                                                        #
# 1. xdmf_writer        : Licensed by Technische Universitaet Muenchen                   #
#                         See 'COPYING_XDMF_WRITER' for more information.                #
#                                                                                        #
# 2. tiny_xml           : This software is provided 'as-is', without any express or      #
#                         implied warranty. In no event will the authors be held         #
#                         liable for any damages arising from the use of this software.  #
#                         See COPYING_TINY_XMLfor more information.                      #
#                                                                                        #
# 3. expression_toolkit : Free use of The C++ Mathematical Expression Toolkit Library is #
#                         permitted under the guidelines and in accordance with the most #
#                         current version of the Common Public License.                  #
#                         http://www.opensource.org/licenses/cpl1.0.php                  #
#                         See COPYING_EXPRESSION_TOOLKITfor more information.            #
#                                                                                        #
##########################################################################################
#                                                                                        #
# AUTHORS                                                                                #
#                                                                                        #
#   Prof. Dr. Nikolaus A. Adams                                                          #
#                                                                                        #
#   Dr. Stefan Adami                                                                     #
#   Vladimir Bogdanov                                                                    #
#   Nico Fleischmann                                                                     #
#   Nils Hoppe                                                                           #
#   Naeimeh Hosseini                                                                     #
#   Jakob Kaiser                                                                         #
#   Aleksandr Lunkov                                                                     #
#   Thomas Paula                                                                         #
#   Josef Winter                                                                         #
#                                                                                        #
##########################################################################################
#                                                                                        #
# CONTACT                                                                                #
#                                                                                        #
#   nanoshock@aer.mw.tum.de                                                              #
#                                                                                        #
##########################################################################################
#                                                                                        #
# Munich, December 15th 2017                                                             #
#                                                                                        #
##########################################################################################

# source directory
SOURCE_DIR := ./src
# define a recursive wildcard function, as Make does not provide one
rwildcard=$(wildcard $1$2) $(foreach d,$(wildcard $1*),$(call rwildcard,$d/,$2))
CCFILES := $(call rwildcard,$(SOURCE_DIR)/,*.cpp)

# Compiler Choosing: GNU is assumed as standard unless working on the LRZ Cluster systems.
ifeq ($(MACHINE), LRZ) # SuperMUC
        CXX := mpiCC
	HDF5_LRZ_FLAGS = $(SZIP_LIB) -lz
	PERF_FLAGS = -ipo -fp-model precise
else
        CXX := mpic++
endif

# Compiler Flags for (64 bit systems only)
CCFLAGS = -std=c++11 -m64

# Setting up HDF5 and MPI includes if not already done by the operating system

ifndef HDF5_BASE
	HDF5_BASE = /group/nanoshock/opt/HDF5
endif

ifndef HDF5_INC
	HDF5_INC = -I$(HDF5_BASE)/include
endif

ifndef HDF5_LIB
	HDF5_LIB = -L$(HDF5_BASE)/lib -lhdf5_hl -lhdf5
endif

ifndef MPI_BASE
	MPI_BASE = /global/mpich-3.1
endif

ifndef MPI_INC
	MPI_INC = -I$(MPI_BASE)/include
endif

ifndef MPI_LIB
	MPI_LIB = -L$(MPI_BASE)/lib
endif

# Enable Debugging
ifeq ($(DBG),1)
        CCFLAGS   += -g -O0
else
        CCFLAGS   += -O3
endif

#-----------------------------------------------------------------------------
# Warning flags
#-----------------------------------------------------------------------------
WARN_FLAGS := -Wall -pedantic -W -Wformat -Wparentheses -Wmultichar -Wtrigraphs -Wpointer-arith -Wreturn-type -Wno-unused-function

#-----------------------------------------------------------------------------
# Programm Settings
#-----------------------------------------------------------------------------

# Default Values:
DIM = 3
# Internal Cells:
IC = 8

GITHASH := $(shell git rev-parse HEAD)
ifeq ($(GITHASH),)
	GITHASH := "NOT A GIT REPOSITORY"
endif

# Executable
EXECUTABLE      := ALPACA_$(DIM)D_IC$(IC)

#-----------------------------------------------------------------------------
# PREPARE AND COMPILE:
#-----------------------------------------------------------------------------

# Preprocessor Flags
DFLAGS = -D DIMENSION=$(DIM) -D GITHASH=$(GITHASH) -D INTERNALCELLS=$(IC)

#Combining the Flags

CXXFLAGS = $(CCFLAGS) $(PERF_FLAGS) $(DFLAGS) $(WARN_FLAGS) -fopenmp

# Include and Linker Setup

INCLUDES      := -I$(SOURCE_DIR) $(MPI_INC) $(HDF5_INC)

LDFLAGS       := $(MPI_LIB) $(HDF5_LIB) $(HDF5_LRZ_FLAGS)

#Setting up Object files away from the source (avaoids mistakenly adding them to git :D)
OBJECTDIR = ./objects

# Set up object files
OBJECTS 	:=  $(patsubst $(SOURCE_DIR)/%.cpp,$(OBJECTDIR)/%.o, $(CCFILES))

#Creating the Executable

$(EXECUTABLE): $(OBJECTS) Makefile
	$(CXX) $(CXXFLAGS) -o $(EXECUTABLE) $(OBJECTS) $(LDFLAGS)

$(OBJECTDIR)/%.o : $(SOURCE_DIR)/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ -c $<
.PHONY: clean softclean doc
clean:
	rm -rf $(EXECUTABLE) $(OBJECTDIR)
softclean:
	find $(OBJECTDIR) -path $(OBJECTDIR)/3rdParty -prune -o -name *.o ! -name user_expression.o  -exec rm -f {} +
	rm -f $(EXECUTABLE)
doc:
	doxygen Doxyfile
