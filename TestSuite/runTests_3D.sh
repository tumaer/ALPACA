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

#!/bin/bash
if [ "$#" -lt 2 ] 
then
	echo "-----------------------------------------"
	echo "This is the 3D testsuite for ALPACA"
	echo "Please launch the testsuite with an argument for the outputfolder and the respective testcase:"
	echo "-----------------------------------------"
	echo ""
	echo "./runTests_3D.sh foldername all"
	echo "     -> Compile binaries and run all Tests"
	echo ""
	echo "./runTests_3D.sh foldername all nomake"
	echo "     -> Run all Tests with pre-existing binaries (for testing only recommended)" 
	echo ""
	echo "./runTests_3D.sh foldername 300 (nomake)"
	echo "     -> Sod shock-tube problem in X: Lmax=0, IC=8"
	echo "     -> 1 rank | 4 ranks | 8 ranks"
	echo "     -> reverse-sod: 1rank | 4 ranks | 8 ranks"
	echo ""
	echo "./runTests_3D.sh foldername 301 (nomake)"
	echo "     -> Sod shock-tube problem in Y: Lmax=0, IC=8"
	echo "     -> 1 rank | 4 ranks | 8 ranks"
	echo "     -> reverse-sod: 1rank | 4 ranks | 8 ranks"
	echo ""
	echo "./runTests_3D.sh foldername 302 (nomake)"
	echo "     -> Sod shock-tube problem in Z: Lmax=0, IC=8"
	echo "     -> 1 rank | 4 ranks | 8 ranks"
	echo "     -> reverse-sod: 1rank | 4 ranks | 8 ranks"
	echo ""

	exit 0
fi

DIRNAME=$1
TESTCASENUMBER=$2
if [ $# -lt 3 ]; then
	MAKEBINARIES="yes"
elif [ "$4" == "nomake" ]; then
	MAKEBINARIES="no"
fi

# Start.
echo "---------------------------------------------"
echo "-                                           -"
echo "-      Starting 3D-Testsuite for ALPACA     -"
echo "-                                           -"
echo "---------------------------------------------"
echo "ParaviewPath: " $PARAVIEW_BASE

# Copy binaries -> 3D
cd ..
if [ "$MAKEBINARIES" == "yes" ]; then
	make DIM=3 IC=8 clean
	make DIM=3 IC=8 -j
fi
cp ALPACA_3D_IC8 ./TestSuite/$DIRNAME/ALPACA_3D_IC8
if [ $? -ne 0 ]; then
	echo "Could not copy binary..."
	exit 1
fi

cd TestSuite

echo "3D Binaries copied."
echo "########################################################################"

numcases=0
failedcases=0

#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------
testnumber=300
if [ "$TESTCASENUMBER" == "$testnumber" ] || [ "$TESTCASENUMBER" == "all" ]; then
	echo "Test "$testnumber" started..."
	# format: dirname 1D/2D/3D mpi-ranks Lmax IC run/postproc
	./Cases/testcase_SOD_X.sh $DIRNAME 3 1 0 8 run
	if [ $? -eq 0 ]; then
		echo "Test "$testnumber" (sim) -> DONE"; 
	else	
		echo "Test "$testnumber" (sim) -> FAILED - FAILED - FAILED - FAILED - FAILED"; 
		failedcases=1
	fi
	./Cases/testcase_SOD_X.sh $DIRNAME 3 1 0 8 postproc
	if [ $? -eq 0 ]; then
		echo "Test "$testnumber" (postproc) -> DONE"; 
	else	
		echo "Test "$testnumber" (postproc) -> FAILED - FAILED - FAILED - FAILED - FAILED"; 
		failedcases=1
	fi
	echo "########################################################################"
	numcases=$(( $numcases + 1 ))

	# format: dirname 1D/2D/3D mpi-ranks Lmax IC run/postproc
	./Cases/testcase_SOD_X.sh $DIRNAME 3 4 0 8 run
	if [ $? -eq 0 ]; then
		echo "Test "$testnumber" (sim) -> DONE"; 
	else	
		echo "Test "$testnumber" (sim) -> FAILED - FAILED - FAILED - FAILED - FAILED"; 
		failedcases=1
	fi
	./Cases/testcase_SOD_X.sh $DIRNAME 3 4 0 8 postproc
	if [ $? -eq 0 ]; then
		echo "Test "$testnumber" (postproc) -> DONE"; 
	else	
		echo "Test "$testnumber" (postproc) -> FAILED - FAILED - FAILED - FAILED - FAILED"; 
		failedcases=1
	fi
	echo "########################################################################"
	numcases=$(( $numcases + 1 ))
	
	# format: dirname 1D/2D/3D mpi-ranks Lmax IC run/postproc
	./Cases/testcase_SOD_X.sh $DIRNAME 3 8 0 8 run
	if [ $? -eq 0 ]; then
		echo "Test "$testnumber" (sim) -> DONE"; 
	else	
		echo "Test "$testnumber" (sim) -> FAILED - FAILED - FAILED - FAILED - FAILED"; 
		failedcases=1
	fi
	./Cases/testcase_SOD_X.sh $DIRNAME 3 8 0 8 postproc
	if [ $? -eq 0 ]; then
		echo "Test "$testnumber" (postproc) -> DONE"; 
	else	
		echo "Test "$testnumber" (postproc) -> FAILED - FAILED - FAILED - FAILED - FAILED"; 
		failedcases=1
	fi
	echo "########################################################################"
	numcases=$(( $numcases + 1 ))

	# format: dirname 1D/2D/3D mpi-ranks Lmax IC run/postproc
	./Cases/testcase_REVERSE_SOD_X.sh $DIRNAME 3 1 0 8 run
	if [ $? -eq 0 ]; then
		echo "Test "$testnumber" (sim) -> DONE"; 
	else	
		echo "Test "$testnumber" (sim) -> FAILED - FAILED - FAILED - FAILED - FAILED"; 
		failedcases=1
	fi
	./Cases/testcase_REVERSE_SOD_X.sh $DIRNAME 3 1 0 8 postproc
	if [ $? -eq 0 ]; then
		echo "Test "$testnumber" (postproc) -> DONE"; 
	else	
		echo "Test "$testnumber" (postproc) -> FAILED - FAILED - FAILED - FAILED - FAILED"; 
		failedcases=1
	fi
	echo "########################################################################"
	numcases=$(( $numcases + 1 ))

	# format: dirname 1D/2D/3D mpi-ranks Lmax IC run/postproc
	./Cases/testcase_REVERSE_SOD_X.sh $DIRNAME 3 4 0 8 run
	if [ $? -eq 0 ]; then
		echo "Test "$testnumber" (sim) -> DONE"; 
	else	
		echo "Test "$testnumber" (sim) -> FAILED - FAILED - FAILED - FAILED - FAILED"; 
		failedcases=1
	fi
	./Cases/testcase_REVERSE_SOD_X.sh $DIRNAME 3 4 0 8 postproc
	if [ $? -eq 0 ]; then
		echo "Test "$testnumber" (postproc) -> DONE"; 
	else	
		echo "Test "$testnumber" (postproc) -> FAILED - FAILED - FAILED - FAILED - FAILED"; 
		failedcases=1
	fi
	echo "########################################################################"
	numcases=$(( $numcases + 1 ))
	
	# format: dirname 1D/2D/3D mpi-ranks Lmax IC run/postproc
	./Cases/testcase_REVERSE_SOD_X.sh $DIRNAME 3 8 0 8 run
	if [ $? -eq 0 ]; then
		echo "Test "$testnumber" (sim) -> DONE"; 
	else	
		echo "Test "$testnumber" (sim) -> FAILED - FAILED - FAILED - FAILED - FAILED"; 
		failedcases=1
	fi
	./Cases/testcase_REVERSE_SOD_X.sh $DIRNAME 3 8 0 8 postproc
	if [ $? -eq 0 ]; then
		echo "Test "$testnumber" (postproc) -> DONE"; 
	else	
		echo "Test "$testnumber" (postproc) -> FAILED - FAILED - FAILED - FAILED - FAILED"; 
		failedcases=1
	fi
	echo "########################################################################"
	numcases=$(( $numcases + 1 ))
fi
	
#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------
testnumber=301
if [ "$TESTCASENUMBER" == "$testnumber" ] || [ "$TESTCASENUMBER" == "all" ]; then
	echo "Test "$testnumber" started..."
	# format: dirname 1D/2D/3D mpi-ranks Lmax IC run/postproc
	./Cases/testcase_SOD_Y.sh $DIRNAME 3 1 0 8 run
	if [ $? -eq 0 ]; then
		echo "Test "$testnumber" (sim) -> DONE"; 
	else	
		echo "Test "$testnumber" (sim) -> FAILED - FAILED - FAILED - FAILED - FAILED"; 
		failedcases=1
	fi
	./Cases/testcase_SOD_Y.sh $DIRNAME 3 1 0 8 postproc
	if [ $? -eq 0 ]; then
		echo "Test "$testnumber" (postproc) -> DONE"; 
	else	
		echo "Test "$testnumber" (postproc) -> FAILED - FAILED - FAILED - FAILED - FAILED"; 
		failedcases=1
	fi
	echo "########################################################################"
	numcases=$(( $numcases + 1 ))

	# format: dirname 1D/2D/3D mpi-ranks Lmax IC run/postproc
	./Cases/testcase_SOD_Y.sh $DIRNAME 3 4 0 8 run
	if [ $? -eq 0 ]; then
		echo "Test "$testnumber" (sim) -> DONE"; 
	else	
		echo "Test "$testnumber" (sim) -> FAILED - FAILED - FAILED - FAILED - FAILED"; 
		failedcases=1
	fi
	./Cases/testcase_SOD_Y.sh $DIRNAME 3 4 0 8 postproc
	if [ $? -eq 0 ]; then
		echo "Test "$testnumber" (postproc) -> DONE"; 
	else	
		echo "Test "$testnumber" (postproc) -> FAILED - FAILED - FAILED - FAILED - FAILED"; 
		failedcases=1
	fi
	echo "########################################################################"
	numcases=$(( $numcases + 1 ))
	
	# format: dirname 1D/2D/3D mpi-ranks Lmax IC run/postproc
	./Cases/testcase_SOD_Y.sh $DIRNAME 3 8 0 8 run
	if [ $? -eq 0 ]; then
		echo "Test "$testnumber" (sim) -> DONE"; 
	else	
		echo "Test "$testnumber" (sim) -> FAILED - FAILED - FAILED - FAILED - FAILED"; 
		failedcases=1
	fi
	./Cases/testcase_SOD_Y.sh $DIRNAME 3 8 0 8 postproc
	if [ $? -eq 0 ]; then
		echo "Test "$testnumber" (postproc) -> DONE"; 
	else	
		echo "Test "$testnumber" (postproc) -> FAILED - FAILED - FAILED - FAILED - FAILED"; 
		failedcases=1
	fi
	echo "########################################################################"
	numcases=$(( $numcases + 1 ))

	# format: dirname 1D/2D/3D mpi-ranks Lmax IC run/postproc
	./Cases/testcase_REVERSE_SOD_Y.sh $DIRNAME 3 1 0 8 run
	if [ $? -eq 0 ]; then
		echo "Test "$testnumber" (sim) -> DONE"; 
	else	
		echo "Test "$testnumber" (sim) -> FAILED - FAILED - FAILED - FAILED - FAILED"; 
		failedcases=1
	fi
	./Cases/testcase_REVERSE_SOD_Y.sh $DIRNAME 3 1 0 8 postproc
	if [ $? -eq 0 ]; then
		echo "Test "$testnumber" (postproc) -> DONE"; 
	else	
		echo "Test "$testnumber" (postproc) -> FAILED - FAILED - FAILED - FAILED - FAILED"; 
		failedcases=1
	fi
	echo "########################################################################"
	numcases=$(( $numcases + 1 ))

	# format: dirname 1D/2D/3D mpi-ranks Lmax IC run/postproc
	./Cases/testcase_REVERSE_SOD_Y.sh $DIRNAME 3 4 0 8 run
	if [ $? -eq 0 ]; then
		echo "Test "$testnumber" (sim) -> DONE"; 
	else	
		echo "Test "$testnumber" (sim) -> FAILED - FAILED - FAILED - FAILED - FAILED"; 
		failedcases=1
	fi
	./Cases/testcase_REVERSE_SOD_Y.sh $DIRNAME 3 4 0 8 postproc
	if [ $? -eq 0 ]; then
		echo "Test "$testnumber" (postproc) -> DONE"; 
	else	
		echo "Test "$testnumber" (postproc) -> FAILED - FAILED - FAILED - FAILED - FAILED"; 
		failedcases=1
	fi
	echo "########################################################################"
	numcases=$(( $numcases + 1 ))
	
	# format: dirname 1D/2D/3D mpi-ranks Lmax IC run/postproc
	./Cases/testcase_REVERSE_SOD_Y.sh $DIRNAME 3 8 0 8 run
	if [ $? -eq 0 ]; then
		echo "Test "$testnumber" (sim) -> DONE"; 
	else	
		echo "Test "$testnumber" (sim) -> FAILED - FAILED - FAILED - FAILED - FAILED"; 
		failedcases=1
	fi
	./Cases/testcase_REVERSE_SOD_Y.sh $DIRNAME 3 8 0 8 postproc
	if [ $? -eq 0 ]; then
		echo "Test "$testnumber" (postproc) -> DONE"; 
	else	
		echo "Test "$testnumber" (postproc) -> FAILED - FAILED - FAILED - FAILED - FAILED"; 
		failedcases=1
	fi
	echo "########################################################################"
	numcases=$(( $numcases + 1 ))
fi	

#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------
testnumber=302
if [ "$TESTCASENUMBER" == "$testnumber" ] || [ "$TESTCASENUMBER" == "all" ]; then
	echo "Test "$testnumber" started..."
	# format: dirname 1D/2D/3D mpi-ranks Lmax IC run/postproc
	./Cases/testcase_SOD_Z.sh $DIRNAME 3 1 0 8 run
	if [ $? -eq 0 ]; then
		echo "Test "$testnumber" (sim) -> DONE"; 
	else	
		echo "Test "$testnumber" (sim) -> FAILED - FAILED - FAILED - FAILED - FAILED"; 
		failedcases=1
	fi
	./Cases/testcase_SOD_Z.sh $DIRNAME 3 1 0 8 postproc
	if [ $? -eq 0 ]; then
		echo "Test "$testnumber" (postproc) -> DONE"; 
	else	
		echo "Test "$testnumber" (postproc) -> FAILED - FAILED - FAILED - FAILED - FAILED"; 
		failedcases=1
	fi
	echo "########################################################################"
	numcases=$(( $numcases + 1 ))

	# format: dirname 1D/2D/3D mpi-ranks Lmax IC run/postproc
	./Cases/testcase_SOD_Z.sh $DIRNAME 3 4 0 8 run
	if [ $? -eq 0 ]; then
		echo "Test "$testnumber" (sim) -> DONE"; 
	else	
		echo "Test "$testnumber" (sim) -> FAILED - FAILED - FAILED - FAILED - FAILED"; 
		failedcases=1
	fi
	./Cases/testcase_SOD_Z.sh $DIRNAME 3 4 0 8 postproc
	if [ $? -eq 0 ]; then
		echo "Test "$testnumber" (postproc) -> DONE"; 
	else	
		echo "Test "$testnumber" (postproc) -> FAILED - FAILED - FAILED - FAILED - FAILED"; 
		failedcases=1
	fi
	echo "########################################################################"
	numcases=$(( $numcases + 1 ))
	
	# format: dirname 1D/2D/3D mpi-ranks Lmax IC run/postproc
	./Cases/testcase_SOD_Z.sh $DIRNAME 3 8 0 8 run
	if [ $? -eq 0 ]; then
		echo "Test "$testnumber" (sim) -> DONE"; 
	else	
		echo "Test "$testnumber" (sim) -> FAILED - FAILED - FAILED - FAILED - FAILED"; 
		failedcases=1
	fi
	./Cases/testcase_SOD_Z.sh $DIRNAME 3 8 0 8 postproc
	if [ $? -eq 0 ]; then
		echo "Test "$testnumber" (postproc) -> DONE"; 
	else	
		echo "Test "$testnumber" (postproc) -> FAILED - FAILED - FAILED - FAILED - FAILED"; 
		failedcases=1
	fi
	echo "########################################################################"
	numcases=$(( $numcases + 1 ))

	# format: dirname 1D/2D/3D mpi-ranks Lmax IC run/postproc
	./Cases/testcase_REVERSE_SOD_Z.sh $DIRNAME 3 1 0 8 run
	if [ $? -eq 0 ]; then
		echo "Test "$testnumber" (sim) -> DONE"; 
	else	
		echo "Test "$testnumber" (sim) -> FAILED - FAILED - FAILED - FAILED - FAILED"; 
		failedcases=1
	fi
	./Cases/testcase_REVERSE_SOD_Z.sh $DIRNAME 3 1 0 8 postproc
	if [ $? -eq 0 ]; then
		echo "Test "$testnumber" (postproc) -> DONE"; 
	else	
		echo "Test "$testnumber" (postproc) -> FAILED - FAILED - FAILED - FAILED - FAILED"; 
		failedcases=1
	fi
	echo "########################################################################"
	numcases=$(( $numcases + 1 ))

	# format: dirname 1D/2D/3D mpi-ranks Lmax IC run/postproc
	./Cases/testcase_REVERSE_SOD_Z.sh $DIRNAME 3 4 0 8 run
	if [ $? -eq 0 ]; then
		echo "Test "$testnumber" (sim) -> DONE"; 
	else	
		echo "Test "$testnumber" (sim) -> FAILED - FAILED - FAILED - FAILED - FAILED"; 
		failedcases=1
	fi
	./Cases/testcase_REVERSE_SOD_Z.sh $DIRNAME 3 4 0 8 postproc
	if [ $? -eq 0 ]; then
		echo "Test "$testnumber" (postproc) -> DONE"; 
	else	
		echo "Test "$testnumber" (postproc) -> FAILED - FAILED - FAILED - FAILED - FAILED"; 
		failedcases=1
	fi
	echo "########################################################################"
	numcases=$(( $numcases + 1 ))
	
	# format: dirname 1D/2D/3D mpi-ranks Lmax IC run/postproc
	./Cases/testcase_REVERSE_SOD_Z.sh $DIRNAME 3 8 0 8 run
	if [ $? -eq 0 ]; then
		echo "Test "$testnumber" (sim) -> DONE"; 
	else	
		echo "Test "$testnumber" (sim) -> FAILED - FAILED - FAILED - FAILED - FAILED"; 
		failedcases=1
	fi
	./Cases/testcase_REVERSE_SOD_Z.sh $DIRNAME 3 8 0 8 postproc
	if [ $? -eq 0 ]; then
		echo "Test "$testnumber" (postproc) -> DONE"; 
	else	
		echo "Test "$testnumber" (postproc) -> FAILED - FAILED - FAILED - FAILED - FAILED"; 
		failedcases=1
	fi
	echo "########################################################################"
	numcases=$(( $numcases + 1 ))
fi	

if [ $numcases -eq 0 ]; then
	echo "Non-existing testcase chosen..."; 
	exit 1
fi

if [ $failedcases -ne 0 ]; then
	echo "WARNING!!!!! Not all testcases passed!"; 
	exit 1
else	
	echo "All testcases passed!"; 
	exit 0
fi

