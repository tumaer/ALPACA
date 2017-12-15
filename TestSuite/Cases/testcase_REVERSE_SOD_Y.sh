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
# Running sod test 

# arguments:
# $1: directoryname
# $2: dimensions
# $3: number of mpi-ranks
# $4: Lmax
# $5: InternalCells
# $6: run or postproc

DIRNAME=$1
DIM=$2
RANKS=$3
MAXLEVEL=$4
IC=$5
RUNORPOSTPROC=$6

FNAME="reverse_sod_y_"$2"D_L"${MAXLEVEL}"_IC"${IC}
BINNAME="ALPACA_"$2"D_IC"${IC}
if [ $RANKS -gt 1 ]; then
	FNAME=$FNAME"_"$RANKS"ranks"
fi

if [ "$RUNORPOSTPROC" == "run" ]; then
	echo "SOD SHOCK-TUBE IN REVERSE Y-DIRECTION using "${BINNAME}
	echo "Dimensions     : "$DIM
	if [ $RANKS -gt 1 ]; then
		echo "MPI            : "$RANKS" ranks"
	fi
	echo "Max. level     : "$MAXLEVEL
	echo "Internal cells : "$IC
else
	echo "Postprocessing..."
fi

if [ "$RUNORPOSTPROC" == "run" ]; then
	# 2 - Copy input.xml and adjust parameter
	cp -r ./Cases/sod_template.xml ./$DIRNAME/$FNAME.xml

	# Substitute Maxlevel in inputfile
	sed -i "s/MYMAXLEVEL/$MAXLEVEL/g" ./$DIRNAME/$FNAME.xml
	sed -i "s/MYXBLOCKRATIO/1/g" ./$DIRNAME/$FNAME.xml
	sed -i "s/MYYBLOCKRATIO/4/g" ./$DIRNAME/$FNAME.xml
	sed -i "s/MYZBLOCKRATIO/1/g" ./$DIRNAME/$FNAME.xml

	# Reverse direction of sod-case in inputfile
	sed -i '0,/if (x &lt; 0.5)/{s/(x/(y/}' ./$DIRNAME/$FNAME.xml
	sed -i '0,/if (y &lt; 0.5)/{s/lt/gt/}' ./$DIRNAME/$FNAME.xml

	# 4 - Start simulation. 
	cd ./$DIRNAME
	rm -rf ./$FNAME
	mpiexec -np $RANKS ./${BINNAME} $FNAME.xml > /dev/null

else 
	# Get blocksize from inputfile
	BLOCKSIZE=$(echo 'cat //blockSize' | xmllint --shell ./$DIRNAME/$FNAME.xml | awk '/<blockSize>/{print $(NF-1)}')

	# this mkdir seems to be required on the cluster only for some buffer-flush reasons
        # I'm not really sure what the problem is, but this command fixes it ;-)
	mkdir -p $DIRNAME

	cd ./$DIRNAME/$FNAME/domain
	RESNAME=$(ls -t data*.xdmf | head -1)
	cd ..

	# 5 - Start testsuite routine.
	# a.) Call Paraview to extract data
	$PARAVIEW_BASE/bin/pvpython ../../../Scripts/sod_ParaviewRead.py $RESNAME

	# b.) Sort, Calculate analytic result, Calculate Error norms
	TIME=${RESNAME%.*}
	TIME=${TIME#*_}
	python ../../../Scripts/reverse_sod_postprocess.py $DIM $TIME $MAXLEVEL $BLOCKSIZE $IC 1

	# check zero's in other momentum components
	xmom=0
	zmom=0
	if [ $DIM -eq 2 ]; then
		xmom=$(awk -F ',' 'NR>0 {total += sqrt($4*$4)} END{if (total > 0.) print "1"; else print "0"}' cellcenterdata.csv)
	elif [ $DIM -eq 3 ]; then
		xmom=$(awk -F ',' 'NR>0 {total += sqrt($4*$4)} END{if (total > 0.) print "1"; else print "0"}' cellcenterdata.csv)
		zmom=$(awk -F ',' 'NR>0 {total += sqrt($6*$6)} END{if (total > 0.) print "1"; else print "0"}' cellcenterdata.csv)
	fi


	# c.) Plot
	gnuplot "../../../Scripts/sod_plot_rho.gnu" 
	gnuplot "../../../Scripts/sod_plot_vel.gnu"

	rhoL1ok=1
	rhoL2ok=1
	rhoLinfok=1
	uL1ok=1
	uL2ok=1
	uLinfok=1

	rhoOK=(1 1 1)
	uOK=(1 1 1)

	# default result values to check the proper result
	err=0.0
	
	if [ $DIM -eq 2 ] && [ $MAXLEVEL -eq 0 ] && [ $IC -eq 8 ]; then
		rhorefs=(1.13815604e-02 5.24661980e-02 5.24822376e-01)
		urefs=(9.03233410e-03 5.00993901e-02 5.32644864e-01)

	elif [ $DIM -eq 2 ] && [ $MAXLEVEL -eq 3 ] && [ $IC -eq 16 ]; then
		rhorefs=(7.63509680e-04 9.24835979e-03 2.22899238e-01)
		urefs=(4.55147257e-04 7.48092001e-03 2.55804507e-01)

	elif [ $DIM -eq 2 ] && [ $MAXLEVEL -eq 2 ] && [ $IC -eq 32 ]; then
		rhorefs=(9.79945201e-04 1.65495149e-02 6.62604079e-01)
		urefs=(6.30265446e-04 1.46440785e-02 6.36290160e-01)

	elif [ $DIM -eq 3 ] && [ $MAXLEVEL -eq 0 ] && [ $IC -eq 8 ]; then
		rhorefs=(3.08498479e-03 3.04191947e-02 6.25914055e-01)
		urefs=(2.44261352e-03 2.85170544e-02 6.09101438e-01)
	fi


	rhoOK[0]=$(awk -v ref="${rhorefs[0]}" -v err="$err" 'NR==2 {if (sqrt(($2/ref-1.)^2) <= err) print 0; else print 1}' errornorms.txt)
	rhoOK[1]=$(awk -v ref="${rhorefs[1]}" -v err="$err" 'NR==2 {if (sqrt(($3/ref-1.)^2) <= err) print 0; else print 1}' errornorms.txt)
	rhoOK[2]=$(awk -v ref="${rhorefs[2]}" -v err="$err" 'NR==2 {if (sqrt(($4/ref-1.)^2) <= err) print 0; else print 1}' errornorms.txt)

	uOK[0]=$(awk -v ref="${urefs[0]}" -v err="$err" 'NR==3 {if (sqrt(($2/ref-1.)^2) <= err) print 0; else print 1}' errornorms.txt)
	uOK[1]=$(awk -v ref="${urefs[1]}" -v err="$err" 'NR==3 {if (sqrt(($3/ref-1.)^2) <= err) print 0; else print 1}' errornorms.txt)
	uOK[2]=$(awk -v ref="${urefs[2]}" -v err="$err" 'NR==3 {if (sqrt(($4/ref-1.)^2) <= err) print 0; else print 1}' errornorms.txt)

	echo "------------------------------------------------------------------------"
	printf "REFERENCE VALUES    (rho) : %s %s %s\n" $(awk -v val1="${rhorefs[0]}" -v val2="${rhorefs[1]}" -v val3="${rhorefs[2]}" 'BEGIN{printf "%.8e %.8e %.8e",val1,val2,val3}')
	printf "REFERENCE VALUES    ( u ) : %s %s %s\n" $(awk -v val1="${urefs[0]}" -v val2="${urefs[1]}" -v val3="${urefs[2]}" 'BEGIN{printf "%.8e %.8e %.8e",val1,val2,val3}')
	echo "------------------------------------------------------------------------"

	if [ "${rhoOK[0]}" == "" ]; then
		echo "Problem in L1-norm of density"
		exit 1
	fi
	if [ "${rhoOK[1]}" == ""  ]; then
		echo "Problem in L2-norm of density"
		exit 1
	fi
	if [ "${rhoOK[2]}" == ""  ]; then
		echo "Problem in Linf-norm of density"
		exit 1
	fi	

	
	if [ ${rhoOK[0]} -eq 1 ]; then
		echo "Problem in L1-norm of density"
		exit 1
	fi	
	if [ ${rhoOK[1]} -eq 1 ]; then
		echo "Problem in L2-norm of density"
		exit 1
	fi	
	if [ ${rhoOK[2]} -eq 1 ]; then
		echo "Problem in Linf-norm of density"
		exit 1
	fi


	if [ "${uOK[0]}" == "" -o ${uOK[0]} -eq 1 ]; then
		echo "Problem in L1-norm of velocity"
		exit 1
	fi
	if [ "${uOK[0]}" == "" -o ${uOK[0]} -eq 1 ]; then
		echo "Problem in L1-norm of velocity"
		exit 1
	fi
	
	if [ "${uOK[1]}" == "" -o ${uOK[1]} -eq 1 ]; then
		echo "Problem in L2-norm of velocity"
		exit 1
	fi

	if [ "${uOK[2]}" == "" -o ${uOK[2]} -eq 1 ]; then
		echo "Problem in Linf-norm of velocity"
		exit 1
	fi

	if [ "$xmom" == "" -o $xmom -ne 0 ]; then
		echo "Problem in y-component of velocity (non-zero)"
		echo $xmom
		exit 1
	fi

	if [ "$zmom" == "" -o $zmom -ne 0 ]; then
		echo "Problem in z-component of velocity (non-zero)"
		echo $zmom
		exit 1
	fi
fi
exit 0

