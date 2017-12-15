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

from sys import argv
import csv
import math
from decimal import *

script, dim, time, maxlevel, blocksize, ic, sod_direction = argv
dim = int(dim)
time = float(time)
maxlevel = int(maxlevel)
blocksize = float(blocksize)
ic = int(ic)
sod_direction = int(sod_direction)

# Open Point Data (x, y, z)
with open('cellcenterdata.csv', 'rb') as f:
	reader = csv.reader(f,delimiter=',')
	position = []
	density  = []
	velocity = []
	next(reader)
	for row in reader:
		pos = row[4+dim-1+sod_direction]
		dens = row[1]
		u = float(row[3+sod_direction])/float(dens)

		position.append(pos)
		density.append(dens)
		velocity.append(u)
f.close()

n = len(position)

# Analytic results

from sod_RiemannSolver import Riemann
[rhoexact,uexact] = Riemann(position, time, 1.4, 1.0,0.125, 0.0,0.0, 1.0,0.1)

dx = []
for j in range(0,maxlevel+1):
	dx.append(float(blocksize/ic * (float(1)/2)**j))

threshold = 1.e-10
cellvolume = []
for j in range(0,n):
	for jj in range(0, maxlevel+1):
		mod = Decimal(float(position[j])+0.5*dx[jj]) % Decimal(float(dx[jj]))	

		if mod < threshold:
			break
	
	cellvolume.append(dx[jj]**dim)
 

# Calculate Errors
from errorNorms import errorNorms_abs
from errorNorms import errorNorms_rel
[rhoL1,rhoL2,rhoLinf] = errorNorms_rel(cellvolume,density,rhoexact)
[uL1,uL2,uLinf] = errorNorms_abs(cellvolume,velocity,uexact)

print "rel. L1 | L2 | Linf (rho) : %.8e %.8e %.8e" %(rhoL1,rhoL2,rhoLinf)
print "abs. L1 | L2 | Linf ( u ) : %.8e %.8e %.8e" %(uL1,uL2,uLinf)
print "volume : %f" %(sum(cellvolume)) 

target = open('errornorms.txt','w')
target.truncate()
target.write("#Error norms: L1, L2, Linf \n")
target.write("%s %.8e %.8e %.8e\n" %("rho",rhoL1,rhoL2,rhoLinf))
target.write("%s %.8e %.8e %.8e\n" %(" u ",uL1,uL2,uLinf))
target.close

target = open('pos_rho_vel_sod.txt','w')
target.truncate()
target.write("# x rho rho-analytical velocity velocity-analytical \n")
for ii in range(0,len(position)):
	target.write("%.8e %.8e %.8e %.8e %.8e\n" %(float(position[ii]),float(density[ii]),float(rhoexact[ii]),float(velocity[ii]),float(uexact[ii])))

target.close()
