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

def Riemann(points, tend, Gamma, rholeft, rhoright, uleft, uright, pleft, pright,):
	#print 'Analytic results:'
	# Assumed structure of exact solution
	#
	#    \         /      |con |       |s|
	#     \   f   /       |tact|       |h|
	# left \  a  /  state |disc| state |o| right
	# state \ n /    2    |cont|   3   |c| state
	#   1    \ /          |tinu|       |k|   4
	#         |           |ity |       | |
	#	    x = 0.5

	# HEADER
	# inputs
	global gamma
	gamma = Gamma

	cellmid = []
	for i in range(0, len(points)):
		cellmid.append(float(points[i]))

	xinit = 1.5 # initial solution for Newton's method
	error = 1e-15 # for Newton's method

	## Riemann solver
	cleft = (gamma * pleft / rholeft)**(0.5)
	cright = (gamma * pright / rhoright)**(0.5)

	global PRL, CRL
	PRL = pright / pleft
	CRL = cright / cleft

	# Call Newton's method to solve shock number
	from sod_RiemannSolverNewton import solverMs
	Ms = solverMs(PRL, CRL, gamma, xinit, error)
	
	#print "Shock number: %f" % Ms

	p34 = 1 + 2 * gamma / (gamma + 1) * (Ms**2 - 1)
	#print "p3/p4: %f" % p34

	p3 = p34 * pright

	alpha = (gamma + 1) / (gamma - 1)
	rho3 = rhoright * (1 + alpha * p34) / (alpha + p34)
	rho2 = rholeft * (p34 * pright / pleft)**(1 / gamma)
	u2 = uleft - uright + (2 / (gamma - 1)) * cleft * (1 - (p34 * pright / pleft)**((gamma - 1) / (2 * gamma)))
	c2 = (gamma * p3 / rho2)**(0.5)

	# Shock position
	spos = 0.5 + tend * cright * ((gamma-1)/(2*gamma) + (gamma+1)/(2*gamma)*p34)**0.5 + tend * uright
	# Position of contact discontinuity 
	conpos = 0.5 + u2 * tend + tend * uright;	
	# Start of expansion fan
	pos1 = 0.5 + (uleft - cleft)*tend;	
	# End of expansion fan
	pos2 = 0.5 + (u2+uright-c2)*tend;	

	#print "Shock position: %f" % spos
	#print "Contact disc. position: %f" % conpos
	#print "Start exp. fan: %f" % pos1
	#print "End exp. fan: %f" % pos2


	## check
	pexact = []
	uexact = []
	rhoexact = []
	machexact = []
	cexact = []

	for ii in range(0, len(cellmid)):
		xx = cellmid[ii]	
		if xx <= pos1:
			pexact.append(pleft)
			rhoexact.append(rholeft)
			uexact.append(uleft)
			cexact.append(((gamma * pexact[ii]) / rhoexact[ii])**0.5)
			machexact.append(uexact[ii] / cexact[ii])
		elif xx <= pos2:
			pexact.append(pleft * (1 + (pos1 - xx) / (cleft * alpha * tend))**(2 * gamma / (gamma - 1)))
			rhoexact.append(rholeft * (1 + (pos1 - xx) / (cleft * alpha * tend))**(2 / (gamma - 1)))
			uexact.append(uleft + (2 /(gamma + 1))*(xx - pos1) / tend)
			cexact.append((gamma * pexact[ii] / rhoexact[ii])**0.5)
			machexact.append(uexact[ii] / cexact[ii])
		elif xx <= conpos:
			pexact.append(p3)
			rhoexact.append(rho2)
			uexact.append(u2 + uright)
			cexact.append((gamma * pexact[ii] / rhoexact[ii])**0.5)
			machexact.append(uexact[ii] / cexact[ii])
		elif xx <= spos:
			pexact.append(p3)
			rhoexact.append(rho3)
			uexact.append(u2 + uright)
			cexact.append((gamma * pexact[ii] / rhoexact[ii])**0.5)
			machexact.append(uexact[ii] / cexact[ii])
		else:
			pexact.append(pright)
			rhoexact.append(rhoright)
			uexact.append(uright)
			cexact.append((gamma * pexact[ii] / rhoexact[ii])**0.5)
			machexact.append(uexact[ii] / cexact[ii])

		# print xx, rhoexact[ii]
		
	### Write out 
#	target = open("dens.txt", 'w')
#	target.truncate()
#	for ii in range(0, len(posvec)):
	#	target.write(str(0.5 * posvec[ii] + 0.5) + ", " + str(rhoexact[ii])  + "\n")
#		target.write(str(posvec[ii]) + ", " + str(rhoexact[ii])  + "\n")

#	target.close()

#	with open('dens_analytic.csv','w') as csvfile:
#		csvfile.write("x, y, z, Density Analytic \n")
#		for ii in range(0, len(posvec) ):
#			csvfile.write(str(posvec[ii]) + ", 0.125, 0.125, " + str(rhoexact[ii])  + "\n")	
			
	return [rhoexact,uexact]
		
