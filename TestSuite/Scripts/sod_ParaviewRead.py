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

# import the simple module from the paraview
from paraview.simple import *
# import sys
import sys

resname = sys.argv[1]
sourcepath = "domain/" + resname

#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'XDMF Reader'
data = XDMFReader(FileNames=[sourcepath])
data.CellArrayStatus = ['density', 'energy', 'partition', 'x_momentum', 'y_momentum', 'z_momentum']

# Properties modified on data
#data.GridStatus = ['step_4']

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
#renderView1.ViewSize = [3004, 1616]

# get color transfer function/color map for 'partition'
partitionLUT = GetColorTransferFunction('partition')

# show data in view
dataDisplay = Show(data, renderView1)
# trace defaults for the display properties.
dataDisplay.ColorArrayName = ['CELLS', 'partition']
dataDisplay.LookupTable = partitionLUT
#dataDisplay.OSPRayScaleArray = 'partition'
#dataDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
#dataDisplay.GlyphType = 'Arrow'
#dataDisplay.ScalarOpacityUnitDistance = 0.029807488393148774
#dataDisplay.SetScaleArray = [None, '']
#dataDisplay.ScaleTransferFunction = 'PiecewiseFunction'
#dataDisplay.OpacityArray = [None, '']
#dataDisplay.OpacityTransferFunction = 'PiecewiseFunction'

# reset view to fit data
renderView1.ResetCamera()

# show color bar/color legend
dataDisplay.SetScalarBarVisibility(renderView1, True)

# get opacity transfer function/opacity map for 'partition'
partitionPWF = GetOpacityTransferFunction('partition')

# set scalar coloring
ColorBy(dataDisplay, ('CELLS', 'density'))

# rescale color and/or opacity maps used to include current data range
dataDisplay.RescaleTransferFunctionToDataRange(True)

# show color bar/color legend
dataDisplay.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'density'
densityLUT = GetColorTransferFunction('density')

# get opacity transfer function/opacity map for 'density'
densityPWF = GetOpacityTransferFunction('density')

# rescale color and/or opacity maps used to exactly fit the current data range
dataDisplay.RescaleTransferFunctionToDataRange(False)

# change representation type
dataDisplay.SetRepresentationType('Surface With Edges')

# current camera placement for renderView1
renderView1.CameraPosition = [-0.2640603688866507, 1.6946332667063804, 1.1978561503111311]
renderView1.CameraFocalPoint = [0.4999999999999996, 0.12499999999999997, 0.12499999999999989]
renderView1.CameraViewUp = [0.2613690381647314, 0.6281547893712759, -0.732876378715245]
renderView1.CameraParallelScale = 0.5303300858899106

# save screenshot
SaveScreenshot('density.png', magnification=2, quality=100, view=renderView1)

# create a new 'Cell Centers'
cellCenters1 = CellCenters(Input=data)

# show data in view
cellCenters1Display = Show(cellCenters1, renderView1)
# trace defaults for the display properties.
cellCenters1Display.ColorArrayName = ['POINTS', 'partition']
cellCenters1Display.LookupTable = partitionLUT
#cellCenters1Display.OSPRayScaleArray = 'partition'
#cellCenters1Display.OSPRayScaleFunction = 'PiecewiseFunction'
#cellCenters1Display.GlyphType = 'Arrow'
#cellCenters1Display.SetScaleArray = ['POINTS', 'partition']
#cellCenters1Display.ScaleTransferFunction = 'PiecewiseFunction'
#cellCenters1Display.OpacityArray = ['POINTS', 'partition']
#cellCenters1Display.OpacityTransferFunction = 'PiecewiseFunction'

# hide data in view
Hide(data, renderView1)

# show color bar/color legend
cellCenters1Display.SetScalarBarVisibility(renderView1, True)

# save data
SaveData('cellcenterdata.csv', proxy=cellCenters1, Precision=10)

#### saving camera placements for all active views

# current camera placement for renderView1
#renderView1.CameraPosition = [0.5, 0.125, 2.174038105676658]
#renderView1.CameraFocalPoint = [0.5, 0.125, 0.125]
#renderView1.CameraParallelScale = 0.5303300858899106

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).

