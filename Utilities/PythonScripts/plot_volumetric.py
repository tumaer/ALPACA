from paraview.simple import *
import sys
import numpy as np


print(str(sys.argv))
time = sys.argv[1]
print(str(time))

resolution_x = 2835
resolution_y = 2304
#resolution_x = 5760
#resolution_y = 4608
factor = 350.0
max_scaled_value = np.log( factor + 1.0 )

specular = 90.0
luminosity = 90.0
ambient = 15.0

#pictureName = './baro.pdf'
#dataFile = './custom_0.104278.xdmf'
#interfaceFile = './interface_0.104278.xdmf'
dataFile = './custom_' + str(time) + '.xdmf'
interfaceFile = './interface_' + str(time) + '.xdmf'
pictureName = './vorti' + str(time) + '.pdf'
#volume_array = 'vortex_stretching'
volume_array = 'vorticity'
#volume_array = 'baroclinicity'

interfaceXdmf = XDMFReader(FileNames=[interfaceFile])
interfaceXdmf.CellArrayStatus = ['levelset']
interfaceXdmf.GridStatus = ['SpatialData']
interfaceXdmf = CleantoGrid(Input=interfaceXdmf)

interfacePointData = CellDatatoPointData(Input=interfaceXdmf)
interfacePointData.UpdatePipeline()

interfaceContour = Contour(Input=interfacePointData)
interfaceContour.ContourBy = ['POINTS', 'levelset']
interfaceContour.Isosurfaces = [0.0]
#rep = GetDisplayProperties(interfaceContour)
#rep.DiffuseColor = [1.0, 1.0, 1.0]


fieldXdmf = XDMFReader(FileNames=[dataFile])
fieldXdmf.CellArrayStatus = [volume_array]
fieldXdmf.GridStatus = ['SpatialData']

scaledVorticity = PythonCalculator(Input=fieldXdmf)
expression_string = 'ln( ' + str(factor) + ' * ' + volume_array + ' / max(' + volume_array + ') + 1.0 )'
scaledVorticity.Expression = expression_string
scaledVorticity.ArrayAssociation = 'Cell Data'
scaledVorticity.ArrayName = 'scaled_vorticity'
scaledVorticity.UpdatePipeline()

reducedDomain = IsoVolume(Input=scaledVorticity)
reducedDomain.InputScalars = ['CELLS', 'scaled_vorticity']
reducedDomain.ThresholdRange = [0.1 * max_scaled_value, 1.1 * max_scaled_value]


renderView1 = GetActiveViewOrCreate('RenderView')
renderView1.ViewSize = [resolution_x, resolution_y]
renderView1.OrientationAxesVisibility = 0
SetViewProperties(Background=[0.05, 0.05, 0.05])

cam = GetActiveCamera()
cam.SetPosition(-0.1, 0.45, 0.4)
cam.SetFocalPoint(0.21, 0.1, 0.1)
#cam.SetPosition(0.8, 0.32, 0.51)
#cam.SetFocalPoint(0.2, 0.085, 0.13)
#cam.SetPosition(-0.23, -0.15, -0.25)
#cam.SetFocalPoint(0.2, 0.085, 0.13)
cam.OrthogonalizeViewUp()

contourDisplay = Show(interfaceContour, renderView1)
contourDisplay.SetRepresentationType('Volume')
contourDisplay.SetScalarBarVisibility(renderView1, False)
contourDisplay.Opacity = 0.25
contourDisplay.Specular = specular / 100.0
contourDisplay.SpecularPower = specular
contourDisplay.Luminosity = luminosity
contourDisplay.Ambient = ambient / 100.0
contourDisplay.Diffuse = 0.1

isoVolumeDisplay = Show(reducedDomain, renderView1)
isoVolumeDisplay.SetScalarBarVisibility(renderView1, False)
isoVolumeDisplay.Opacity = 0.5
isoVolumeDisplay.Specular = specular / 100.0
isoVolumeDisplay.SpecularPower = specular
isoVolumeDisplay.Luminosity = luminosity
isoVolumeDisplay.Ambient = ambient / 100.0
isoVolumeDisplay.Diffuse = 0.1

ColorBy(isoVolumeDisplay, ('CELLS', 'scaled_vorticity'))
isoVolumeDisplay.RescaleTransferFunctionToDataRange(True, False)
scaled_vorticityLUT = GetColorTransferFunction('scaled_vorticity')
scaled_vorticityPWF = GetOpacityTransferFunction('scaled_vorticity')

scaled_vorticityLUT.ApplyPreset('erdc_rainbow_bright', True)
scaled_vorticityPWF.Points = [0.0, 0.0, 0.5, 0.0,
                              0.6 * max_scaled_value, 0.1, 0.5, 0.0,
                              0.8 * max_scaled_value, 0.81, 0.5, 0.0,
                              1.0 * max_scaled_value, 1.0, 0.5, 0.0]

#ResetCamera()
isoVolumeDisplay.SetRepresentationType('Volume')
#SaveScreenshot(pictureName, renderView1, ImageResolution=[1280, 1000])
ExportView(pictureName, view=renderView1)


#renderView1.ResetCamera()




