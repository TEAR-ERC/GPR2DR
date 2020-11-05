# trace generated using paraview version 5.6.0
#
# To ensure correct image size when batch processing, please search 
# for and uncomment the line `# renderView*.ViewSize = [*,*]`

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'PVD Reader'
obvdampvd = PVDReader(FileName='/import/freenas-m-05-seissol/dli/ExaHyPE/Obv/test2//obv-dam.pvd')
obvdampvd.CellArrays = ['time', 'Q']

# get animation scene
animationScene1 = GetAnimationScene()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1482, 799]

# show data in view
obvdampvdDisplay = Show(obvdampvd, renderView1)

# trace defaults for the display properties.
obvdampvdDisplay.Representation = 'Surface'
obvdampvdDisplay.ColorArrayName = [None, '']
obvdampvdDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
obvdampvdDisplay.SelectOrientationVectors = 'None'
obvdampvdDisplay.ScaleFactor = 4000.0
obvdampvdDisplay.SelectScaleArray = 'None'
obvdampvdDisplay.GlyphType = 'Arrow'
obvdampvdDisplay.GlyphTableIndexArray = 'None'
obvdampvdDisplay.GaussianRadius = 200.0
obvdampvdDisplay.SetScaleArray = [None, '']
obvdampvdDisplay.ScaleTransferFunction = 'PiecewiseFunction'
obvdampvdDisplay.OpacityArray = [None, '']
obvdampvdDisplay.OpacityTransferFunction = 'PiecewiseFunction'
obvdampvdDisplay.DataAxesGrid = 'GridAxesRepresentation'
obvdampvdDisplay.SelectionCellLabelFontFile = ''
obvdampvdDisplay.SelectionPointLabelFontFile = ''
obvdampvdDisplay.PolarAxes = 'PolarAxesRepresentation'
obvdampvdDisplay.ScalarOpacityUnitDistance = 1162.835927399163

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
obvdampvdDisplay.DataAxesGrid.XTitleFontFile = ''
obvdampvdDisplay.DataAxesGrid.YTitleFontFile = ''
obvdampvdDisplay.DataAxesGrid.ZTitleFontFile = ''
obvdampvdDisplay.DataAxesGrid.XLabelFontFile = ''
obvdampvdDisplay.DataAxesGrid.YLabelFontFile = ''
obvdampvdDisplay.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
obvdampvdDisplay.PolarAxes.PolarAxisTitleFontFile = ''
obvdampvdDisplay.PolarAxes.PolarAxisLabelFontFile = ''
obvdampvdDisplay.PolarAxes.LastRadialAxisTextFontFile = ''
obvdampvdDisplay.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# reset view to fit data
renderView1.ResetCamera()

#changing interaction mode based on data extents
renderView1.CameraPosition = [0.0, 0.0, 10000.0]

# get the material library
materialLibrary1 = GetMaterialLibrary()

# update the view to ensure updated data information
renderView1.Update()

# set scalar coloring
ColorBy(obvdampvdDisplay, ('CELLS', 'Q', '1'))

# rescale color and/or opacity maps used to include current data range
obvdampvdDisplay.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
obvdampvdDisplay.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'Q'
qLUT = GetColorTransferFunction('Q')

# get opacity transfer function/opacity map for 'Q'
qPWF = GetOpacityTransferFunction('Q')

# set scalar coloring
ColorBy(obvdampvdDisplay, ('CELLS', 'Q', '20'))

# rescale color and/or opacity maps used to exactly fit the current data range
obvdampvdDisplay.RescaleTransferFunctionToDataRange(False, False)

# Update a scalar bar component title.
UpdateScalarBarsComponentTitle(qLUT, obvdampvdDisplay)

animationScene1.GoToLast()

# Rescale transfer function
qLUT.RescaleTransferFunction(0.0, 0.1)

# Rescale transfer function
qPWF.RescaleTransferFunction(0.0, 0.1)

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.0, 0.0, 10000.0]
renderView1.CameraParallelScale = 28284.2712474619

# save screenshot
SaveScreenshot('/import/freenas-m-05-seissol/dli/ExaHyPE/Obv/dam-xi.jpg', renderView1, ImageResolution=[1482, 799])

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.0, 0.0, 10000.0]
renderView1.CameraParallelScale = 28284.2712474619

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).
