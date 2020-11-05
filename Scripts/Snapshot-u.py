# trace generated using paraview version 5.6.0
#
# To ensure correct image size when batch processing, please search 
# for and uncomment the line `# renderView*.ViewSize = [*,*]`

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'PVD Reader'
obvnodampvd = PVDReader(FileName='/import/freenas-m-05-seissol/dli/ExaHyPE/Obv/test2/obv-dam.pvd')
obvnodampvd.CellArrays = ['time', 'Q']

# get animation scene
animationScene1 = GetAnimationScene()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1482, 799]

# show data in view
obvnodampvdDisplay = Show(obvnodampvd, renderView1)

# trace defaults for the display properties.
obvnodampvdDisplay.Representation = 'Surface'
obvnodampvdDisplay.ColorArrayName = [None, '']
obvnodampvdDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
obvnodampvdDisplay.SelectOrientationVectors = 'None'
obvnodampvdDisplay.ScaleFactor = 8000.0
obvnodampvdDisplay.SelectScaleArray = 'None'
obvnodampvdDisplay.GlyphType = 'Arrow'
obvnodampvdDisplay.GlyphTableIndexArray = 'None'
obvnodampvdDisplay.GaussianRadius = 400.0
obvnodampvdDisplay.SetScaleArray = [None, '']
obvnodampvdDisplay.ScaleTransferFunction = 'PiecewiseFunction'
obvnodampvdDisplay.OpacityArray = [None, '']
obvnodampvdDisplay.OpacityTransferFunction = 'PiecewiseFunction'
obvnodampvdDisplay.DataAxesGrid = 'GridAxesRepresentation'
obvnodampvdDisplay.SelectionCellLabelFontFile = ''
obvnodampvdDisplay.SelectionPointLabelFontFile = ''
obvnodampvdDisplay.PolarAxes = 'PolarAxesRepresentation'
obvnodampvdDisplay.ScalarOpacityUnitDistance = 2564.382088535499

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
obvnodampvdDisplay.DataAxesGrid.XTitleFontFile = ''
obvnodampvdDisplay.DataAxesGrid.YTitleFontFile = ''
obvnodampvdDisplay.DataAxesGrid.ZTitleFontFile = ''
obvnodampvdDisplay.DataAxesGrid.XLabelFontFile = ''
obvnodampvdDisplay.DataAxesGrid.YLabelFontFile = ''
obvnodampvdDisplay.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
obvnodampvdDisplay.PolarAxes.PolarAxisTitleFontFile = ''
obvnodampvdDisplay.PolarAxes.PolarAxisLabelFontFile = ''
obvnodampvdDisplay.PolarAxes.LastRadialAxisTextFontFile = ''
obvnodampvdDisplay.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# reset view to fit data
renderView1.ResetCamera()

# get the material library
materialLibrary1 = GetMaterialLibrary()

# update the view to ensure updated data information
renderView1.Update()

# set scalar coloring
ColorBy(obvnodampvdDisplay, ('CELLS', 'Q', '20'))

# rescale color and/or opacity maps used to include current data range
obvnodampvdDisplay.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
obvnodampvdDisplay.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'Q'
qLUT = GetColorTransferFunction('Q')

# get opacity transfer function/opacity map for 'Q'
qPWF = GetOpacityTransferFunction('Q')

animationScene1.GoToLast()

# set scalar coloring
ColorBy(obvnodampvdDisplay, ('CELLS', 'Q', '1'))

# rescale color and/or opacity maps used to exactly fit the current data range
obvnodampvdDisplay.RescaleTransferFunctionToDataRange(False, False)

# Update a scalar bar component title.
UpdateScalarBarsComponentTitle(qLUT, obvnodampvdDisplay)

# rescale color and/or opacity maps used to exactly fit the current data range
obvnodampvdDisplay.RescaleTransferFunctionToDataRange(False, True)

# rescale color and/or opacity maps used to exactly fit the current data range
obvnodampvdDisplay.RescaleTransferFunctionToDataRange(False, True)

# rescale color and/or opacity maps used to exactly fit the current data range
obvnodampvdDisplay.RescaleTransferFunctionToDataRange(False, True)

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.0, 0.0, 218564.06460551018]
renderView1.CameraParallelScale = 56568.5424949238

# save screenshot
SaveScreenshot('/import/freenas-m-05-seissol/dli/ExaHyPE/Obv/dam-u.jpg', renderView1, ImageResolution=[1482, 799])

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.0, 0.0, 218564.06460551018]
renderView1.CameraParallelScale = 56568.5424949238

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).
