# trace generated using paraview version 5.6.0
#
# To ensure correct image size when batch processing, please search 
# for and uncomment the line `# renderView*.ViewSize = [*,*]`

#### import the simple module from the paraview
from paraview.simple import *

#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'PVD Reader'
conservedpvd = PVDReader(FileName='/import/freenas-m-05-seissol/dli/ExaHyPE/Obv/nodam/obv-nodam.pvd')
conservedpvd = PVDReader(FileName='/import/freenas-m-05-seissol/dli/ExaHyPE/Obv/test2/obv-dam.pvd')
conservedpvd.CellArrays = ['time', 'Q']

# find source
# conservedpvd = FindSource('obv-nodam.pvd')

# create a new 'Plot Over Line'
plotOverLine1 = PlotOverLine(Input=conservedpvd,
    Source='High Resolution Line Source')

# init the 'High Resolution Line Source' selected for 'Source'
plotOverLine1.Source.Point1 = [-9746.83984375, -6746.83984375, 0.0]
plotOverLine1.Source.Point2 = [10000.0, 7177.22021484375, 0.0]

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1482, 799]

# show data in view
plotOverLine1Display = Show(plotOverLine1, renderView1)

# trace defaults for the display properties.
plotOverLine1Display.Representation = 'Surface'
plotOverLine1Display.ColorArrayName = [None, '']
plotOverLine1Display.OSPRayScaleArray = 'PreviousRefinementStatus'
plotOverLine1Display.OSPRayScaleFunction = 'PiecewiseFunction'
plotOverLine1Display.SelectOrientationVectors = 'None'
plotOverLine1Display.ScaleFactor = 1974.6839843750001
plotOverLine1Display.SelectScaleArray = 'None'
plotOverLine1Display.GlyphType = 'Arrow'
plotOverLine1Display.GlyphTableIndexArray = 'None'
plotOverLine1Display.GaussianRadius = 98.73419921875
plotOverLine1Display.SetScaleArray = ['POINTS', 'PreviousRefinementStatus']
plotOverLine1Display.ScaleTransferFunction = 'PiecewiseFunction'
plotOverLine1Display.OpacityArray = ['POINTS', 'PreviousRefinementStatus']
plotOverLine1Display.OpacityTransferFunction = 'PiecewiseFunction'
plotOverLine1Display.DataAxesGrid = 'GridAxesRepresentation'
plotOverLine1Display.SelectionCellLabelFontFile = ''
plotOverLine1Display.SelectionPointLabelFontFile = ''
plotOverLine1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
plotOverLine1Display.DataAxesGrid.XTitleFontFile = ''
plotOverLine1Display.DataAxesGrid.YTitleFontFile = ''
plotOverLine1Display.DataAxesGrid.ZTitleFontFile = ''
plotOverLine1Display.DataAxesGrid.XLabelFontFile = ''
plotOverLine1Display.DataAxesGrid.YLabelFontFile = ''
plotOverLine1Display.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
plotOverLine1Display.PolarAxes.PolarAxisTitleFontFile = ''
plotOverLine1Display.PolarAxes.PolarAxisLabelFontFile = ''
plotOverLine1Display.PolarAxes.LastRadialAxisTextFontFile = ''
plotOverLine1Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# Create a new 'Line Chart View'
lineChartView1 = CreateView('XYChartView')
lineChartView1.ViewSize = [736, 799]
lineChartView1.ChartTitleFontFile = ''
lineChartView1.LeftAxisTitleFontFile = ''
lineChartView1.LeftAxisRangeMaximum = 6.66
lineChartView1.LeftAxisLabelFontFile = ''
lineChartView1.BottomAxisTitleFontFile = ''
lineChartView1.BottomAxisRangeMaximum = 6.66
lineChartView1.BottomAxisLabelFontFile = ''
lineChartView1.RightAxisRangeMaximum = 6.66
lineChartView1.RightAxisLabelFontFile = ''
lineChartView1.TopAxisTitleFontFile = ''
lineChartView1.TopAxisRangeMaximum = 6.66
lineChartView1.TopAxisLabelFontFile = ''

# get layout
layout1 = GetLayout()

# place view in the layout
layout1.AssignView(2, lineChartView1)

# show data in view
plotOverLine1Display_1 = Show(plotOverLine1, lineChartView1)

# trace defaults for the display properties.
plotOverLine1Display_1.CompositeDataSetIndex = [0]
plotOverLine1Display_1.UseIndexForXAxis = 0
plotOverLine1Display_1.XArrayName = 'time'
plotOverLine1Display_1.SeriesVisibility = ['PreviousRefinementStatus', 'Q_Magnitude', 'RefinementStatus', 'time']
plotOverLine1Display_1.SeriesLabel = ['arc_length', 'arc_length', 'PreviousRefinementStatus', 'PreviousRefinementStatus', 'Q_0', 'Q_0', 'Q_1', 'Q_1', 'Q_2', 'Q_2', 'Q_3', 'Q_3', 'Q_4', 'Q_4', 'Q_5', 'Q_5', 'Q_6', 'Q_6', 'Q_7', 'Q_7', 'Q_8', 'Q_8', 'Q_9', 'Q_9', 'Q_10', 'Q_10', 'Q_11', 'Q_11', 'Q_12', 'Q_12', 'Q_13', 'Q_13', 'Q_14', 'Q_14', 'Q_15', 'Q_15', 'Q_16', 'Q_16', 'Q_17', 'Q_17', 'Q_18', 'Q_18', 'Q_19', 'Q_19', 'Q_20', 'Q_20', 'Q_21', 'Q_21', 'Q_22', 'Q_22', 'Q_23', 'Q_23', 'Q_24', 'Q_24', 'Q_25', 'Q_25', 'Q_26', 'Q_26', 'Q_27', 'Q_27', 'Q_28', 'Q_28', 'Q_29', 'Q_29', 'Q_30', 'Q_30', 'Q_31', 'Q_31', 'Q_32', 'Q_32', 'Q_33', 'Q_33', 'Q_34', 'Q_34', 'Q_35', 'Q_35', 'Q_36', 'Q_36', 'Q_37', 'Q_37', 'Q_38', 'Q_38', 'Q_39', 'Q_39', 'Q_40', 'Q_40', 'Q_41', 'Q_41', 'Q_42', 'Q_42', 'Q_Magnitude', 'Q_Magnitude', 'RefinementStatus', 'RefinementStatus', 'time', 'time', 'vtkValidPointMask', 'vtkValidPointMask', 'Points_X', 'Points_X', 'Points_Y', 'Points_Y', 'Points_Z', 'Points_Z', 'Points_Magnitude', 'Points_Magnitude']
plotOverLine1Display_1.SeriesColor = ['arc_length', '0', '0', '0', 'PreviousRefinementStatus', '0.89', '0.1', '0.11', 'Q_0', '0.22', '0.49', '0.72', 'Q_1', '0.3', '0.69', '0.29', 'Q_2', '0.6', '0.31', '0.64', 'Q_3', '1', '0.5', '0', 'Q_4', '0.65', '0.34', '0.16', 'Q_5', '0', '0', '0', 'Q_6', '0.89', '0.1', '0.11', 'Q_7', '0.22', '0.49', '0.72', 'Q_8', '0.3', '0.69', '0.29', 'Q_9', '0.6', '0.31', '0.64', 'Q_10', '1', '0.5', '0', 'Q_11', '0.65', '0.34', '0.16', 'Q_12', '0', '0', '0', 'Q_13', '0.89', '0.1', '0.11', 'Q_14', '0.22', '0.49', '0.72', 'Q_15', '0.3', '0.69', '0.29', 'Q_16', '0.6', '0.31', '0.64', 'Q_17', '1', '0.5', '0', 'Q_18', '0.65', '0.34', '0.16', 'Q_19', '0', '0', '0', 'Q_20', '0.89', '0.1', '0.11', 'Q_21', '0.22', '0.49', '0.72', 'Q_22', '0.3', '0.69', '0.29', 'Q_23', '0.6', '0.31', '0.64', 'Q_24', '1', '0.5', '0', 'Q_25', '0.65', '0.34', '0.16', 'Q_26', '0', '0', '0', 'Q_27', '0.89', '0.1', '0.11', 'Q_28', '0.22', '0.49', '0.72', 'Q_29', '0.3', '0.69', '0.29', 'Q_30', '0.6', '0.31', '0.64', 'Q_31', '1', '0.5', '0', 'Q_32', '0.65', '0.34', '0.16', 'Q_33', '0', '0', '0', 'Q_34', '0.89', '0.1', '0.11', 'Q_35', '0.22', '0.49', '0.72', 'Q_36', '0.3', '0.69', '0.29', 'Q_37', '0.6', '0.31', '0.64', 'Q_38', '1', '0.5', '0', 'Q_39', '0.65', '0.34', '0.16', 'Q_40', '0', '0', '0', 'Q_41', '0.89', '0.1', '0.11', 'Q_42', '0.22', '0.49', '0.72', 'Q_Magnitude', '0.3', '0.69', '0.29', 'RefinementStatus', '0.6', '0.31', '0.64', 'time', '1', '0.5', '0', 'vtkValidPointMask', '0.65', '0.34', '0.16', 'Points_X', '0', '0', '0', 'Points_Y', '0.89', '0.1', '0.11', 'Points_Z', '0.22', '0.49', '0.72', 'Points_Magnitude', '0.3', '0.69', '0.29']
plotOverLine1Display_1.SeriesPlotCorner = ['arc_length', '0', 'PreviousRefinementStatus', '0', 'Q_0', '0', 'Q_1', '0', 'Q_2', '0', 'Q_3', '0', 'Q_4', '0', 'Q_5', '0', 'Q_6', '0', 'Q_7', '0', 'Q_8', '0', 'Q_9', '0', 'Q_10', '0', 'Q_11', '0', 'Q_12', '0', 'Q_13', '0', 'Q_14', '0', 'Q_15', '0', 'Q_16', '0', 'Q_17', '0', 'Q_18', '0', 'Q_19', '0', 'Q_20', '0', 'Q_21', '0', 'Q_22', '0', 'Q_23', '0', 'Q_24', '0', 'Q_25', '0', 'Q_26', '0', 'Q_27', '0', 'Q_28', '0', 'Q_29', '0', 'Q_30', '0', 'Q_31', '0', 'Q_32', '0', 'Q_33', '0', 'Q_34', '0', 'Q_35', '0', 'Q_36', '0', 'Q_37', '0', 'Q_38', '0', 'Q_39', '0', 'Q_40', '0', 'Q_41', '0', 'Q_42', '0', 'Q_Magnitude', '0', 'RefinementStatus', '0', 'time', '0', 'vtkValidPointMask', '0', 'Points_X', '0', 'Points_Y', '0', 'Points_Z', '0', 'Points_Magnitude', '0']
plotOverLine1Display_1.SeriesLabelPrefix = ''
plotOverLine1Display_1.SeriesLineStyle = ['arc_length', '1', 'PreviousRefinementStatus', '1', 'Q_0', '1', 'Q_1', '1', 'Q_2', '1', 'Q_3', '1', 'Q_4', '1', 'Q_5', '1', 'Q_6', '1', 'Q_7', '1', 'Q_8', '1', 'Q_9', '1', 'Q_10', '1', 'Q_11', '1', 'Q_12', '1', 'Q_13', '1', 'Q_14', '1', 'Q_15', '1', 'Q_16', '1', 'Q_17', '1', 'Q_18', '1', 'Q_19', '1', 'Q_20', '1', 'Q_21', '1', 'Q_22', '1', 'Q_23', '1', 'Q_24', '1', 'Q_25', '1', 'Q_26', '1', 'Q_27', '1', 'Q_28', '1', 'Q_29', '1', 'Q_30', '1', 'Q_31', '1', 'Q_32', '1', 'Q_33', '1', 'Q_34', '1', 'Q_35', '1', 'Q_36', '1', 'Q_37', '1', 'Q_38', '1', 'Q_39', '1', 'Q_40', '1', 'Q_41', '1', 'Q_42', '1', 'Q_Magnitude', '1', 'RefinementStatus', '1', 'time', '1', 'vtkValidPointMask', '1', 'Points_X', '1', 'Points_Y', '1', 'Points_Z', '1', 'Points_Magnitude', '1']
plotOverLine1Display_1.SeriesLineThickness = ['arc_length', '2', 'PreviousRefinementStatus', '2', 'Q_0', '2', 'Q_1', '2', 'Q_2', '2', 'Q_3', '2', 'Q_4', '2', 'Q_5', '2', 'Q_6', '2', 'Q_7', '2', 'Q_8', '2', 'Q_9', '2', 'Q_10', '2', 'Q_11', '2', 'Q_12', '2', 'Q_13', '2', 'Q_14', '2', 'Q_15', '2', 'Q_16', '2', 'Q_17', '2', 'Q_18', '2', 'Q_19', '2', 'Q_20', '2', 'Q_21', '2', 'Q_22', '2', 'Q_23', '2', 'Q_24', '2', 'Q_25', '2', 'Q_26', '2', 'Q_27', '2', 'Q_28', '2', 'Q_29', '2', 'Q_30', '2', 'Q_31', '2', 'Q_32', '2', 'Q_33', '2', 'Q_34', '2', 'Q_35', '2', 'Q_36', '2', 'Q_37', '2', 'Q_38', '2', 'Q_39', '2', 'Q_40', '2', 'Q_41', '2', 'Q_42', '2', 'Q_Magnitude', '2', 'RefinementStatus', '2', 'time', '2', 'vtkValidPointMask', '2', 'Points_X', '2', 'Points_Y', '2', 'Points_Z', '2', 'Points_Magnitude', '2']
plotOverLine1Display_1.SeriesMarkerStyle = ['arc_length', '0', 'PreviousRefinementStatus', '0', 'Q_0', '0', 'Q_1', '0', 'Q_2', '0', 'Q_3', '0', 'Q_4', '0', 'Q_5', '0', 'Q_6', '0', 'Q_7', '0', 'Q_8', '0', 'Q_9', '0', 'Q_10', '0', 'Q_11', '0', 'Q_12', '0', 'Q_13', '0', 'Q_14', '0', 'Q_15', '0', 'Q_16', '0', 'Q_17', '0', 'Q_18', '0', 'Q_19', '0', 'Q_20', '0', 'Q_21', '0', 'Q_22', '0', 'Q_23', '0', 'Q_24', '0', 'Q_25', '0', 'Q_26', '0', 'Q_27', '0', 'Q_28', '0', 'Q_29', '0', 'Q_30', '0', 'Q_31', '0', 'Q_32', '0', 'Q_33', '0', 'Q_34', '0', 'Q_35', '0', 'Q_36', '0', 'Q_37', '0', 'Q_38', '0', 'Q_39', '0', 'Q_40', '0', 'Q_41', '0', 'Q_42', '0', 'Q_Magnitude', '0', 'RefinementStatus', '0', 'time', '0', 'vtkValidPointMask', '0', 'Points_X', '0', 'Points_Y', '0', 'Points_Z', '0', 'Points_Magnitude', '0']

# update the view to ensure updated data information
renderView1.Update()

# update the view to ensure updated data information
lineChartView1.Update()

# Properties modified on plotOverLine1.Source
plotOverLine1.Source.Point1 = [-9746.83984375, 215.190185546875, 0.0]
plotOverLine1.Source.Point2 = [10000.0, 215.190185546875, 0.0]

# update the view to ensure updated data information
renderView1.Update()

# update the view to ensure updated data information
lineChartView1.Update()

# Properties modified on plotOverLine1.Source
plotOverLine1.Source.Point2 = [0.0, 215.190185546875, 0.0]

# update the view to ensure updated data information
renderView1.Update()

# update the view to ensure updated data information
lineChartView1.Update()

# Properties modified on plotOverLine1.Source
plotOverLine1.Source.Point1 = [-9746.83984375, 100.0, 0.0]

# update the view to ensure updated data information
renderView1.Update()

# update the view to ensure updated data information
lineChartView1.Update()

# Properties modified on plotOverLine1.Source
plotOverLine1.Source.Point2 = [0.0, 100.0, 0.0]

# update the view to ensure updated data information
renderView1.Update()

# update the view to ensure updated data information
lineChartView1.Update()

# Properties modified on plotOverLine1.Source
plotOverLine1.Source.Resolution = 3

# update the view to ensure updated data information
renderView1.Update()

# update the view to ensure updated data information
lineChartView1.Update()

# Properties modified on plotOverLine1.Source
plotOverLine1.Source.Resolution = 30

# update the view to ensure updated data information
renderView1.Update()

# update the view to ensure updated data information
lineChartView1.Update()

# Properties modified on plotOverLine1.Source
plotOverLine1.Source.Resolution = 300

# update the view to ensure updated data information
renderView1.Update()

# update the view to ensure updated data information
lineChartView1.Update()

# Properties modified on plotOverLine1Display_1
plotOverLine1Display_1.SeriesVisibility = []
plotOverLine1Display_1.SeriesColor = ['arc_length', '0', '0', '0', 'PreviousRefinementStatus', '0.889998', '0.100008', '0.110002', 'Q_0', '0.220005', '0.489998', '0.719997', 'Q_1', '0.300008', '0.689998', '0.289998', 'Q_2', '0.6', '0.310002', '0.639994', 'Q_3', '1', '0.500008', '0', 'Q_4', '0.650004', '0.340002', '0.160006', 'Q_5', '0', '0', '0', 'Q_6', '0.889998', '0.100008', '0.110002', 'Q_7', '0.220005', '0.489998', '0.719997', 'Q_8', '0.300008', '0.689998', '0.289998', 'Q_9', '0.6', '0.310002', '0.639994', 'Q_10', '1', '0.500008', '0', 'Q_11', '0.650004', '0.340002', '0.160006', 'Q_12', '0', '0', '0', 'Q_13', '0.889998', '0.100008', '0.110002', 'Q_14', '0.220005', '0.489998', '0.719997', 'Q_15', '0.300008', '0.689998', '0.289998', 'Q_16', '0.6', '0.310002', '0.639994', 'Q_17', '1', '0.500008', '0', 'Q_18', '0.650004', '0.340002', '0.160006', 'Q_19', '0', '0', '0', 'Q_20', '0.889998', '0.100008', '0.110002', 'Q_21', '0.220005', '0.489998', '0.719997', 'Q_22', '0.300008', '0.689998', '0.289998', 'Q_23', '0.6', '0.310002', '0.639994', 'Q_24', '1', '0.500008', '0', 'Q_25', '0.650004', '0.340002', '0.160006', 'Q_26', '0', '0', '0', 'Q_27', '0.889998', '0.100008', '0.110002', 'Q_28', '0.220005', '0.489998', '0.719997', 'Q_29', '0.300008', '0.689998', '0.289998', 'Q_30', '0.6', '0.310002', '0.639994', 'Q_31', '1', '0.500008', '0', 'Q_32', '0.650004', '0.340002', '0.160006', 'Q_33', '0', '0', '0', 'Q_34', '0.889998', '0.100008', '0.110002', 'Q_35', '0.220005', '0.489998', '0.719997', 'Q_36', '0.300008', '0.689998', '0.289998', 'Q_37', '0.6', '0.310002', '0.639994', 'Q_38', '1', '0.500008', '0', 'Q_39', '0.650004', '0.340002', '0.160006', 'Q_40', '0', '0', '0', 'Q_41', '0.889998', '0.100008', '0.110002', 'Q_42', '0.220005', '0.489998', '0.719997', 'Q_Magnitude', '0.300008', '0.689998', '0.289998', 'RefinementStatus', '0.6', '0.310002', '0.639994', 'time', '1', '0.500008', '0', 'vtkValidPointMask', '0.650004', '0.340002', '0.160006', 'Points_X', '0', '0', '0', 'Points_Y', '0.889998', '0.100008', '0.110002', 'Points_Z', '0.220005', '0.489998', '0.719997', 'Points_Magnitude', '0.300008', '0.689998', '0.289998']
plotOverLine1Display_1.SeriesPlotCorner = ['Points_Magnitude', '0', 'Points_X', '0', 'Points_Y', '0', 'Points_Z', '0', 'PreviousRefinementStatus', '0', 'Q_0', '0', 'Q_1', '0', 'Q_10', '0', 'Q_11', '0', 'Q_12', '0', 'Q_13', '0', 'Q_14', '0', 'Q_15', '0', 'Q_16', '0', 'Q_17', '0', 'Q_18', '0', 'Q_19', '0', 'Q_2', '0', 'Q_20', '0', 'Q_21', '0', 'Q_22', '0', 'Q_23', '0', 'Q_24', '0', 'Q_25', '0', 'Q_26', '0', 'Q_27', '0', 'Q_28', '0', 'Q_29', '0', 'Q_3', '0', 'Q_30', '0', 'Q_31', '0', 'Q_32', '0', 'Q_33', '0', 'Q_34', '0', 'Q_35', '0', 'Q_36', '0', 'Q_37', '0', 'Q_38', '0', 'Q_39', '0', 'Q_4', '0', 'Q_40', '0', 'Q_41', '0', 'Q_42', '0', 'Q_5', '0', 'Q_6', '0', 'Q_7', '0', 'Q_8', '0', 'Q_9', '0', 'Q_Magnitude', '0', 'RefinementStatus', '0', 'arc_length', '0', 'time', '0', 'vtkValidPointMask', '0']
plotOverLine1Display_1.SeriesLineStyle = ['Points_Magnitude', '1', 'Points_X', '1', 'Points_Y', '1', 'Points_Z', '1', 'PreviousRefinementStatus', '1', 'Q_0', '1', 'Q_1', '1', 'Q_10', '1', 'Q_11', '1', 'Q_12', '1', 'Q_13', '1', 'Q_14', '1', 'Q_15', '1', 'Q_16', '1', 'Q_17', '1', 'Q_18', '1', 'Q_19', '1', 'Q_2', '1', 'Q_20', '1', 'Q_21', '1', 'Q_22', '1', 'Q_23', '1', 'Q_24', '1', 'Q_25', '1', 'Q_26', '1', 'Q_27', '1', 'Q_28', '1', 'Q_29', '1', 'Q_3', '1', 'Q_30', '1', 'Q_31', '1', 'Q_32', '1', 'Q_33', '1', 'Q_34', '1', 'Q_35', '1', 'Q_36', '1', 'Q_37', '1', 'Q_38', '1', 'Q_39', '1', 'Q_4', '1', 'Q_40', '1', 'Q_41', '1', 'Q_42', '1', 'Q_5', '1', 'Q_6', '1', 'Q_7', '1', 'Q_8', '1', 'Q_9', '1', 'Q_Magnitude', '1', 'RefinementStatus', '1', 'arc_length', '1', 'time', '1', 'vtkValidPointMask', '1']
plotOverLine1Display_1.SeriesLineThickness = ['Points_Magnitude', '2', 'Points_X', '2', 'Points_Y', '2', 'Points_Z', '2', 'PreviousRefinementStatus', '2', 'Q_0', '2', 'Q_1', '2', 'Q_10', '2', 'Q_11', '2', 'Q_12', '2', 'Q_13', '2', 'Q_14', '2', 'Q_15', '2', 'Q_16', '2', 'Q_17', '2', 'Q_18', '2', 'Q_19', '2', 'Q_2', '2', 'Q_20', '2', 'Q_21', '2', 'Q_22', '2', 'Q_23', '2', 'Q_24', '2', 'Q_25', '2', 'Q_26', '2', 'Q_27', '2', 'Q_28', '2', 'Q_29', '2', 'Q_3', '2', 'Q_30', '2', 'Q_31', '2', 'Q_32', '2', 'Q_33', '2', 'Q_34', '2', 'Q_35', '2', 'Q_36', '2', 'Q_37', '2', 'Q_38', '2', 'Q_39', '2', 'Q_4', '2', 'Q_40', '2', 'Q_41', '2', 'Q_42', '2', 'Q_5', '2', 'Q_6', '2', 'Q_7', '2', 'Q_8', '2', 'Q_9', '2', 'Q_Magnitude', '2', 'RefinementStatus', '2', 'arc_length', '2', 'time', '2', 'vtkValidPointMask', '2']
plotOverLine1Display_1.SeriesMarkerStyle = ['Points_Magnitude', '0', 'Points_X', '0', 'Points_Y', '0', 'Points_Z', '0', 'PreviousRefinementStatus', '0', 'Q_0', '0', 'Q_1', '0', 'Q_10', '0', 'Q_11', '0', 'Q_12', '0', 'Q_13', '0', 'Q_14', '0', 'Q_15', '0', 'Q_16', '0', 'Q_17', '0', 'Q_18', '0', 'Q_19', '0', 'Q_2', '0', 'Q_20', '0', 'Q_21', '0', 'Q_22', '0', 'Q_23', '0', 'Q_24', '0', 'Q_25', '0', 'Q_26', '0', 'Q_27', '0', 'Q_28', '0', 'Q_29', '0', 'Q_3', '0', 'Q_30', '0', 'Q_31', '0', 'Q_32', '0', 'Q_33', '0', 'Q_34', '0', 'Q_35', '0', 'Q_36', '0', 'Q_37', '0', 'Q_38', '0', 'Q_39', '0', 'Q_4', '0', 'Q_40', '0', 'Q_41', '0', 'Q_42', '0', 'Q_5', '0', 'Q_6', '0', 'Q_7', '0', 'Q_8', '0', 'Q_9', '0', 'Q_Magnitude', '0', 'RefinementStatus', '0', 'arc_length', '0', 'time', '0', 'vtkValidPointMask', '0']

# Properties modified on plotOverLine1Display_1
plotOverLine1Display_1.SeriesVisibility = ['Q_1']

# Properties modified on plotOverLine1Display_1
plotOverLine1Display_1.XArrayName = 'Points_X'

# save data
SaveData('/import/freenas-m-05-seissol/dli/ExaHyPE/Obv/test2/spatiotemporal_data.csv', proxy=plotOverLine1, WriteTimeSteps=1)

# Properties modified on plotOverLine1.Source
plotOverLine1.Source.Point1 = [-9746.83984375, -100.0, 0.0]

# update the view to ensure updated data information
renderView1.Update()

# update the view to ensure updated data information
lineChartView1.Update()

# Properties modified on plotOverLine1.Source
plotOverLine1.Source.Point2 = [0.0, -100.0, 0.0]

# update the view to ensure updated data information
renderView1.Update()

# update the view to ensure updated data information
lineChartView1.Update()

# save data
SaveData('/import/freenas-m-05-seissol/dli/ExaHyPE/Obv/test2/spatiotemporal_data2.csv', proxy=plotOverLine1, WriteTimeSteps=1)

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [126.580078125, 215.190185546875, 10000.0]
renderView1.CameraFocalPoint = [126.580078125, 215.190185546875, 0.0]
renderView1.CameraParallelScale = 12081.154045972164

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).
