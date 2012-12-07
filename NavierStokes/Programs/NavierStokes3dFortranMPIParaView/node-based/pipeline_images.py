
try: paraview.simple
except: from paraview.simple import *

def RequestDataDescription(datadescription):
    "Callback to populate the request for current timestep"
    if datadescription.GetForceOutput() == True:
        for i in range(datadescription.GetNumberOfInputDescriptions()):
            datadescription.GetInputDescription(i).AllFieldsOn()
            datadescription.GetInputDescription(i).GenerateMeshOn()
        return

    timestep = datadescription.GetTimeStep()

    input_name = 'input'
    if (timestep % 1 == 0) :
        datadescription.GetInputDescriptionByName(input_name).AllFieldsOn()
        datadescription.GetInputDescriptionByName(input_name).GenerateMeshOn()
    else:
        datadescription.GetInputDescriptionByName(input_name).AllFieldsOff()
        datadescription.GetInputDescriptionByName(input_name).GenerateMeshOff()


def DoCoProcessing(datadescription):
    "Callback to do co-processing for current timestep"
    cp_writers = []
    cp_views = []
    timestep = datadescription.GetTimeStep()

    RenderView1 = CreateView( CreateRenderView, "np-2_%t.png", 1, 0, 1, cp_views )
    RenderView1.LightSpecularColor = [1.0, 1.0, 1.0]
    RenderView1.KeyLightAzimuth = 10.0
    RenderView1.UseTexturedBackground = 0
    RenderView1.UseLight = 1
    RenderView1.CameraPosition = [217.11071921506289, 369.10961277586267, 116.87368903495434]
    RenderView1.FillLightKFRatio = 3.0
    RenderView1.Background2 = [0.0, 0.0, 0.16470588235294117]
    RenderView1.FillLightAzimuth = -10.0
    RenderView1.LODResolution = 50.0
    RenderView1.BackgroundTexture = []
    RenderView1.InteractionMode = '3D'
    RenderView1.StencilCapable = 1
    RenderView1.LightIntensity = 1.0
    RenderView1.CameraFocalPoint = [23.825137927181057, 0.0099511400572356041, 33.274249608756051]
    RenderView1.ImageReductionFactor = 2
    RenderView1.CameraViewAngle = 30.0
    RenderView1.CameraParallelScale = 109.98522628062371
    RenderView1.EyeAngle = 2.0
    RenderView1.HeadLightKHRatio = 3.0
    RenderView1.StereoRender = 0
    RenderView1.KeyLightIntensity = 0.75
    RenderView1.BackLightAzimuth = 110.0
    RenderView1.OrientationAxesInteractivity = 0
    RenderView1.UseInteractiveRenderingForSceenshots = 0
    RenderView1.UseOffscreenRendering = 0
    RenderView1.Background = [0.0, 0.0, 0.0]
    RenderView1.UseOffscreenRenderingForScreenshots = 0
    RenderView1.NonInteractiveRenderDelay = 2
    RenderView1.CenterOfRotation = [63.5, 63.5, 63.5]
    RenderView1.CameraParallelProjection = 0
    RenderView1.CompressorConfig = 'vtkSquirtCompressor 0 3'
    RenderView1.HeadLightWarmth = 0.5
    RenderView1.MaximumNumberOfPeels = 4
    RenderView1.LightDiffuseColor = [1.0, 1.0, 1.0]
    RenderView1.StereoType = 'Red-Blue'
    RenderView1.DepthPeeling = 1
    RenderView1.BackLightKBRatio = 3.5
    RenderView1.StereoCapableWindow = 1
    RenderView1.CameraViewUp = [-0.10216114394387266, -0.16855412123856731, 0.98038390892628069]
    RenderView1.LightType = 'HeadLight'
    RenderView1.LightAmbientColor = [1.0, 1.0, 1.0]
    RenderView1.RemoteRenderThreshold = 3.0
    RenderView1.CacheKey = 0.0
    RenderView1.UseCache = 0
    RenderView1.KeyLightElevation = 50.0
    RenderView1.CenterAxesVisibility = 0
    RenderView1.MaintainLuminance = 0
    RenderView1.StillRenderImageReductionFactor = 1
    RenderView1.BackLightWarmth = 0.5
    RenderView1.FillLightElevation = -75.0
    RenderView1.MultiSamples = 0
    RenderView1.FillLightWarmth = 0.40000000000000002
    RenderView1.AlphaBitPlanes = 1
    RenderView1.LightSwitch = 1
    RenderView1.OrientationAxesVisibility = 1
    RenderView1.CameraClippingRange = [150.26173427523366, 592.80491988083509]
    RenderView1.BackLightElevation = 0.0
    RenderView1.ViewTime = 0.0
    RenderView1.OrientationAxesOutlineColor = [1.0, 1.0, 1.0]
    RenderView1.LODThreshold = 18.199999999999999
    RenderView1.CollectGeometryThreshold = 100.0
    RenderView1.UseGradientBackground = 0
    RenderView1.KeyLightWarmth = 0.59999999999999998
    RenderView1.OrientationAxesLabelColor = [1.0, 1.0, 1.0]
    
    test2_0_0_vti = CreateProducer( datadescription, "input" )
    
    Calculator1 = Calculator( guiName="Calculator1", Function='mag(realtempx*iHat+realtempy*jHat+realtempz*kHat)', ReplacementValue=0.0, ResultArrayName='Result', ReplaceInvalidResults=1, AttributeMode='point_data', CoordinateResults=0 )
    
    a1_realtempx_PiecewiseFunction = CreatePiecewiseFunction( Points=[0.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0] )
    
    a1_realtempy_PiecewiseFunction = CreatePiecewiseFunction( Points=[0.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0] )
    
    a1_realtempz_PiecewiseFunction = CreatePiecewiseFunction( Points=[0.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0] )
    
    a1_Result_PiecewiseFunction = CreatePiecewiseFunction( Points=[5.1558018820268853e-15, 0.58799999952316284, 0.5, 0.60829490423202515, 0.76975244283676147, 0.0, 0.5, 0.0, 0.86602500000000004, 0.0, 0.5, 0.0, 0.95262794416286289, 0.0, 0.5, 0.39072847366333008, 1.7320507764816284, 0.71999996900558472, 0.5, 0.0] )
    
    a1_realtempx_PVLookupTable = GetLookupTableForArray( "realtempx", 1, Discretize=1, RGBPoints=[-1.5000000000000338, 0.23000000000000001, 0.29899999999999999, 0.754, 1.5000000000000338, 0.70599999999999996, 0.016, 0.14999999999999999], UseLogScale=0, VectorComponent=0, NanColor=[0.25, 0.0, 0.0], NumberOfTableValues=256, ColorSpace='Diverging', VectorMode='Magnitude', HSVWrap=0, ScalarRangeInitialized=1.0, LockScalarRange=0 )
    
    a1_realtempy_PVLookupTable = GetLookupTableForArray( "realtempy", 1, Discretize=1, RGBPoints=[-1.5000000000000393, 0.23000000000000001, 0.29899999999999999, 0.754, 1.5000000000000409, 0.70599999999999996, 0.016, 0.14999999999999999], UseLogScale=0, VectorComponent=0, NanColor=[0.25, 0.0, 0.0], NumberOfTableValues=256, ColorSpace='Diverging', VectorMode='Magnitude', HSVWrap=0, ScalarRangeInitialized=1.0, LockScalarRange=0 )
    
    a1_realtempz_PVLookupTable = GetLookupTableForArray( "realtempz", 1, Discretize=1, RGBPoints=[-1.7320508075688363, 0.23000000000000001, 0.29899999999999999, 0.754, 1.7320508075688372, 0.70599999999999996, 0.016, 0.14999999999999999], UseLogScale=0, VectorComponent=0, NanColor=[0.25, 0.0, 0.0], NumberOfTableValues=256, ColorSpace='Diverging', VectorMode='Magnitude', HSVWrap=0, ScalarRangeInitialized=1.0, LockScalarRange=0 )
    
    a1_Result_PVLookupTable = GetLookupTableForArray( "Result", 1, Discretize=1, RGBPoints=[5.1558017812958882e-15, 0.0, 1.0, 1.0, 0.77942286340597955, 0.015686274509803921, 0.0, 0.92156862745098034, 0.86602540378442117, 0.22745098039215686, 0.0, 0.38823529411764707, 0.95262794416286289, 0.88627450980392153, 0.0, 0.058823529411764705, 1.7320508075688372, 1.0, 1.0, 0.0], UseLogScale=0, VectorComponent=0, NanColor=[1.0, 1.0, 0.0], NumberOfTableValues=256, ColorSpace='RGB', VectorMode='Magnitude', HSVWrap=0, ScalarRangeInitialized=1.0, LockScalarRange=0 )
    
    SetActiveSource(test2_0_0_vti)
    DataRepresentation1 = Show()
    DataRepresentation1.CubeAxesZAxisVisibility = 1
    DataRepresentation1.SelectionPointLabelColor = [0.5, 0.5, 0.5]
    DataRepresentation1.SelectionPointFieldDataArrayName = 'realtempx'
    DataRepresentation1.SuppressLOD = 0
    DataRepresentation1.CubeAxesXGridLines = 0
    DataRepresentation1.CubeAxesYAxisTickVisibility = 1
    DataRepresentation1.CubeAxesColor = [1.0, 1.0, 1.0]
    DataRepresentation1.Position = [0.0, 0.0, 0.0]
    DataRepresentation1.BackfaceRepresentation = 'Follow Frontface'
    DataRepresentation1.SelectionOpacity = 1.0
    DataRepresentation1.SelectionPointLabelShadow = 0
    DataRepresentation1.CubeAxesYGridLines = 0
    DataRepresentation1.CubeAxesZAxisRange = [0.0, 1.0]
    DataRepresentation1.OrientationMode = 'Direction'
    DataRepresentation1.Source.TipResolution = 6
    DataRepresentation1.ScaleMode = 'No Data Scaling Off'
    DataRepresentation1.Diffuse = 1.0
    DataRepresentation1.SelectionUseOutline = 0
    DataRepresentation1.CubeAxesZTitle = 'Z-Axis'
    DataRepresentation1.Specular = 0.10000000000000001
    DataRepresentation1.SelectionVisibility = 1
    DataRepresentation1.InterpolateScalarsBeforeMapping = 1
    DataRepresentation1.CubeAxesZAxisTickVisibility = 1
    DataRepresentation1.Origin = [0.0, 0.0, 0.0]
    DataRepresentation1.CubeAxesVisibility = 0
    DataRepresentation1.Scale = [1.0, 1.0, 1.0]
    DataRepresentation1.SelectionCellLabelJustification = 'Left'
    DataRepresentation1.DiffuseColor = [1.0, 1.0, 1.0]
    DataRepresentation1.Shade = 0
    DataRepresentation1.SelectionCellLabelOpacity = 1.0
    DataRepresentation1.CubeAxesInertia = 1
    DataRepresentation1.Source = "Arrow"
    DataRepresentation1.Source.Invert = 0
    DataRepresentation1.Masking = 0
    DataRepresentation1.Opacity = 1.0
    DataRepresentation1.LineWidth = 1.0
    DataRepresentation1.MeshVisibility = 0
    DataRepresentation1.Visibility = 0
    DataRepresentation1.SelectionCellLabelFontSize = 18
    DataRepresentation1.CubeAxesCornerOffset = 0.0
    DataRepresentation1.SelectionPointLabelJustification = 'Left'
    DataRepresentation1.SelectionPointLabelVisibility = 0
    DataRepresentation1.SelectOrientationVectors = ''
    DataRepresentation1.CubeAxesTickLocation = 'Inside'
    DataRepresentation1.BackfaceDiffuseColor = [1.0, 1.0, 1.0]
    DataRepresentation1.CubeAxesYAxisVisibility = 1
    DataRepresentation1.SelectionPointLabelFontFamily = 'Arial'
    DataRepresentation1.Source.ShaftResolution = 6
    DataRepresentation1.CubeAxesFlyMode = 'Closest Triad'
    DataRepresentation1.SelectScaleArray = ''
    DataRepresentation1.CubeAxesYTitle = 'Y-Axis'
    DataRepresentation1.ColorAttributeType = 'POINT_DATA'
    DataRepresentation1.SpecularPower = 100.0
    DataRepresentation1.Texture = []
    DataRepresentation1.SelectionCellLabelShadow = 0
    DataRepresentation1.AmbientColor = [1.0, 1.0, 1.0]
    DataRepresentation1.MapScalars = 1
    DataRepresentation1.PointSize = 2.0
    DataRepresentation1.Source.TipLength = 0.34999999999999998
    DataRepresentation1.SelectionCellLabelFormat = ''
    DataRepresentation1.Scaling = 0
    DataRepresentation1.StaticMode = 0
    DataRepresentation1.SelectionCellLabelColor = [0.0, 1.0, 0.0]
    DataRepresentation1.SliceMode = 'XY Plane'
    DataRepresentation1.Source.TipRadius = 0.10000000000000001
    DataRepresentation1.EdgeColor = [0.0, 0.0, 0.50000762951094835]
    DataRepresentation1.CubeAxesXAxisTickVisibility = 1
    DataRepresentation1.SelectionCellLabelVisibility = 0
    DataRepresentation1.NonlinearSubdivisionLevel = 1
    DataRepresentation1.CubeAxesXAxisRange = [0.0, 1.0]
    DataRepresentation1.Representation = 'Surface'
    DataRepresentation1.CubeAxesYAxisRange = [0.0, 1.0]
    DataRepresentation1.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
    DataRepresentation1.Orientation = [0.0, 0.0, 0.0]
    DataRepresentation1.CubeAxesEnableCustomAxisRange = 0
    DataRepresentation1.CubeAxesXTitle = 'X-Axis'
    DataRepresentation1.ScalarOpacityUnitDistance = 1.7320508075688776
    DataRepresentation1.BackfaceOpacity = 1.0
    DataRepresentation1.SelectionPointLabelFontSize = 18
    DataRepresentation1.SelectionCellFieldDataArrayName = 'vtkOriginalCellIds'
    DataRepresentation1.SelectionColor = [1.0, 0.0, 1.0]
    DataRepresentation1.Ambient = 0.0
    DataRepresentation1.VolumeRenderingMode = 'Smart'
    DataRepresentation1.CubeAxesXAxisMinorTickVisibility = 1
    DataRepresentation1.ScaleFactor = 12.700000000000001
    DataRepresentation1.BackfaceAmbientColor = [1.0, 1.0, 1.0]
    DataRepresentation1.Slice = 0
    DataRepresentation1.Source.ShaftRadius = 0.029999999999999999
    DataRepresentation1.ScalarOpacityFunction = a1_realtempz_PiecewiseFunction
    DataRepresentation1.SelectMaskArray = ''
    DataRepresentation1.SelectionLineWidth = 2.0
    DataRepresentation1.CubeAxesZAxisMinorTickVisibility = 1
    DataRepresentation1.CubeAxesXAxisVisibility = 1
    DataRepresentation1.Interpolation = 'Gouraud'
    DataRepresentation1.SelectionCellLabelFontFamily = 'Arial'
    DataRepresentation1.SelectionCellLabelItalic = 0
    DataRepresentation1.CubeAxesYAxisMinorTickVisibility = 1
    DataRepresentation1.InterpolationType = 'Linear'
    DataRepresentation1.CubeAxesZGridLines = 0
    DataRepresentation1.SelectionPointLabelFormat = ''
    DataRepresentation1.SelectionPointLabelOpacity = 1.0
    DataRepresentation1.Pickable = 1
    DataRepresentation1.CustomBoundsActive = [0, 0, 0]
    DataRepresentation1.SelectionRepresentation = 'Wireframe'
    DataRepresentation1.SelectionPointLabelBold = 0
    DataRepresentation1.ColorArrayName = 'realtempz'
    DataRepresentation1.SelectionPointLabelItalic = 0
    DataRepresentation1.AllowSpecularHighlightingWithScalarColoring = 0
    DataRepresentation1.SpecularColor = [1.0, 1.0, 1.0]
    DataRepresentation1.LookupTable = a1_realtempz_PVLookupTable
    DataRepresentation1.SelectionPointSize = 5.0
    DataRepresentation1.SelectionCellLabelBold = 0
    DataRepresentation1.Orient = 0
    
    SetActiveSource(Calculator1)
    DataRepresentation2 = Show()
    DataRepresentation2.CubeAxesZAxisVisibility = 1
    DataRepresentation2.SelectionPointLabelColor = [0.5, 0.5, 0.5]
    DataRepresentation2.SelectionPointFieldDataArrayName = 'Result'
    DataRepresentation2.SuppressLOD = 0
    DataRepresentation2.CubeAxesXGridLines = 0
    DataRepresentation2.CubeAxesYAxisTickVisibility = 1
    DataRepresentation2.CubeAxesColor = [1.0, 1.0, 1.0]
    DataRepresentation2.Position = [0.0, 0.0, 0.0]
    DataRepresentation2.BackfaceRepresentation = 'Follow Frontface'
    DataRepresentation2.SelectionOpacity = 1.0
    DataRepresentation2.SelectionPointLabelShadow = 0
    DataRepresentation2.CubeAxesYGridLines = 0
    DataRepresentation2.CubeAxesZAxisRange = [0.0, 1.0]
    DataRepresentation2.OrientationMode = 'Direction'
    DataRepresentation2.Source.TipResolution = 6
    DataRepresentation2.ScaleMode = 'No Data Scaling Off'
    DataRepresentation2.Diffuse = 1.0
    DataRepresentation2.SelectionUseOutline = 0
    DataRepresentation2.CubeAxesZTitle = 'Z-Axis'
    DataRepresentation2.Specular = 0.10000000000000001
    DataRepresentation2.SelectionVisibility = 1
    DataRepresentation2.InterpolateScalarsBeforeMapping = 1
    DataRepresentation2.CubeAxesZAxisTickVisibility = 1
    DataRepresentation2.Origin = [0.0, 0.0, 0.0]
    DataRepresentation2.CubeAxesVisibility = 0
    DataRepresentation2.Scale = [1.0, 1.0, 1.0]
    DataRepresentation2.SelectionCellLabelJustification = 'Left'
    DataRepresentation2.DiffuseColor = [1.0, 1.0, 1.0]
    DataRepresentation2.Shade = 0
    DataRepresentation2.SelectionCellLabelOpacity = 1.0
    DataRepresentation2.CubeAxesInertia = 1
    DataRepresentation2.Source = "Arrow"
    DataRepresentation2.Source.Invert = 0
    DataRepresentation2.Masking = 0
    DataRepresentation2.Opacity = 1.0
    DataRepresentation2.LineWidth = 1.0
    DataRepresentation2.MeshVisibility = 0
    DataRepresentation2.Visibility = 1
    DataRepresentation2.SelectionCellLabelFontSize = 18
    DataRepresentation2.CubeAxesCornerOffset = 0.0
    DataRepresentation2.SelectionPointLabelJustification = 'Left'
    DataRepresentation2.SelectionPointLabelVisibility = 0
    DataRepresentation2.SelectOrientationVectors = ''
    DataRepresentation2.CubeAxesTickLocation = 'Inside'
    DataRepresentation2.BackfaceDiffuseColor = [1.0, 1.0, 1.0]
    DataRepresentation2.CubeAxesYAxisVisibility = 1
    DataRepresentation2.SelectionPointLabelFontFamily = 'Arial'
    DataRepresentation2.Source.ShaftResolution = 6
    DataRepresentation2.CubeAxesFlyMode = 'Closest Triad'
    DataRepresentation2.SelectScaleArray = ''
    DataRepresentation2.CubeAxesYTitle = 'Y-Axis'
    DataRepresentation2.ColorAttributeType = 'POINT_DATA'
    DataRepresentation2.SpecularPower = 100.0
    DataRepresentation2.Texture = []
    DataRepresentation2.SelectionCellLabelShadow = 0
    DataRepresentation2.AmbientColor = [1.0, 1.0, 1.0]
    DataRepresentation2.MapScalars = 1
    DataRepresentation2.PointSize = 2.0
    DataRepresentation2.Source.TipLength = 0.34999999999999998
    DataRepresentation2.SelectionCellLabelFormat = ''
    DataRepresentation2.Scaling = 0
    DataRepresentation2.StaticMode = 0
    DataRepresentation2.SelectionCellLabelColor = [0.0, 1.0, 0.0]
    DataRepresentation2.SliceMode = 'XY Plane'
    DataRepresentation2.Source.TipRadius = 0.10000000000000001
    DataRepresentation2.EdgeColor = [0.0, 0.0, 0.50000762951094835]
    DataRepresentation2.CubeAxesXAxisTickVisibility = 1
    DataRepresentation2.SelectionCellLabelVisibility = 0
    DataRepresentation2.NonlinearSubdivisionLevel = 1
    DataRepresentation2.CubeAxesXAxisRange = [0.0, 1.0]
    DataRepresentation2.Representation = 'Volume'
    DataRepresentation2.CubeAxesYAxisRange = [0.0, 1.0]
    DataRepresentation2.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
    DataRepresentation2.Orientation = [0.0, 0.0, 0.0]
    DataRepresentation2.CubeAxesEnableCustomAxisRange = 0
    DataRepresentation2.CubeAxesXTitle = 'X-Axis'
    DataRepresentation2.ScalarOpacityUnitDistance = 1.7320508075688776
    DataRepresentation2.BackfaceOpacity = 1.0
    DataRepresentation2.SelectionPointLabelFontSize = 18
    DataRepresentation2.SelectionCellFieldDataArrayName = 'vtkOriginalCellIds'
    DataRepresentation2.SelectionColor = [1.0, 0.0, 1.0]
    DataRepresentation2.Ambient = 0.0
    DataRepresentation2.VolumeRenderingMode = 'Smart'
    DataRepresentation2.CubeAxesXAxisMinorTickVisibility = 1
    DataRepresentation2.ScaleFactor = 12.700000000000001
    DataRepresentation2.BackfaceAmbientColor = [1.0, 1.0, 1.0]
    DataRepresentation2.Slice = 0
    DataRepresentation2.Source.ShaftRadius = 0.029999999999999999
    DataRepresentation2.ScalarOpacityFunction = a1_Result_PiecewiseFunction
    DataRepresentation2.SelectMaskArray = ''
    DataRepresentation2.SelectionLineWidth = 2.0
    DataRepresentation2.CubeAxesZAxisMinorTickVisibility = 1
    DataRepresentation2.CubeAxesXAxisVisibility = 1
    DataRepresentation2.Interpolation = 'Gouraud'
    DataRepresentation2.SelectionCellLabelFontFamily = 'Arial'
    DataRepresentation2.SelectionCellLabelItalic = 0
    DataRepresentation2.CubeAxesYAxisMinorTickVisibility = 1
    DataRepresentation2.InterpolationType = 'Linear'
    DataRepresentation2.CubeAxesZGridLines = 0
    DataRepresentation2.SelectionPointLabelFormat = ''
    DataRepresentation2.SelectionPointLabelOpacity = 1.0
    DataRepresentation2.Pickable = 1
    DataRepresentation2.CustomBoundsActive = [0, 0, 0]
    DataRepresentation2.SelectionRepresentation = 'Wireframe'
    DataRepresentation2.SelectionPointLabelBold = 0
    DataRepresentation2.ColorArrayName = 'Result'
    DataRepresentation2.SelectionPointLabelItalic = 0
    DataRepresentation2.AllowSpecularHighlightingWithScalarColoring = 0
    DataRepresentation2.SpecularColor = [1.0, 1.0, 1.0]
    DataRepresentation2.LookupTable = a1_Result_PVLookupTable
    DataRepresentation2.SelectionPointSize = 5.0
    DataRepresentation2.SelectionCellLabelBold = 0
    DataRepresentation2.Orient = 0
    

    for writer in cp_writers:
        if timestep % writer.cpFrequency == 0 or datadescription.GetForceOutput() == True:
            writer.FileName = writer.cpFileName.replace("%t", str(timestep))
            writer.UpdatePipeline()

    if False : # rescale data range
        import math
        for view in cp_views:
            if timestep % view.cpFrequency == 0 or datadescription.GetForceOutput() == True:
                reps = view.Representations
                for rep in reps:
                    if hasattr(rep, 'Visibility') and rep.Visibility == 1 and hasattr(rep, 'MapScalars') and rep.MapScalars != '':
                        input = rep.Input
                        input.UpdatePipeline() #make sure range is up-to-date
                        lut = rep.LookupTable
                        if lut == None:
                            continue
                        if rep.ColorAttributeType == 'POINT_DATA':
                            datainformation = input.GetPointDataInformation()
                        elif rep.ColorAttributeType == 'CELL_DATA':
                            datainformation = input.GetCellDataInformation()
                        else:
                            print 'something strange with color attribute type', rep.ColorAttributeType

                        if lut.VectorMode != 'Magnitude' or                            datainformation.GetArray(rep.ColorArrayName).GetNumberOfComponents() == 1:
                            datarange = datainformation.GetArray(rep.ColorArrayName).GetRange(lut.VectorComponent)
                        else:
                            datarange = [0,0]
                            for i in range(datainformation.GetArray(rep.ColorArrayName).GetNumberOfComponents()):
                                for j in range(2):
                                    datarange[j] += datainformation.GetArray(rep.ColorArrayName).GetRange(i)[j]*datainformation.GetArray(rep.ColorArrayName).GetRange(i)[j]
                            datarange[0] = math.sqrt(datarange[0])
                            datarange[1] = math.sqrt(datarange[1])

                        rgbpoints = lut.RGBPoints.GetData()
                        numpts = len(rgbpoints)/4
                        minvalue = min(datarange[0], rgbpoints[0])
                        maxvalue = max(datarange[1], rgbpoints[(numpts-1)*4])
                        if minvalue != rgbpoints[0] or maxvalue != rgbpoints[(numpts-1)*4]:
                            # rescale all of the points
                            oldrange = rgbpoints[(numpts-1)*4] - rgbpoints[0]
                            newrange = maxvalue - minvalue
                            newrgbpoints = list(rgbpoints)
                            for v in range(numpts):
                                newrgbpoints[v*4] = minvalue+(rgbpoints[v*4] - rgbpoints[0])*newrange/oldrange

                            lut.RGBPoints.SetData(newrgbpoints)

    for view in cp_views:
        if timestep % view.cpFrequency == 0 or datadescription.GetForceOutput() == True:
            fname = view.cpFileName
            fname = fname.replace("%t", str(timestep))
            if view.cpFitToScreen != 0:
                if view.IsA("vtkSMRenderViewProxy") == True:
                    view.ResetCamera()
                elif view.IsA("vtkSMContextViewProxy") == True:
                    view.ResetDisplay()
                else:
                    print ' do not know what to do with a ', view.GetClassName()

            WriteImage(fname, view, Magnification=view.cpMagnification)


    # explicitly delete the proxies -- we do it this way to avoid problems with prototypes
    tobedeleted = GetNextProxyToDelete()
    while tobedeleted != None:
        Delete(tobedeleted)
        tobedeleted = GetNextProxyToDelete()

def GetNextProxyToDelete():
    proxyiterator = servermanager.ProxyIterator()
    for proxy in proxyiterator:
        group = proxyiterator.GetGroup()
        if group.find("prototypes") != -1:
            continue
        if group != 'timekeeper' and group.find("pq_helper_proxies") == -1 :
            return proxy
    return None

def CreateProducer(datadescription, gridname):
    "Creates a producer proxy for the grid"
    if not datadescription.GetInputDescriptionByName(gridname):
        raise RuntimeError, "Simulation input name '%s' does not exist" % gridname
    grid = datadescription.GetInputDescriptionByName(gridname).GetGrid()
    producer = PVTrivialProducer()
    producer.GetClientSideObject().SetOutput(grid)
    if grid.IsA("vtkImageData") == True or grid.IsA("vtkStructuredGrid") == True or grid.IsA("vtkRectilinearGrid") == True:
        extent = datadescription.GetInputDescriptionByName(gridname).GetWholeExtent()
        producer.WholeExtent= [ extent[0], extent[1], extent[2], extent[3], extent[4], extent[5] ]

    producer.UpdatePipeline()
    return producer


def CreateWriter(proxy_ctor, filename, freq, cp_writers):
    writer = proxy_ctor()
    writer.FileName = filename
    writer.add_attribute("cpFrequency", freq)
    writer.add_attribute("cpFileName", filename)
    cp_writers.append(writer)
    return writer

def CreateView(proxy_ctor, filename, freq, fittoscreen, magnification, cp_views):
    view = proxy_ctor()
    view.add_attribute("cpFileName", filename)
    view.add_attribute("cpFrequency", freq)
    view.add_attribute("cpFileName", filename)
    view.add_attribute("cpFitToScreen", fittoscreen)
    view.add_attribute("cpMagnification", magnification)
    cp_views.append(view)
    return view

