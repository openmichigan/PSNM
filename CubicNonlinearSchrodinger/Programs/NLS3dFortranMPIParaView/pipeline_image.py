
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

    RenderView1 = CreateView( CreateRenderView, "image_%t.png", 1, 0, 1, cp_views )
    RenderView1.LightSpecularColor = [1.0, 1.0, 1.0]
    RenderView1.KeyLightAzimuth = 10.0
    RenderView1.UseTexturedBackground = 0
    RenderView1.UseLight = 1
    RenderView1.CameraPosition = [3.3432480226218373, 1.4248299445682202, 2.0022591914156083]
    RenderView1.FillLightKFRatio = 3.0
    RenderView1.Background2 = [0.0, 0.0, 0.16500000000000001]
    RenderView1.FillLightAzimuth = -10.0
    RenderView1.LODResolution = 50.0
    RenderView1.BackgroundTexture = []
    RenderView1.InteractionMode = '3D'
    RenderView1.StencilCapable = 1
    RenderView1.LightIntensity = 1.0
    RenderView1.CameraFocalPoint = [0.50000000000000089, 0.49999999999999911, 0.49999999999999983]
    RenderView1.ImageReductionFactor = 2
    RenderView1.CameraViewAngle = 30.0
    RenderView1.CameraParallelScale = 0.8660254037844386
    RenderView1.EyeAngle = 2.0
    RenderView1.HeadLightKHRatio = 3.0
    RenderView1.StereoRender = 0
    RenderView1.KeyLightIntensity = 0.75
    RenderView1.BackLightAzimuth = 110.0
    RenderView1.OrientationAxesInteractivity = 0
    RenderView1.UseInteractiveRenderingForSceenshots = 0
    RenderView1.UseOffscreenRendering = 0
    RenderView1.Background = [0.31999694819562063, 0.34000152590218968, 0.42999923704890519]
    RenderView1.UseOffscreenRenderingForScreenshots = 0
    RenderView1.NonInteractiveRenderDelay = 2
    RenderView1.CenterOfRotation = [0.5, 0.5, 0.5]
    RenderView1.CameraParallelProjection = 0
    RenderView1.CompressorConfig = 'vtkSquirtCompressor 0 3'
    RenderView1.HeadLightWarmth = 0.40000000000000002
    RenderView1.MaximumNumberOfPeels = 4
    RenderView1.LightDiffuseColor = [1.0, 1.0, 1.0]
    RenderView1.StereoType = 'Red-Blue'
    RenderView1.DepthPeeling = 1
    RenderView1.BackLightKBRatio = 3.5
    RenderView1.StereoCapableWindow = 1
    RenderView1.CameraViewUp = [-0.24173636424191886, 0.96102928228546003, -0.13411282113575904]
    RenderView1.LightType = 'HeadLight'
    RenderView1.LightAmbientColor = [1.0, 1.0, 1.0]
    RenderView1.RemoteRenderThreshold = 3.0
    RenderView1.CacheKey = 0.0
    RenderView1.UseCache = 0
    RenderView1.KeyLightElevation = 50.0
    RenderView1.CenterAxesVisibility = 0
    RenderView1.MaintainLuminance = 1
    RenderView1.StillRenderImageReductionFactor = 1
    RenderView1.BackLightWarmth = 0.5
    RenderView1.FillLightElevation = -75.0
    RenderView1.MultiSamples = 0
    RenderView1.FillLightWarmth = 0.40000000000000002
    RenderView1.AlphaBitPlanes = 1
    RenderView1.LightSwitch = 0
    RenderView1.OrientationAxesVisibility = 1
    RenderView1.CameraClippingRange = [1.7453949790250534, 5.3690501667232082]
    RenderView1.BackLightElevation = 0.0
    RenderView1.ViewTime = 0.0
    RenderView1.OrientationAxesOutlineColor = [1.0, 1.0, 1.0]
    RenderView1.LODThreshold = 18.199999999999999
    RenderView1.CollectGeometryThreshold = 100.0
    RenderView1.UseGradientBackground = 0
    RenderView1.KeyLightWarmth = 0.59999999999999998
    RenderView1.OrientationAxesLabelColor = [1.0, 1.0, 1.0]
    
    filename_ = CreateProducer( datadescription, "input" )
    
    CellDatatoPointData3 = CellDatatoPointData( guiName="CellDatatoPointData3", PieceInvariant=0, PassCellData=1 )
    
    Calculator3 = Calculator( guiName="Calculator3", Function='mag(u_complex_X*iHat+u_complex_Y*jHat)', ReplacementValue=0.0, ResultArrayName='Result', ReplaceInvalidResults=1, AttributeMode='point_data', CoordinateResults=0 )
    
    SetActiveSource(filename_)
    OutlineCorners1 = OutlineCorners( guiName="OutlineCorners1", CornerFactor=0.20000000000000001 )
    
    a2_u_PiecewiseFunction = CreatePiecewiseFunction( Points=[0.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0] )
    
    a1_Result_PiecewiseFunction = CreatePiecewiseFunction( Points=[1.2717775183445376e-20, 0.0, 0.5, 0.0, 0.80058018366957484, 1.0, 0.5, 0.0] )
    
    a2_u_PVLookupTable = GetLookupTableForArray( "u", 2, Discretize=1, RGBPoints=[6.1893147673285288e-20, 0.23000000000000001, 0.29899999999999999, 0.754, 0.99998800036605051, 0.70599999999999996, 0.016, 0.14999999999999999], UseLogScale=0, VectorComponent=0, NanColor=[0.25, 0.0, 0.0], NumberOfTableValues=256, ColorSpace='Diverging', VectorMode='Magnitude', HSVWrap=0, ScalarRangeInitialized=1.0, LockScalarRange=0 )
    
    a1_Result_PVLookupTable = GetLookupTableForArray( "Result", 1, Discretize=1, RGBPoints=[1.2717775183445376e-20, 0.23000000000000001, 0.29899999999999999, 0.754, 0.80058018366957484, 0.70599999999999996, 0.016, 0.14999999999999999], UseLogScale=0, VectorComponent=0, NanColor=[0.25, 0.0, 0.0], NumberOfTableValues=256, ColorSpace='Diverging', VectorMode='Magnitude', HSVWrap=0, ScalarRangeInitialized=1.0, LockScalarRange=0 )
    
    SetActiveSource(filename_)
    DataRepresentation11 = Show()
    DataRepresentation11.CubeAxesZAxisVisibility = 1
    DataRepresentation11.SelectionPointLabelColor = [0.5, 0.5, 0.5]
    DataRepresentation11.SelectionPointFieldDataArrayName = 'vtkOriginalPointIds'
    DataRepresentation11.SuppressLOD = 0
    DataRepresentation11.CubeAxesXGridLines = 0
    DataRepresentation11.CubeAxesYAxisTickVisibility = 1
    DataRepresentation11.CubeAxesColor = [1.0, 1.0, 1.0]
    DataRepresentation11.Position = [0.0, 0.0, 0.0]
    DataRepresentation11.BackfaceRepresentation = 'Follow Frontface'
    DataRepresentation11.SelectionOpacity = 1.0
    DataRepresentation11.SelectionPointLabelShadow = 0
    DataRepresentation11.CubeAxesYGridLines = 0
    DataRepresentation11.CubeAxesZAxisRange = [0.0, 1.0]
    DataRepresentation11.OrientationMode = 'Direction'
    DataRepresentation11.Source.TipResolution = 6
    DataRepresentation11.ScaleMode = 'No Data Scaling Off'
    DataRepresentation11.Diffuse = 1.0
    DataRepresentation11.SelectionUseOutline = 0
    DataRepresentation11.CubeAxesZTitle = 'Z-Axis'
    DataRepresentation11.Specular = 0.10000000000000001
    DataRepresentation11.SelectionVisibility = 1
    DataRepresentation11.InterpolateScalarsBeforeMapping = 1
    DataRepresentation11.CubeAxesZAxisTickVisibility = 1
    DataRepresentation11.Origin = [0.0, 0.0, 0.0]
    DataRepresentation11.CubeAxesVisibility = 0
    DataRepresentation11.Scale = [1.0, 1.0, 1.0]
    DataRepresentation11.SelectionCellLabelJustification = 'Left'
    DataRepresentation11.DiffuseColor = [1.0, 1.0, 1.0]
    DataRepresentation11.Shade = 0
    DataRepresentation11.SelectionCellLabelOpacity = 1.0
    DataRepresentation11.CubeAxesInertia = 1
    DataRepresentation11.Source = "Arrow"
    DataRepresentation11.Source.Invert = 0
    DataRepresentation11.Masking = 0
    DataRepresentation11.Opacity = 1.0
    DataRepresentation11.LineWidth = 1.0
    DataRepresentation11.MeshVisibility = 0
    DataRepresentation11.Visibility = 0
    DataRepresentation11.SelectionCellLabelFontSize = 18
    DataRepresentation11.CubeAxesCornerOffset = 0.0
    DataRepresentation11.SelectionPointLabelJustification = 'Left'
    DataRepresentation11.SelectionPointLabelVisibility = 0
    DataRepresentation11.SelectOrientationVectors = ''
    DataRepresentation11.CubeAxesTickLocation = 'Inside'
    DataRepresentation11.BackfaceDiffuseColor = [1.0, 1.0, 1.0]
    DataRepresentation11.CubeAxesYAxisVisibility = 1
    DataRepresentation11.SelectionPointLabelFontFamily = 'Arial'
    DataRepresentation11.Source.ShaftResolution = 6
    DataRepresentation11.CubeAxesFlyMode = 'Closest Triad'
    DataRepresentation11.SelectScaleArray = ''
    DataRepresentation11.CubeAxesYTitle = 'Y-Axis'
    DataRepresentation11.ColorAttributeType = 'POINT_DATA'
    DataRepresentation11.SpecularPower = 100.0
    DataRepresentation11.Texture = []
    DataRepresentation11.SelectionCellLabelShadow = 0
    DataRepresentation11.AmbientColor = [1.0, 1.0, 1.0]
    DataRepresentation11.MapScalars = 1
    DataRepresentation11.PointSize = 2.0
    DataRepresentation11.Source.TipLength = 0.34999999999999998
    DataRepresentation11.SelectionCellLabelFormat = ''
    DataRepresentation11.Scaling = 0
    DataRepresentation11.StaticMode = 0
    DataRepresentation11.SelectionCellLabelColor = [0.0, 1.0, 0.0]
    DataRepresentation11.SliceMode = 'XY Plane'
    DataRepresentation11.Source.TipRadius = 0.10000000000000001
    DataRepresentation11.EdgeColor = [0.0, 0.0, 0.50000762951094835]
    DataRepresentation11.CubeAxesXAxisTickVisibility = 1
    DataRepresentation11.SelectionCellLabelVisibility = 0
    DataRepresentation11.NonlinearSubdivisionLevel = 1
    DataRepresentation11.CubeAxesXAxisRange = [0.0, 1.0]
    DataRepresentation11.Representation = 'Outline'
    DataRepresentation11.CubeAxesYAxisRange = [0.0, 1.0]
    DataRepresentation11.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
    DataRepresentation11.Orientation = [0.0, 0.0, 0.0]
    DataRepresentation11.CubeAxesEnableCustomAxisRange = 0
    DataRepresentation11.CubeAxesXTitle = 'X-Axis'
    DataRepresentation11.ScalarOpacityUnitDistance = 0.054126587736527426
    DataRepresentation11.BackfaceOpacity = 1.0
    DataRepresentation11.SelectionPointLabelFontSize = 18
    DataRepresentation11.SelectionCellFieldDataArrayName = 'u_complex'
    DataRepresentation11.SelectionColor = [1.0, 0.0, 1.0]
    DataRepresentation11.Ambient = 0.0
    DataRepresentation11.VolumeRenderingMode = 'Smart'
    DataRepresentation11.CubeAxesXAxisMinorTickVisibility = 1
    DataRepresentation11.ScaleFactor = 0.10000000000000001
    DataRepresentation11.BackfaceAmbientColor = [1.0, 1.0, 1.0]
    DataRepresentation11.Slice = 0
    DataRepresentation11.Source.ShaftRadius = 0.029999999999999999
    DataRepresentation11.ScalarOpacityFunction = []
    DataRepresentation11.SelectMaskArray = ''
    DataRepresentation11.SelectionLineWidth = 2.0
    DataRepresentation11.CubeAxesZAxisMinorTickVisibility = 1
    DataRepresentation11.CubeAxesXAxisVisibility = 1
    DataRepresentation11.Interpolation = 'Gouraud'
    DataRepresentation11.SelectionCellLabelFontFamily = 'Arial'
    DataRepresentation11.SelectionCellLabelItalic = 0
    DataRepresentation11.CubeAxesYAxisMinorTickVisibility = 1
    DataRepresentation11.InterpolationType = 'Linear'
    DataRepresentation11.CubeAxesZGridLines = 0
    DataRepresentation11.SelectionPointLabelFormat = ''
    DataRepresentation11.SelectionPointLabelOpacity = 1.0
    DataRepresentation11.Pickable = 1
    DataRepresentation11.CustomBoundsActive = [0, 0, 0]
    DataRepresentation11.SelectionRepresentation = 'Wireframe'
    DataRepresentation11.SelectionPointLabelBold = 0
    DataRepresentation11.ColorArrayName = ''
    DataRepresentation11.SelectionPointLabelItalic = 0
    DataRepresentation11.AllowSpecularHighlightingWithScalarColoring = 0
    DataRepresentation11.SpecularColor = [1.0, 1.0, 1.0]
    DataRepresentation11.LookupTable = []
    DataRepresentation11.SelectionPointSize = 5.0
    DataRepresentation11.SelectionCellLabelBold = 0
    DataRepresentation11.Orient = 0
    
    SetActiveSource(CellDatatoPointData3)
    DataRepresentation12 = Show()
    DataRepresentation12.CubeAxesZAxisVisibility = 1
    DataRepresentation12.SelectionPointLabelColor = [0.5, 0.5, 0.5]
    DataRepresentation12.SelectionPointFieldDataArrayName = 'u_complex'
    DataRepresentation12.SuppressLOD = 0
    DataRepresentation12.CubeAxesXGridLines = 0
    DataRepresentation12.CubeAxesYAxisTickVisibility = 1
    DataRepresentation12.CubeAxesColor = [1.0, 1.0, 1.0]
    DataRepresentation12.Position = [0.0, 0.0, 0.0]
    DataRepresentation12.BackfaceRepresentation = 'Follow Frontface'
    DataRepresentation12.SelectionOpacity = 1.0
    DataRepresentation12.SelectionPointLabelShadow = 0
    DataRepresentation12.CubeAxesYGridLines = 0
    DataRepresentation12.CubeAxesZAxisRange = [0.0, 1.0]
    DataRepresentation12.OrientationMode = 'Direction'
    DataRepresentation12.Source.TipResolution = 6
    DataRepresentation12.ScaleMode = 'No Data Scaling Off'
    DataRepresentation12.Diffuse = 1.0
    DataRepresentation12.SelectionUseOutline = 0
    DataRepresentation12.CubeAxesZTitle = 'Z-Axis'
    DataRepresentation12.Specular = 0.10000000000000001
    DataRepresentation12.SelectionVisibility = 1
    DataRepresentation12.InterpolateScalarsBeforeMapping = 1
    DataRepresentation12.CubeAxesZAxisTickVisibility = 1
    DataRepresentation12.Origin = [0.0, 0.0, 0.0]
    DataRepresentation12.CubeAxesVisibility = 0
    DataRepresentation12.Scale = [1.0, 1.0, 1.0]
    DataRepresentation12.SelectionCellLabelJustification = 'Left'
    DataRepresentation12.DiffuseColor = [1.0, 1.0, 1.0]
    DataRepresentation12.Shade = 0
    DataRepresentation12.SelectionCellLabelOpacity = 1.0
    DataRepresentation12.CubeAxesInertia = 1
    DataRepresentation12.Source = "Arrow"
    DataRepresentation12.Source.Invert = 0
    DataRepresentation12.Masking = 0
    DataRepresentation12.Opacity = 1.0
    DataRepresentation12.LineWidth = 1.0
    DataRepresentation12.MeshVisibility = 0
    DataRepresentation12.Visibility = 0
    DataRepresentation12.SelectionCellLabelFontSize = 18
    DataRepresentation12.CubeAxesCornerOffset = 0.0
    DataRepresentation12.SelectionPointLabelJustification = 'Left'
    DataRepresentation12.SelectionPointLabelVisibility = 0
    DataRepresentation12.SelectOrientationVectors = ''
    DataRepresentation12.CubeAxesTickLocation = 'Inside'
    DataRepresentation12.BackfaceDiffuseColor = [1.0, 1.0, 1.0]
    DataRepresentation12.CubeAxesYAxisVisibility = 1
    DataRepresentation12.SelectionPointLabelFontFamily = 'Arial'
    DataRepresentation12.Source.ShaftResolution = 6
    DataRepresentation12.CubeAxesFlyMode = 'Closest Triad'
    DataRepresentation12.SelectScaleArray = ''
    DataRepresentation12.CubeAxesYTitle = 'Y-Axis'
    DataRepresentation12.ColorAttributeType = 'POINT_DATA'
    DataRepresentation12.SpecularPower = 100.0
    DataRepresentation12.Texture = []
    DataRepresentation12.SelectionCellLabelShadow = 0
    DataRepresentation12.AmbientColor = [1.0, 1.0, 1.0]
    DataRepresentation12.MapScalars = 1
    DataRepresentation12.PointSize = 2.0
    DataRepresentation12.Source.TipLength = 0.34999999999999998
    DataRepresentation12.SelectionCellLabelFormat = ''
    DataRepresentation12.Scaling = 0
    DataRepresentation12.StaticMode = 0
    DataRepresentation12.SelectionCellLabelColor = [0.0, 1.0, 0.0]
    DataRepresentation12.SliceMode = 'XY Plane'
    DataRepresentation12.Source.TipRadius = 0.10000000000000001
    DataRepresentation12.EdgeColor = [0.0, 0.0, 0.50000762951094835]
    DataRepresentation12.CubeAxesXAxisTickVisibility = 1
    DataRepresentation12.SelectionCellLabelVisibility = 0
    DataRepresentation12.NonlinearSubdivisionLevel = 1
    DataRepresentation12.CubeAxesXAxisRange = [0.0, 1.0]
    DataRepresentation12.Representation = 'Outline'
    DataRepresentation12.CubeAxesYAxisRange = [0.0, 1.0]
    DataRepresentation12.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
    DataRepresentation12.Orientation = [0.0, 0.0, 0.0]
    DataRepresentation12.CubeAxesEnableCustomAxisRange = 0
    DataRepresentation12.CubeAxesXTitle = 'X-Axis'
    DataRepresentation12.ScalarOpacityUnitDistance = 0.054126587736527426
    DataRepresentation12.BackfaceOpacity = 1.0
    DataRepresentation12.SelectionPointLabelFontSize = 18
    DataRepresentation12.SelectionCellFieldDataArrayName = 'u_complex'
    DataRepresentation12.SelectionColor = [1.0, 0.0, 1.0]
    DataRepresentation12.Ambient = 0.0
    DataRepresentation12.VolumeRenderingMode = 'Smart'
    DataRepresentation12.CubeAxesXAxisMinorTickVisibility = 1
    DataRepresentation12.ScaleFactor = 0.10000000000000001
    DataRepresentation12.BackfaceAmbientColor = [1.0, 1.0, 1.0]
    DataRepresentation12.Slice = 0
    DataRepresentation12.Source.ShaftRadius = 0.029999999999999999
    DataRepresentation12.ScalarOpacityFunction = []
    DataRepresentation12.SelectMaskArray = ''
    DataRepresentation12.SelectionLineWidth = 2.0
    DataRepresentation12.CubeAxesZAxisMinorTickVisibility = 1
    DataRepresentation12.CubeAxesXAxisVisibility = 1
    DataRepresentation12.Interpolation = 'Gouraud'
    DataRepresentation12.SelectionCellLabelFontFamily = 'Arial'
    DataRepresentation12.SelectionCellLabelItalic = 0
    DataRepresentation12.CubeAxesYAxisMinorTickVisibility = 1
    DataRepresentation12.InterpolationType = 'Linear'
    DataRepresentation12.CubeAxesZGridLines = 0
    DataRepresentation12.SelectionPointLabelFormat = ''
    DataRepresentation12.SelectionPointLabelOpacity = 1.0
    DataRepresentation12.Pickable = 1
    DataRepresentation12.CustomBoundsActive = [0, 0, 0]
    DataRepresentation12.SelectionRepresentation = 'Wireframe'
    DataRepresentation12.SelectionPointLabelBold = 0
    DataRepresentation12.ColorArrayName = ''
    DataRepresentation12.SelectionPointLabelItalic = 0
    DataRepresentation12.AllowSpecularHighlightingWithScalarColoring = 0
    DataRepresentation12.SpecularColor = [1.0, 1.0, 1.0]
    DataRepresentation12.LookupTable = []
    DataRepresentation12.SelectionPointSize = 5.0
    DataRepresentation12.SelectionCellLabelBold = 0
    DataRepresentation12.Orient = 0
    
    SetActiveSource(Calculator3)
    DataRepresentation13 = Show()
    DataRepresentation13.CubeAxesZAxisVisibility = 1
    DataRepresentation13.SelectionPointLabelColor = [0.5, 0.5, 0.5]
    DataRepresentation13.SelectionPointFieldDataArrayName = 'Result'
    DataRepresentation13.SuppressLOD = 0
    DataRepresentation13.CubeAxesXGridLines = 0
    DataRepresentation13.CubeAxesYAxisTickVisibility = 1
    DataRepresentation13.CubeAxesColor = [1.0, 1.0, 1.0]
    DataRepresentation13.Position = [0.0, 0.0, 0.0]
    DataRepresentation13.BackfaceRepresentation = 'Follow Frontface'
    DataRepresentation13.SelectionOpacity = 1.0
    DataRepresentation13.SelectionPointLabelShadow = 0
    DataRepresentation13.CubeAxesYGridLines = 0
    DataRepresentation13.CubeAxesZAxisRange = [0.0, 1.0]
    DataRepresentation13.OrientationMode = 'Direction'
    DataRepresentation13.Source.TipResolution = 6
    DataRepresentation13.ScaleMode = 'No Data Scaling Off'
    DataRepresentation13.Diffuse = 1.0
    DataRepresentation13.SelectionUseOutline = 0
    DataRepresentation13.CubeAxesZTitle = 'Z-Axis'
    DataRepresentation13.Specular = 0.10000000000000001
    DataRepresentation13.SelectionVisibility = 1
    DataRepresentation13.InterpolateScalarsBeforeMapping = 1
    DataRepresentation13.CubeAxesZAxisTickVisibility = 1
    DataRepresentation13.Origin = [0.0, 0.0, 0.0]
    DataRepresentation13.CubeAxesVisibility = 0
    DataRepresentation13.Scale = [1.0, 1.0, 1.0]
    DataRepresentation13.SelectionCellLabelJustification = 'Left'
    DataRepresentation13.DiffuseColor = [1.0, 1.0, 1.0]
    DataRepresentation13.Shade = 0
    DataRepresentation13.SelectionCellLabelOpacity = 1.0
    DataRepresentation13.CubeAxesInertia = 1
    DataRepresentation13.Source = "Arrow"
    DataRepresentation13.Source.Invert = 0
    DataRepresentation13.Masking = 0
    DataRepresentation13.Opacity = 1.0
    DataRepresentation13.LineWidth = 1.0
    DataRepresentation13.MeshVisibility = 0
    DataRepresentation13.Visibility = 1
    DataRepresentation13.SelectionCellLabelFontSize = 18
    DataRepresentation13.CubeAxesCornerOffset = 0.0
    DataRepresentation13.SelectionPointLabelJustification = 'Left'
    DataRepresentation13.SelectionPointLabelVisibility = 0
    DataRepresentation13.SelectOrientationVectors = ''
    DataRepresentation13.CubeAxesTickLocation = 'Inside'
    DataRepresentation13.BackfaceDiffuseColor = [1.0, 1.0, 1.0]
    DataRepresentation13.CubeAxesYAxisVisibility = 1
    DataRepresentation13.SelectionPointLabelFontFamily = 'Arial'
    DataRepresentation13.Source.ShaftResolution = 6
    DataRepresentation13.CubeAxesFlyMode = 'Closest Triad'
    DataRepresentation13.SelectScaleArray = ''
    DataRepresentation13.CubeAxesYTitle = 'Y-Axis'
    DataRepresentation13.ColorAttributeType = 'POINT_DATA'
    DataRepresentation13.SpecularPower = 100.0
    DataRepresentation13.Texture = []
    DataRepresentation13.SelectionCellLabelShadow = 0
    DataRepresentation13.AmbientColor = [1.0, 1.0, 1.0]
    DataRepresentation13.MapScalars = 1
    DataRepresentation13.PointSize = 2.0
    DataRepresentation13.Source.TipLength = 0.34999999999999998
    DataRepresentation13.SelectionCellLabelFormat = ''
    DataRepresentation13.Scaling = 0
    DataRepresentation13.StaticMode = 0
    DataRepresentation13.SelectionCellLabelColor = [0.0, 1.0, 0.0]
    DataRepresentation13.SliceMode = 'XY Plane'
    DataRepresentation13.Source.TipRadius = 0.10000000000000001
    DataRepresentation13.EdgeColor = [0.0, 0.0, 0.50000762951094835]
    DataRepresentation13.CubeAxesXAxisTickVisibility = 1
    DataRepresentation13.SelectionCellLabelVisibility = 0
    DataRepresentation13.NonlinearSubdivisionLevel = 1
    DataRepresentation13.CubeAxesXAxisRange = [0.0, 1.0]
    DataRepresentation13.Representation = 'Volume'
    DataRepresentation13.CubeAxesYAxisRange = [0.0, 1.0]
    DataRepresentation13.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
    DataRepresentation13.Orientation = [0.0, 0.0, 0.0]
    DataRepresentation13.CubeAxesEnableCustomAxisRange = 0
    DataRepresentation13.CubeAxesXTitle = 'X-Axis'
    DataRepresentation13.ScalarOpacityUnitDistance = 0.054126587736527426
    DataRepresentation13.BackfaceOpacity = 1.0
    DataRepresentation13.SelectionPointLabelFontSize = 18
    DataRepresentation13.SelectionCellFieldDataArrayName = 'u_complex'
    DataRepresentation13.SelectionColor = [1.0, 0.0, 1.0]
    DataRepresentation13.Ambient = 0.0
    DataRepresentation13.VolumeRenderingMode = 'Smart'
    DataRepresentation13.CubeAxesXAxisMinorTickVisibility = 1
    DataRepresentation13.ScaleFactor = 0.10000000000000001
    DataRepresentation13.BackfaceAmbientColor = [1.0, 1.0, 1.0]
    DataRepresentation13.Slice = 0
    DataRepresentation13.Source.ShaftRadius = 0.029999999999999999
    DataRepresentation13.ScalarOpacityFunction = a1_Result_PiecewiseFunction
    DataRepresentation13.SelectMaskArray = ''
    DataRepresentation13.SelectionLineWidth = 2.0
    DataRepresentation13.CubeAxesZAxisMinorTickVisibility = 1
    DataRepresentation13.CubeAxesXAxisVisibility = 1
    DataRepresentation13.Interpolation = 'Gouraud'
    DataRepresentation13.SelectionCellLabelFontFamily = 'Arial'
    DataRepresentation13.SelectionCellLabelItalic = 0
    DataRepresentation13.CubeAxesYAxisMinorTickVisibility = 1
    DataRepresentation13.InterpolationType = 'Linear'
    DataRepresentation13.CubeAxesZGridLines = 0
    DataRepresentation13.SelectionPointLabelFormat = ''
    DataRepresentation13.SelectionPointLabelOpacity = 1.0
    DataRepresentation13.Pickable = 1
    DataRepresentation13.CustomBoundsActive = [0, 0, 0]
    DataRepresentation13.SelectionRepresentation = 'Wireframe'
    DataRepresentation13.SelectionPointLabelBold = 0
    DataRepresentation13.ColorArrayName = 'Result'
    DataRepresentation13.SelectionPointLabelItalic = 0
    DataRepresentation13.AllowSpecularHighlightingWithScalarColoring = 0
    DataRepresentation13.SpecularColor = [1.0, 1.0, 1.0]
    DataRepresentation13.LookupTable = a1_Result_PVLookupTable
    DataRepresentation13.SelectionPointSize = 5.0
    DataRepresentation13.SelectionCellLabelBold = 0
    DataRepresentation13.Orient = 0
    
    SetActiveSource(OutlineCorners1)
    DataRepresentation14 = Show()
    DataRepresentation14.CubeAxesZAxisVisibility = 1
    DataRepresentation14.SelectionPointLabelColor = [0.5, 0.5, 0.5]
    DataRepresentation14.SelectionPointFieldDataArrayName = 'vtkOriginalPointIds'
    DataRepresentation14.SuppressLOD = 0
    DataRepresentation14.CubeAxesXGridLines = 0
    DataRepresentation14.CubeAxesYAxisTickVisibility = 1
    DataRepresentation14.CubeAxesColor = [1.0, 1.0, 1.0]
    DataRepresentation14.Position = [0.0, 0.0, 0.0]
    DataRepresentation14.BackfaceRepresentation = 'Follow Frontface'
    DataRepresentation14.SelectionOpacity = 1.0
    DataRepresentation14.SelectionPointLabelShadow = 0
    DataRepresentation14.CubeAxesYGridLines = 0
    DataRepresentation14.CubeAxesZAxisRange = [0.0, 1.0]
    DataRepresentation14.OrientationMode = 'Direction'
    DataRepresentation14.Source.TipResolution = 6
    DataRepresentation14.ScaleMode = 'No Data Scaling Off'
    DataRepresentation14.Diffuse = 1.0
    DataRepresentation14.SelectionUseOutline = 0
    DataRepresentation14.CubeAxesZTitle = 'Z-Axis'
    DataRepresentation14.Specular = 0.10000000000000001
    DataRepresentation14.SelectionVisibility = 1
    DataRepresentation14.InterpolateScalarsBeforeMapping = 1
    DataRepresentation14.CubeAxesZAxisTickVisibility = 1
    DataRepresentation14.Origin = [0.0, 0.0, 0.0]
    DataRepresentation14.CubeAxesVisibility = 0
    DataRepresentation14.Scale = [1.0, 1.0, 1.0]
    DataRepresentation14.SelectionCellLabelJustification = 'Left'
    DataRepresentation14.DiffuseColor = [1.0, 1.0, 1.0]
    DataRepresentation14.SelectionCellLabelOpacity = 1.0
    DataRepresentation14.Source = "Arrow"
    DataRepresentation14.Source.Invert = 0
    DataRepresentation14.Masking = 0
    DataRepresentation14.Opacity = 1.0
    DataRepresentation14.LineWidth = 1.0
    DataRepresentation14.MeshVisibility = 0
    DataRepresentation14.Visibility = 1
    DataRepresentation14.SelectionCellLabelFontSize = 18
    DataRepresentation14.CubeAxesCornerOffset = 0.0
    DataRepresentation14.SelectionPointLabelJustification = 'Left'
    DataRepresentation14.SelectionPointLabelVisibility = 0
    DataRepresentation14.SelectOrientationVectors = ''
    DataRepresentation14.CubeAxesTickLocation = 'Inside'
    DataRepresentation14.CubeAxesXAxisMinorTickVisibility = 1
    DataRepresentation14.CubeAxesYAxisVisibility = 1
    DataRepresentation14.SelectionPointLabelFontFamily = 'Arial'
    DataRepresentation14.Source.ShaftResolution = 6
    DataRepresentation14.CubeAxesFlyMode = 'Closest Triad'
    DataRepresentation14.SelectScaleArray = ''
    DataRepresentation14.CubeAxesYTitle = 'Y-Axis'
    DataRepresentation14.ColorAttributeType = 'POINT_DATA'
    DataRepresentation14.SpecularPower = 100.0
    DataRepresentation14.Texture = []
    DataRepresentation14.SelectionCellLabelShadow = 0
    DataRepresentation14.AmbientColor = [1.0, 1.0, 1.0]
    DataRepresentation14.MapScalars = 1
    DataRepresentation14.PointSize = 2.0
    DataRepresentation14.Source.TipLength = 0.34999999999999998
    DataRepresentation14.SelectionCellLabelFormat = ''
    DataRepresentation14.Scaling = 0
    DataRepresentation14.StaticMode = 0
    DataRepresentation14.SelectionCellLabelColor = [0.0, 1.0, 0.0]
    DataRepresentation14.Source.TipRadius = 0.10000000000000001
    DataRepresentation14.EdgeColor = [0.0, 0.0, 0.50000762951094835]
    DataRepresentation14.CubeAxesXAxisTickVisibility = 1
    DataRepresentation14.SelectionCellLabelVisibility = 0
    DataRepresentation14.NonlinearSubdivisionLevel = 1
    DataRepresentation14.CubeAxesXAxisRange = [0.0, 1.0]
    DataRepresentation14.Representation = 'Surface'
    DataRepresentation14.CubeAxesYAxisRange = [0.0, 1.0]
    DataRepresentation14.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
    DataRepresentation14.Orientation = [0.0, 0.0, 0.0]
    DataRepresentation14.CubeAxesEnableCustomAxisRange = 0
    DataRepresentation14.CubeAxesXTitle = 'X-Axis'
    DataRepresentation14.CubeAxesInertia = 1
    DataRepresentation14.BackfaceOpacity = 1.0
    DataRepresentation14.SelectionCellFieldDataArrayName = 'u_complex'
    DataRepresentation14.SelectionColor = [1.0, 0.0, 1.0]
    DataRepresentation14.Ambient = 0.0
    DataRepresentation14.SelectionPointLabelFontSize = 18
    DataRepresentation14.ScaleFactor = 0.10000000000000001
    DataRepresentation14.BackfaceAmbientColor = [1.0, 1.0, 1.0]
    DataRepresentation14.Source.ShaftRadius = 0.029999999999999999
    DataRepresentation14.SelectMaskArray = ''
    DataRepresentation14.SelectionLineWidth = 2.0
    DataRepresentation14.CubeAxesZAxisMinorTickVisibility = 1
    DataRepresentation14.CubeAxesXAxisVisibility = 1
    DataRepresentation14.Interpolation = 'Gouraud'
    DataRepresentation14.SelectionCellLabelFontFamily = 'Arial'
    DataRepresentation14.SelectionCellLabelItalic = 0
    DataRepresentation14.CubeAxesYAxisMinorTickVisibility = 1
    DataRepresentation14.CubeAxesZGridLines = 0
    DataRepresentation14.SelectionPointLabelFormat = ''
    DataRepresentation14.SelectionPointLabelOpacity = 1.0
    DataRepresentation14.BackfaceDiffuseColor = [1.0, 1.0, 1.0]
    DataRepresentation14.Pickable = 1
    DataRepresentation14.CustomBoundsActive = [0, 0, 0]
    DataRepresentation14.SelectionRepresentation = 'Wireframe'
    DataRepresentation14.SelectionPointLabelBold = 0
    DataRepresentation14.ColorArrayName = ''
    DataRepresentation14.SelectionPointLabelItalic = 0
    DataRepresentation14.AllowSpecularHighlightingWithScalarColoring = 0
    DataRepresentation14.SpecularColor = [1.0, 1.0, 1.0]
    DataRepresentation14.LookupTable = []
    DataRepresentation14.SelectionPointSize = 5.0
    DataRepresentation14.SelectionCellLabelBold = 0
    DataRepresentation14.Orient = 0
    

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

