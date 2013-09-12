
try: paraview.simple
except: from paraview.simple import *

from paraview import coprocessing


#--------------------------------------------------------------
# Code generated from cpstate.py to create the CoProcessor.


# ----------------------- CoProcessor definition -----------------------

def CreateCoProcessor():
  def _CreatePipeline(coprocessor, datadescription):
    class Pipeline:
      RenderView1 = coprocessor.CreateView( CreateRenderView, "image_%t.png", 1, 0, 1, 600, 600 )
      RenderView1.LightSpecularColor = [1.0, 1.0, 1.0]
      RenderView1.UseOutlineForLODRendering = 0
      RenderView1.KeyLightAzimuth = 10.0
      RenderView1.UseTexturedBackground = 0
      RenderView1.UseLight = 1
      RenderView1.CameraPosition = [0.5162288830998716, 0.5298241735361152, 4.548596535110469]
      RenderView1.FillLightKFRatio = 3.0
      RenderView1.Background2 = [0.0, 0.0, 0.16470588235294117]
      RenderView1.FillLightAzimuth = -10.0
      RenderView1.LODResolution = 0.5
      RenderView1.BackgroundTexture = []
      RenderView1.InteractionMode = '3D'
      RenderView1.StencilCapable = 1
      RenderView1.LightIntensity = 0.11
      RenderView1.CameraFocalPoint = [0.49999999999999956, 0.5000000000000012, 0.4999999999999978]
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
      RenderView1.Background = [0.0, 0.0, 0.0]
      RenderView1.UseOffscreenRenderingForScreenshots = 0
      RenderView1.NonInteractiveRenderDelay = 0.0
      RenderView1.CenterOfRotation = [0.5, 0.5, 0.5]
      RenderView1.CameraParallelProjection = 0
      RenderView1.CompressorConfig = 'vtkSquirtCompressor 0 3'
      RenderView1.HeadLightWarmth = 0.5
      RenderView1.MaximumNumberOfPeels = 4
      RenderView1.LightDiffuseColor = [1.0, 1.0, 1.0]
      RenderView1.StereoType = 'Red-Blue'
      RenderView1.DepthPeeling = 1
      RenderView1.BackLightKBRatio = 3.5
      RenderView1.StereoCapableWindow = 1
      RenderView1.CameraViewUp = [0.0032624074714670973, 0.9999674501658008, -0.007379383876510121]
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
      RenderView1.FillLightWarmth = 0.4
      RenderView1.AlphaBitPlanes = 1
      RenderView1.LightSwitch = 1
      RenderView1.OrientationAxesVisibility = 0
      RenderView1.CameraClippingRange = [3.001968716822371, 5.376172719590956]
      RenderView1.BackLightElevation = 0.0
      RenderView1.ViewTime = 0.0
      RenderView1.OrientationAxesOutlineColor = [1.0, 1.0, 1.0]
      RenderView1.LODThreshold = 5.0
      RenderView1.CollectGeometryThreshold = 100.0
      RenderView1.UseGradientBackground = 0
      RenderView1.KeyLightWarmth = 0.6
      RenderView1.OrientationAxesLabelColor = [1.0, 1.0, 1.0]
      
      filename_1_0_vti = coprocessor.CreateProducer( datadescription, "input" )
      
      CellDatatoPointData1 = CellDatatoPointData( guiName="CellDatatoPointData1", PieceInvariant=0, PassCellData=0 )
      
      a1_u_PiecewiseFunction = CreatePiecewiseFunction( Points=[9.299771875392746e-117, 0.0, 0.5, 0.0, 0.5, 1.0, 0.5, 0.0] )
      
      a1_u_PVLookupTable = GetLookupTableForArray( "u", 1, Discretize=1, RGBPoints=[9.299771875392746e-117, 0.0, 1.0, 1.0, 0.225, 0.0, 0.0, 1.0, 0.25, 0.0, 0.0, 0.5019607843137255, 0.275, 1.0, 0.0, 0.0, 0.5, 1.0, 1.0, 0.0], UseLogScale=0, VectorComponent=0, NanColor=[1.0, 1.0, 0.0], NumberOfTableValues=256, EnableOpacityMapping=0, ColorSpace='RGB', IndexedLookup=0, VectorMode='Magnitude', ScalarOpacityFunction=a1_u_PiecewiseFunction, HSVWrap=0, ScalarRangeInitialized=1.0, AllowDuplicateScalars=1, Annotations=['1', '1', '0', '0'], LockScalarRange=0 )
      
      ScalarBarWidgetRepresentation1 = CreateScalarBar( TextPosition=1, Title='u', Position2=[0.13, 0.5], TitleOpacity=1.0, TitleFontSize=12, NanAnnotation='NaN', TitleShadow=0, AutomaticLabelFormat=1, DrawAnnotations=1, TitleColor=[1.0, 1.0, 1.0], AspectRatio=20.0, NumberOfLabels=5, ComponentTitle='', Resizable=1, DrawNanAnnotation=0, TitleFontFamily='Arial', Visibility=0, LabelFontSize=12, LabelFontFamily='Arial', TitleItalic=0, LabelBold=0, LabelItalic=0, Enabled=0, LabelColor=[1.0, 1.0, 1.0], Position=[0.87, 0.25], Selectable=0, UseNonCompositedRenderer=1, LabelOpacity=1.0, TitleBold=0, LabelFormat='%-#6.3g', Orientation='Vertical', LabelShadow=0, LookupTable=a1_u_PVLookupTable, Repositionable=1 )
      GetRenderView().Representations.append(ScalarBarWidgetRepresentation1)
      
      SetActiveSource(filename_1_0_vti)
      DataRepresentation1 = Show()
      DataRepresentation1.CubeAxesZAxisVisibility = 1
      DataRepresentation1.SelectionPointLabelColor = [0.5, 0.5, 0.5]
      DataRepresentation1.SelectionPointFieldDataArrayName = 'vtkOriginalPointIds'
      DataRepresentation1.SuppressLOD = 0
      DataRepresentation1.CubeAxesXGridLines = 0
      DataRepresentation1.BlockVisibility = []
      DataRepresentation1.CubeAxesYAxisTickVisibility = 1
      DataRepresentation1.Position = [0.0, 0.0, 0.0]
      DataRepresentation1.BackfaceRepresentation = 'Follow Frontface'
      DataRepresentation1.SelectionOpacity = 1.0
      DataRepresentation1.SelectionPointLabelShadow = 0
      DataRepresentation1.CubeAxesYGridLines = 0
      DataRepresentation1.CubeAxesZAxisTickVisibility = 1
      DataRepresentation1.OrientationMode = 'Direction'
      DataRepresentation1.ScaleMode = 'No Data Scaling Off'
      DataRepresentation1.Diffuse = 1.0
      DataRepresentation1.SelectionUseOutline = 0
      DataRepresentation1.CubeAxesZTitle = 'Z-Axis'
      DataRepresentation1.Specular = 0.1
      DataRepresentation1.SelectionVisibility = 1
      DataRepresentation1.InterpolateScalarsBeforeMapping = 1
      DataRepresentation1.CustomRangeActive = [0, 0, 0]
      DataRepresentation1.Origin = [0.0, 0.0, 0.0]
      DataRepresentation1.CubeAxesVisibility = 0
      DataRepresentation1.Scale = [1.0, 1.0, 1.0]
      DataRepresentation1.SelectionCellLabelJustification = 'Left'
      DataRepresentation1.DiffuseColor = [1.0, 1.0, 1.0]
      DataRepresentation1.Shade = 0
      DataRepresentation1.SelectionCellLabelOpacity = 1.0
      DataRepresentation1.CubeAxesInertia = 1
      DataRepresentation1.Source = []
      DataRepresentation1.Masking = 0
      DataRepresentation1.Opacity = 1.0
      DataRepresentation1.LineWidth = 1.0
      DataRepresentation1.MeshVisibility = 0
      DataRepresentation1.Visibility = 0
      DataRepresentation1.SelectionCellLabelFontSize = 18
      DataRepresentation1.CubeAxesCornerOffset = 0.0
      DataRepresentation1.SelectionPointLabelJustification = 'Left'
      DataRepresentation1.OriginalBoundsRangeActive = [0, 0, 0]
      DataRepresentation1.SelectionPointLabelVisibility = 0
      DataRepresentation1.SelectOrientationVectors = ''
      DataRepresentation1.CubeAxesTickLocation = 'Inside'
      DataRepresentation1.BackfaceDiffuseColor = [1.0, 1.0, 1.0]
      DataRepresentation1.CubeAxesYLabelFormat = '%-#6.3g'
      DataRepresentation1.CubeAxesYAxisVisibility = 1
      DataRepresentation1.SelectionPointLabelFontFamily = 'Arial'
      DataRepresentation1.CubeAxesUseDefaultYTitle = 1
      DataRepresentation1.SelectScaleArray = ''
      DataRepresentation1.CubeAxesYTitle = 'Y-Axis'
      DataRepresentation1.ColorAttributeType = 'CELL_DATA'
      DataRepresentation1.AxesOrigin = [0.0, 0.0, 0.0]
      DataRepresentation1.UserTransform = [1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0]
      DataRepresentation1.SpecularPower = 100.0
      DataRepresentation1.Texture = []
      DataRepresentation1.SelectionCellLabelShadow = 0
      DataRepresentation1.AmbientColor = [1.0, 1.0, 1.0]
      DataRepresentation1.BlockOpacity = {}
      DataRepresentation1.MapScalars = 1
      DataRepresentation1.PointSize = 2.0
      DataRepresentation1.CubeAxesUseDefaultXTitle = 1
      DataRepresentation1.SelectionCellLabelFormat = ''
      DataRepresentation1.Scaling = 0
      DataRepresentation1.StaticMode = 0
      DataRepresentation1.SelectionCellLabelColor = [0.0, 1.0, 0.0]
      DataRepresentation1.SliceMode = 'XY Plane'
      DataRepresentation1.EdgeColor = [0.0, 0.0, 0.5000076295109483]
      DataRepresentation1.CubeAxesXAxisTickVisibility = 1
      DataRepresentation1.SelectionCellLabelVisibility = 0
      DataRepresentation1.NonlinearSubdivisionLevel = 1
      DataRepresentation1.CubeAxesColor = [1.0, 1.0, 1.0]
      DataRepresentation1.Representation = 'Outline'
      DataRepresentation1.CustomRange = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
      DataRepresentation1.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
      DataRepresentation1.Orientation = [0.0, 0.0, 0.0]
      DataRepresentation1.CubeAxesXTitle = 'X-Axis'
      DataRepresentation1.ScalarOpacityUnitDistance = 0.027063293868263713
      DataRepresentation1.BackfaceOpacity = 1.0
      DataRepresentation1.SelectionPointLabelFontSize = 18
      DataRepresentation1.SelectionCellFieldDataArrayName = 'u'
      DataRepresentation1.SelectionColor = [1.0, 0.0, 1.0]
      DataRepresentation1.BlockColor = {}
      DataRepresentation1.Ambient = 0.0
      DataRepresentation1.VolumeRenderingMode = 'Smart'
      DataRepresentation1.CubeAxesXAxisMinorTickVisibility = 1
      DataRepresentation1.ScaleFactor = 0.1
      DataRepresentation1.BackfaceAmbientColor = [1.0, 1.0, 1.0]
      DataRepresentation1.Slice = 0
      DataRepresentation1.ScalarOpacityFunction = a1_u_PiecewiseFunction
      DataRepresentation1.SelectMaskArray = ''
      DataRepresentation1.SelectionLineWidth = 2.0
      DataRepresentation1.CubeAxesZAxisMinorTickVisibility = 1
      DataRepresentation1.CubeAxesXAxisVisibility = 1
      DataRepresentation1.CubeAxesXLabelFormat = '%-#6.3g'
      DataRepresentation1.Interpolation = 'Gouraud'
      DataRepresentation1.CubeAxesZLabelFormat = '%-#6.3g'
      DataRepresentation1.SelectionCellLabelFontFamily = 'Arial'
      DataRepresentation1.SelectionCellLabelItalic = 0
      DataRepresentation1.CubeAxesYAxisMinorTickVisibility = 1
      DataRepresentation1.InterpolationType = 'Cubic'
      DataRepresentation1.CubeAxesZGridLines = 0
      DataRepresentation1.SelectionPointLabelFormat = ''
      DataRepresentation1.SelectionPointLabelOpacity = 1.0
      DataRepresentation1.UseAxesOrigin = 0
      DataRepresentation1.CubeAxesFlyMode = 'Closest Triad'
      DataRepresentation1.Pickable = 1
      DataRepresentation1.CustomBoundsActive = [0, 0, 0]
      DataRepresentation1.CubeAxesGridLineLocation = 'All Faces'
      DataRepresentation1.SelectionRepresentation = 'Wireframe'
      DataRepresentation1.SelectionPointLabelBold = 0
      DataRepresentation1.ColorArrayName = ('CELL_DATA', 'u')
      DataRepresentation1.SelectionPointLabelItalic = 0
      DataRepresentation1.AllowSpecularHighlightingWithScalarColoring = 0
      DataRepresentation1.SpecularColor = [1.0, 1.0, 1.0]
      DataRepresentation1.CubeAxesUseDefaultZTitle = 1
      DataRepresentation1.LookupTable = a1_u_PVLookupTable
      DataRepresentation1.SelectionPointSize = 5.0
      DataRepresentation1.SelectionCellLabelBold = 0
      DataRepresentation1.Orient = 0
      
      SetActiveSource(CellDatatoPointData1)
      DataRepresentation2 = Show()
      DataRepresentation2.CubeAxesZAxisVisibility = 1
      DataRepresentation2.SelectionPointLabelColor = [0.5, 0.5, 0.5]
      DataRepresentation2.SelectionPointFieldDataArrayName = 'u'
      DataRepresentation2.SuppressLOD = 0
      DataRepresentation2.CubeAxesXGridLines = 0
      DataRepresentation2.BlockVisibility = []
      DataRepresentation2.CubeAxesYAxisTickVisibility = 1
      DataRepresentation2.Position = [0.0, 0.0, 0.0]
      DataRepresentation2.BackfaceRepresentation = 'Follow Frontface'
      DataRepresentation2.SelectionOpacity = 1.0
      DataRepresentation2.SelectionPointLabelShadow = 0
      DataRepresentation2.CubeAxesYGridLines = 0
      DataRepresentation2.CubeAxesZAxisTickVisibility = 1
      DataRepresentation2.OrientationMode = 'Direction'
      DataRepresentation2.ScaleMode = 'No Data Scaling Off'
      DataRepresentation2.Diffuse = 1.0
      DataRepresentation2.SelectionUseOutline = 0
      DataRepresentation2.CubeAxesZTitle = 'Z-Axis'
      DataRepresentation2.Specular = 0.1
      DataRepresentation2.SelectionVisibility = 1
      DataRepresentation2.InterpolateScalarsBeforeMapping = 1
      DataRepresentation2.CustomRangeActive = [0, 0, 0]
      DataRepresentation2.Origin = [0.0, 0.0, 0.0]
      DataRepresentation2.CubeAxesVisibility = 0
      DataRepresentation2.Scale = [1.0, 1.0, 1.0]
      DataRepresentation2.SelectionCellLabelJustification = 'Left'
      DataRepresentation2.DiffuseColor = [1.0, 1.0, 1.0]
      DataRepresentation2.Shade = 0
      DataRepresentation2.SelectionCellLabelOpacity = 1.0
      DataRepresentation2.CubeAxesInertia = 1
      DataRepresentation2.Source = []
      DataRepresentation2.Masking = 0
      DataRepresentation2.Opacity = 1.0
      DataRepresentation2.LineWidth = 1.0
      DataRepresentation2.MeshVisibility = 0
      DataRepresentation2.Visibility = 1
      DataRepresentation2.SelectionCellLabelFontSize = 18
      DataRepresentation2.CubeAxesCornerOffset = 0.0
      DataRepresentation2.SelectionPointLabelJustification = 'Left'
      DataRepresentation2.OriginalBoundsRangeActive = [0, 0, 0]
      DataRepresentation2.SelectionPointLabelVisibility = 0
      DataRepresentation2.SelectOrientationVectors = ''
      DataRepresentation2.CubeAxesTickLocation = 'Inside'
      DataRepresentation2.BackfaceDiffuseColor = [1.0, 1.0, 1.0]
      DataRepresentation2.CubeAxesYLabelFormat = '%-#6.3g'
      DataRepresentation2.CubeAxesYAxisVisibility = 1
      DataRepresentation2.SelectionPointLabelFontFamily = 'Arial'
      DataRepresentation2.CubeAxesUseDefaultYTitle = 1
      DataRepresentation2.SelectScaleArray = ''
      DataRepresentation2.CubeAxesYTitle = 'Y-Axis'
      DataRepresentation2.ColorAttributeType = 'POINT_DATA'
      DataRepresentation2.AxesOrigin = [0.0, 0.0, 0.0]
      DataRepresentation2.UserTransform = [1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0]
      DataRepresentation2.SpecularPower = 100.0
      DataRepresentation2.Texture = []
      DataRepresentation2.SelectionCellLabelShadow = 0
      DataRepresentation2.AmbientColor = [1.0, 1.0, 1.0]
      DataRepresentation2.BlockOpacity = {}
      DataRepresentation2.MapScalars = 1
      DataRepresentation2.PointSize = 2.0
      DataRepresentation2.CubeAxesUseDefaultXTitle = 1
      DataRepresentation2.SelectionCellLabelFormat = ''
      DataRepresentation2.Scaling = 0
      DataRepresentation2.StaticMode = 0
      DataRepresentation2.SelectionCellLabelColor = [0.0, 1.0, 0.0]
      DataRepresentation2.SliceMode = 'XY Plane'
      DataRepresentation2.EdgeColor = [0.0, 0.0, 0.5000076295109483]
      DataRepresentation2.CubeAxesXAxisTickVisibility = 1
      DataRepresentation2.SelectionCellLabelVisibility = 0
      DataRepresentation2.NonlinearSubdivisionLevel = 1
      DataRepresentation2.CubeAxesColor = [1.0, 1.0, 1.0]
      DataRepresentation2.Representation = 'Volume'
      DataRepresentation2.CustomRange = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
      DataRepresentation2.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
      DataRepresentation2.Orientation = [0.0, 0.0, 0.0]
      DataRepresentation2.CubeAxesXTitle = 'X-Axis'
      DataRepresentation2.ScalarOpacityUnitDistance = 0.027063293868263713
      DataRepresentation2.BackfaceOpacity = 1.0
      DataRepresentation2.SelectionPointLabelFontSize = 18
      DataRepresentation2.SelectionCellFieldDataArrayName = 'u'
      DataRepresentation2.SelectionColor = [1.0, 0.0, 1.0]
      DataRepresentation2.BlockColor = {}
      DataRepresentation2.Ambient = 0.0
      DataRepresentation2.VolumeRenderingMode = 'Smart'
      DataRepresentation2.CubeAxesXAxisMinorTickVisibility = 1
      DataRepresentation2.ScaleFactor = 0.1
      DataRepresentation2.BackfaceAmbientColor = [1.0, 1.0, 1.0]
      DataRepresentation2.Slice = 0
      DataRepresentation2.ScalarOpacityFunction = a1_u_PiecewiseFunction
      DataRepresentation2.SelectMaskArray = ''
      DataRepresentation2.SelectionLineWidth = 2.0
      DataRepresentation2.CubeAxesZAxisMinorTickVisibility = 1
      DataRepresentation2.CubeAxesXAxisVisibility = 1
      DataRepresentation2.CubeAxesXLabelFormat = '%-#6.3g'
      DataRepresentation2.Interpolation = 'Gouraud'
      DataRepresentation2.CubeAxesZLabelFormat = '%-#6.3g'
      DataRepresentation2.SelectionCellLabelFontFamily = 'Arial'
      DataRepresentation2.SelectionCellLabelItalic = 0
      DataRepresentation2.CubeAxesYAxisMinorTickVisibility = 1
      DataRepresentation2.InterpolationType = 'Cubic'
      DataRepresentation2.CubeAxesZGridLines = 0
      DataRepresentation2.SelectionPointLabelFormat = ''
      DataRepresentation2.SelectionPointLabelOpacity = 1.0
      DataRepresentation2.UseAxesOrigin = 0
      DataRepresentation2.CubeAxesFlyMode = 'Closest Triad'
      DataRepresentation2.Pickable = 1
      DataRepresentation2.CustomBoundsActive = [0, 0, 0]
      DataRepresentation2.CubeAxesGridLineLocation = 'All Faces'
      DataRepresentation2.SelectionRepresentation = 'Wireframe'
      DataRepresentation2.SelectionPointLabelBold = 0
      DataRepresentation2.ColorArrayName = ('POINT_DATA', 'u')
      DataRepresentation2.SelectionPointLabelItalic = 0
      DataRepresentation2.AllowSpecularHighlightingWithScalarColoring = 0
      DataRepresentation2.SpecularColor = [1.0, 1.0, 1.0]
      DataRepresentation2.CubeAxesUseDefaultZTitle = 1
      DataRepresentation2.LookupTable = a1_u_PVLookupTable
      DataRepresentation2.SelectionPointSize = 5.0
      DataRepresentation2.SelectionCellLabelBold = 0
      DataRepresentation2.Orient = 0
      
    return Pipeline()

  class CoProcessor(coprocessing.CoProcessor):
    def CreatePipeline(self, datadescription):
      self.Pipeline = _CreatePipeline(self, datadescription)

  coprocessor = CoProcessor()
  freqs = {'input': [1]}
  coprocessor.SetUpdateFrequencies(freqs)
  return coprocessor

#--------------------------------------------------------------
# Global variables that will hold the pipeline for each timestep
# Creating the CoProcessor object, doesn't actually create the ParaView pipeline.
# It will be automatically setup when coprocessor.UpdateProducers() is called the
# first time.
coprocessor = CreateCoProcessor()

#--------------------------------------------------------------
# Enable Live-Visualizaton with ParaView
coprocessor.EnableLiveVisualization(False)


# ---------------------- Data Selection method ----------------------

def RequestDataDescription(datadescription):
    "Callback to populate the request for current timestep"
    global coprocessor
    if datadescription.GetForceOutput() == True:
        # We are just going to request all fields and meshes from the simulation
        # code/adaptor.
        for i in range(datadescription.GetNumberOfInputDescriptions()):
            datadescription.GetInputDescription(i).AllFieldsOn()
            datadescription.GetInputDescription(i).GenerateMeshOn()
        return

    # setup requests for all inputs based on the requirements of the
    # pipeline.
    coprocessor.LoadRequestedData(datadescription)

# ------------------------ Processing method ------------------------

def DoCoProcessing(datadescription):
    "Callback to do co-processing for current timestep"
    global coprocessor

    # Update the coprocessor by providing it the newly generated simulation data.
    # If the pipeline hasn't been setup yet, this will setup the pipeline.
    coprocessor.UpdateProducers(datadescription)

    # Write output data, if appropriate.
    coprocessor.WriteData(datadescription);

    # Write image capture (Last arg: rescale lookup table), if appropriate.
    coprocessor.WriteImages(datadescription, rescale_lookuptable=False)

    # Live Visualization, if enabled.
    coprocessor.DoLiveVisualization(datadescription, "localhost", 22222)
