
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

    a2DRenderView1 = CreateView( Create2DRenderView, "image_%t.png", 1, 0, 1, cp_views )
     
    # mvm: PV bug workaround, client doesn't export viewsize
    # mvm: set to Lock View Size Custom settings
    a2DRenderView1.ViewSize = [640,480]

    a2DRenderView1.LightSpecularColor = [1.0, 1.0, 1.0]
    a2DRenderView1.KeyLightAzimuth = 10.0
    a2DRenderView1.UseTexturedBackground = 0
    a2DRenderView1.UseLight = 0
    a2DRenderView1.CameraPosition = [316.21445673906658, 513.7761642761908, 2038.7425740279175]
    a2DRenderView1.FillLightKFRatio = 3.0
    a2DRenderView1.Background2 = [0.0, 0.0, 0.16500000000000001]
    a2DRenderView1.FillLightAzimuth = -10.0
    a2DRenderView1.LODResolution = 50.0
    a2DRenderView1.BackgroundTexture = []
    a2DRenderView1.InteractionMode = '3D'
    a2DRenderView1.StencilCapable = 1
    a2DRenderView1.LightIntensity = 1.0
    a2DRenderView1.CameraFocalPoint = [316.21445673906658, 513.7761642761908, -756.14540211504584]
    a2DRenderView1.ImageReductionFactor = 2
    a2DRenderView1.CameraViewAngle = 30.0
    a2DRenderView1.CameraParallelScale = 723.37023715383816
    a2DRenderView1.EyeAngle = 2.0
    a2DRenderView1.HeadLightKHRatio = 3.0
    a2DRenderView1.StereoRender = 0
    a2DRenderView1.KeyLightIntensity = 0.75
    a2DRenderView1.BackLightAzimuth = 110.0
    a2DRenderView1.AxesVisibility = 0
    a2DRenderView1.OrientationAxesInteractivity = 0
    a2DRenderView1.UseInteractiveRenderingForSceenshots = 0
    a2DRenderView1.UseOffscreenRendering = 0
    a2DRenderView1.Background = [0.18431372549019609, 0.18431372549019609, 0.18431372549019609]
    a2DRenderView1.UseOffscreenRenderingForScreenshots = 0
    a2DRenderView1.NonInteractiveRenderDelay = 2
    a2DRenderView1.CenterOfRotation = [0.0, 0.0, 0.0]
    a2DRenderView1.CameraParallelProjection = 0
    a2DRenderView1.CompressorConfig = 'vtkSquirtCompressor 0 3'
    a2DRenderView1.HeadLightWarmth = 0.5
    a2DRenderView1.MaximumNumberOfPeels = 4
    a2DRenderView1.LightDiffuseColor = [1.0, 1.0, 1.0]
    a2DRenderView1.StereoType = 'Red-Blue'
    a2DRenderView1.DepthPeeling = 1
    a2DRenderView1.BackLightKBRatio = 3.5
    a2DRenderView1.StereoCapableWindow = 1
    a2DRenderView1.CameraViewUp = [0.0, 1.0, 0.0]
    a2DRenderView1.LightType = 'HeadLight'
    a2DRenderView1.LightAmbientColor = [1.0, 1.0, 1.0]
    a2DRenderView1.RemoteRenderThreshold = 3.0
    a2DRenderView1.CacheKey = 0.0
    a2DRenderView1.UseCache = 0
    a2DRenderView1.KeyLightElevation = 50.0
    a2DRenderView1.CenterAxesVisibility = 1
    a2DRenderView1.MaintainLuminance = 0
    a2DRenderView1.StillRenderImageReductionFactor = 1
    a2DRenderView1.BackLightWarmth = 0.5
    a2DRenderView1.FillLightElevation = -75.0
    a2DRenderView1.MultiSamples = 0
    a2DRenderView1.FillLightWarmth = 0.40000000000000002
    a2DRenderView1.AlphaBitPlanes = 1
    a2DRenderView1.LightSwitch = 1
    a2DRenderView1.OrientationAxesVisibility = 1
    a2DRenderView1.CameraClippingRange = [2018.3551482876383, 2069.3237126383365]
    a2DRenderView1.BackLightElevation = 0.0
    a2DRenderView1.ViewTime = 0.0
    a2DRenderView1.OrientationAxesOutlineColor = [1.0, 1.0, 1.0]
    a2DRenderView1.LODThreshold = 5.0
    a2DRenderView1.CollectGeometryThreshold = 100.0
    a2DRenderView1.UseGradientBackground = 0
    a2DRenderView1.KeyLightWarmth = 0.59999999999999998
    a2DRenderView1.OrientationAxesLabelColor = [1.0, 1.0, 1.0]
    
    ns2dcn = CreateProducer( datadescription, "input" )
    
    a1_omeg_PiecewiseFunction = CreatePiecewiseFunction( Points=[0.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0] )
    
    a1_omeg_PVLookupTable = GetLookupTableForArray( "omeg", 1, Discretize=1, RGBPoints=[-13166394.436894085, 0.0, 1.0, 1.0, -1316639.4436894078, 0.0, 0.0, 0.99215686274509807, 0.0, 0.0, 0.0, 0.52156862745098043, 1316639.4436894096, 0.98039215686274506, 0.0, 0.011764705882352941, 13166394.436894085, 1.0, 1.0, 0.0], UseLogScale=0, VectorComponent=0, NanColor=[1.0, 1.0, 0.0], NumberOfTableValues=256, ColorSpace='RGB', VectorMode='Magnitude', HSVWrap=0, ScalarRangeInitialized=1.0, LockScalarRange=0 )
    
    # mvm: PV bug, Enabled=1 means 'active gui component' which draws a box around it
    # mvm: set to 0
    ScalarBarWidgetRepresentation1 = CreateScalarBar( Title='omeg', Position2=[0.12999999999999995, 0.96555323590814335], TitleOpacity=1.0, TitleShadow=0, AutomaticLabelFormat=0, TitleFontSize=14, TitleColor=[1.0, 1.0, 1.0], AspectRatio=20.0, NumberOfLabels=5, ComponentTitle='', Resizable=1, TitleFontFamily='Arial', Visibility=1, LabelFontSize=12, LabelFontFamily='Arial', TitleItalic=0, Selectable=0, LabelItalic=0, Enabled=0, LabelColor=[1.0, 1.0, 1.0], Position=[0.076572769953051609, 0.029227557411273364], LabelBold=0, UseNonCompositedRenderer=1, LabelOpacity=1.0, TitleBold=0, LabelFormat='%-#6.3g', Orientation='Vertical', LabelShadow=0, LookupTable=a1_omeg_PVLookupTable, Repositionable=1 )
    GetRenderView().Representations.append(ScalarBarWidgetRepresentation1)
    
    DataRepresentation1 = Show()
    DataRepresentation1.Opacity = 1.0
    DataRepresentation1.Origin = [0.0, 0.0, 0.0]
    DataRepresentation1.Scale = [1.0, 1.0, 1.0]
    DataRepresentation1.UseXYPlane = 0
    DataRepresentation1.ColorAttributeType = 'POINT_DATA'
    DataRepresentation1.Position = [0.0, 0.0, 0.0]
    DataRepresentation1.ColorArrayName = 'omeg'
    DataRepresentation1.Visibility = 1
    DataRepresentation1.Slice = 0
    DataRepresentation1.LookupTable = a1_omeg_PVLookupTable
    DataRepresentation1.MapScalars = 1
    DataRepresentation1.SliceMode = 'XY Plane'
    DataRepresentation1.Pickable = 1
    DataRepresentation1.Orientation = [0.0, 0.0, 0.0]
    

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

