
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


def DoCoProcessing(datadescription):
    "Callback to do co-processing for current timestep"
    cp_writers = []
    cp_views = []
    timestep = datadescription.GetTimeStep()

    test2_ = FindSource("test2_*")
    DataRepresentation4 = GetDisplayProperties(test2_)
    Calculator2 = GetActiveSource()
    DataRepresentation5 = GetDisplayProperties(Calculator2)
    RenderView2 = CreateView( CreateRenderView, "image_%t.png", 1, 0, 1, cp_views )
    RenderView2.LightSpecularColor = [1.0, 1.0, 1.0]
    RenderView2.KeyLightAzimuth = 10.0
    RenderView2.UseTexturedBackground = 0
    RenderView2.UseLight = 1
    RenderView2.CameraPosition = [675.48569037837376, 592.93151768488212, 587.15955621025546]
    RenderView2.FillLightKFRatio = 3.0
    RenderView2.Background2 = [0.0, 0.0, 0.16500000000000001]
    RenderView2.FillLightAzimuth = -10.0
    RenderView2.LODResolution = 50.0
    RenderView2.BackgroundTexture = []
    RenderView2.InteractionMode = '3D'
    RenderView2.StencilCapable = 1
    RenderView2.LightIntensity = 1.0
    RenderView2.CameraFocalPoint = [137.63418266030334, 120.77716310145153, 122.59972972958906]
    RenderView2.ImageReductionFactor = 2
    RenderView2.CameraViewAngle = 30.0
    RenderView2.CameraParallelScale = 220.83647796503186
    RenderView2.EyeAngle = 2.0
    RenderView2.HeadLightKHRatio = 3.0
    RenderView2.StereoRender = 0
    RenderView2.KeyLightIntensity = 0.75
    RenderView2.BackLightAzimuth = 110.0
    RenderView2.OrientationAxesInteractivity = 0
    RenderView2.UseInteractiveRenderingForSceenshots = 0
    RenderView2.UseOffscreenRendering = 0
    RenderView2.Background = [0.31999694819562063, 0.34000152590218968, 0.42999923704890519]
    RenderView2.UseOffscreenRenderingForScreenshots = 0
    RenderView2.NonInteractiveRenderDelay = 2
    RenderView2.CenterOfRotation = [127.5, 127.5, 127.5]
    RenderView2.CameraParallelProjection = 0
    RenderView2.CompressorConfig = 'vtkSquirtCompressor 0 3'
    RenderView2.HeadLightWarmth = 0.40000000000000002
    RenderView2.MaximumNumberOfPeels = 4
    RenderView2.LightDiffuseColor = [1.0, 1.0, 1.0]
    RenderView2.StereoType = 'Red-Blue'
    RenderView2.DepthPeeling = 1
    RenderView2.BackLightKBRatio = 3.5
    RenderView2.StereoCapableWindow = 1
    RenderView2.CameraViewUp = [-0.39037432619867685, 0.83217528468112856, -0.39381744630070881]
    RenderView2.LightType = 'HeadLight'
    RenderView2.LightAmbientColor = [1.0, 1.0, 1.0]
    RenderView2.RemoteRenderThreshold = 3.0
    RenderView2.CacheKey = 0.0
    RenderView2.UseCache = 0
    RenderView2.KeyLightElevation = 50.0
    RenderView2.CenterAxesVisibility = 0
    RenderView2.MaintainLuminance = 1
    RenderView2.StillRenderImageReductionFactor = 1
    RenderView2.BackLightWarmth = 0.5
    RenderView2.FillLightElevation = -75.0
    RenderView2.MultiSamples = 0
    RenderView2.FillLightWarmth = 0.40000000000000002
    RenderView2.AlphaBitPlanes = 1
    RenderView2.LightSwitch = 0
    RenderView2.OrientationAxesVisibility = 1
    RenderView2.CameraClippingRange = [406.23111967460096, 1418.0051407264114]
    RenderView2.BackLightElevation = 0.0
    RenderView2.ViewTime = 0.0
    RenderView2.OrientationAxesOutlineColor = [1.0, 1.0, 1.0]
    RenderView2.LODThreshold = 18.199999999999999
    RenderView2.CollectGeometryThreshold = 100.0
    RenderView2.UseGradientBackground = 0
    RenderView2.KeyLightWarmth = 0.59999999999999998
    RenderView2.OrientationAxesLabelColor = [1.0, 1.0, 1.0]
    
    a1_realtempx_PiecewiseFunction = CreatePiecewiseFunction( Points=[-1.5000000000002331, 0.0, 0.5, 0.0, 1.5000000000002318, 1.0, 0.5, 0.0] )
    
    a1_Result_PiecewiseFunction = CreatePiecewiseFunction( Points=[2.3316233587421763e-15, 1.0, 0.5, 0.47546851634979248, 0.87359398603439331, 0.0049999980255961418, 0.69853609800338745, 0.35839599370956421, 1.7320507764816284, 0.40495997667312622, 0.5, 0.0] )
    
    a1_realtempx_PVLookupTable = GetLookupTableForArray( "realtempx", 1, Discretize=1, RGBPoints=[-1.5000000000002331, 0.23000000000000001, 0.29899999999999999, 0.754, 1.5000000000002318, 0.70599999999999996, 0.016, 0.14999999999999999], UseLogScale=0, VectorComponent=0, NanColor=[0.25, 0.0, 0.0], NumberOfTableValues=256, ColorSpace='Diverging', VectorMode='Magnitude', HSVWrap=0, ScalarRangeInitialized=1.0, LockScalarRange=0 )
    
    a1_Result_PVLookupTable = GetLookupTableForArray( "Result", 1, Discretize=1, RGBPoints=[2.3316233520252134e-15, 0.0, 1.0, 1.0, 0.77942286340605649, 0.011764705882352941, 0.0, 0.93333333333333335, 0.86602540378450699, 0.20000000000000001, 0.0, 0.40000000000000002, 0.95262794416295749, 0.89803921568627454, 0.0, 0.050980392156862744, 1.7320508075690118, 1.0, 1.0, 0.0], UseLogScale=0, VectorComponent=0, NanColor=[1.0, 1.0, 0.0], NumberOfTableValues=256, ColorSpace='RGB', VectorMode='Magnitude', HSVWrap=0, ScalarRangeInitialized=1.0, LockScalarRange=0 )
    
    DataRepresentation4.Source.ShaftResolution = 6
    DataRepresentation4.Source.TipLength = 0.34999999999999998
    DataRepresentation4.Source.ShaftRadius = 0.029999999999999999
    DataRepresentation4.Source = "Arrow"
    DataRepresentation4.Source.Invert = 0
    DataRepresentation4.Source.TipResolution = 6
    DataRepresentation4.Source.TipRadius = 0.10000000000000001
    
    DataRepresentation5.Source.ShaftResolution = 6
    DataRepresentation5.Source.TipLength = 0.34999999999999998
    DataRepresentation5.Source.ShaftRadius = 0.029999999999999999
    DataRepresentation5.Source = "Arrow"
    DataRepresentation5.Source.Invert = 0
    DataRepresentation5.Source.TipResolution = 6
    DataRepresentation5.Source.TipRadius = 0.10000000000000001
    

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

