#!/usr/bin/env pvpython

import sys
from paraview.simple import *

def renderPorpoiseTracks(domainFile, trackFile, outputFile):
    #### disable automatic camera reset on 'Show'
    paraview.simple._DisableFirstRenderCameraReset()

    domain = XMLUnstructuredGridReader(FileName=[domainFile])
    tracks = PVDReader(FileName=trackFile)

    # Set up animation
    animationScene1 = GetAnimationScene()
    animationScene1.UpdateAnimationUsingDataTimeSteps()

    # get active view
    renderView1 = GetActiveViewOrCreate('RenderView')

    # show data in view
    domainDisplay = Show(domain, renderView1)
    domainDisplay.ScalarOpacityUnitDistance = 8718.547699279605

    # reset view to fit data
    renderView1.ViewSize=[1200,800]
    renderView1.ResetCamera()

    # set active source
    SetActiveSource(domain)

    # set scalar coloring
    ColorBy(domainDisplay, ('POINTS', 'Depth'))

    # get color transfer function/color map for 'Depth'
    depthLUT = GetColorTransferFunction('Depth')
    depthLUT.RGBPoints = [0.0, 0.968627, 0.984314, 1.0, 0.01, 0.701427, 0.826932, 0.910994, 3.0, 0.62269, 0.793469, 0.883436, 4.0, 0.52314, 0.739193, 0.86154, 5.0, 0.422751, 0.684077, 0.839887, 8.0, 0.341421, 0.628962, 0.808698, 10.0, 0.260716, 0.573846, 0.777203, 20.0, 0.195392, 0.509117, 0.743786, 30.0, 0.130434, 0.44416, 0.710323, 40.0, 0.0809644, 0.381125, 0.661357, 50.0, 0.031754, 0.318135, 0.612146, 79.9, 0.0313725, 0.253193, 0.51606, 80.0, 0.0299992, 0.250004, 0.500008]
    depthLUT.ColorSpace = 'Lab'
    depthLUT.NanColor = [0.498039, 0.0, 0.0]
    depthLUT.ScalarRangeInitialized = 1.0

    domainDisplay.SetScalarBarVisibility(renderView1, False)

    tracksDisplay = Show(tracks, renderView1)

    # create a new 'Temporal Particles To Pathlines'
    temporalParticlesToPathlines1 = TemporalParticlesToPathlines(Input=tracks, Selection=None)
    temporalParticlesToPathlines1.MaskPoints = 1
    temporalParticlesToPathlines1.MaxTrackLength = 250000
    temporalParticlesToPathlines1.MaxStepDistance = [10000.0, 10000.0, 100.0]

    # hide data in view
    Hide(tracks, renderView1)

    # show data in view
    pathlines = Show(servermanager.OutputPort(temporalParticlesToPathlines1, 1), renderView1)
    # set active source
    SetActiveSource(temporalParticlesToPathlines1)

    pathlines = Show(temporalParticlesToPathlines1, renderView1)
    pathlines.Position = [0.0, 0.0, 10.0]
    ColorBy(pathlines, ('POINTS', 'Internal Status'))

    # get color transfer function/color map for 'InternalStatus'
    internalStatusLUT = GetColorTransferFunction('InternalStatus')
    internalStatusLUT.RGBPoints = [0.0, 0.6, 0.6, 0.6, 0.5, 0.6, 0.6, 0.6, 0.501, 0.0, 1.0, 0.0, 1.5, 0.0, 1.0, 0.0, 1.501, 0.0, 1.0, 0.0, 1.501, 1.0, 0.0, 0.0, 2.0, 1.0, 0.0, 0.0]
    internalStatusLUT.ColorSpace = 'HSV'
    internalStatusLUT.NanColor = [0.498039, 0.0, 0.0]
    internalStatusLUT.ScalarRangeInitialized = 1.0

    # get opacity transfer function/opacity map for 'InternalStatus'
    internalStatusPWF = GetOpacityTransferFunction('InternalStatus')
    internalStatusPWF.Points = [1.0, 0.0, 0.5, 0.0, 2.0, 1.0, 0.5, 0.0]
    internalStatusPWF.ScalarRangeInitialized = 1

    pathlines.ColorArrayName = ['POINTS', 'Internal Status']
    pathlines.SetScalarBarVisibility(renderView1, False)
    pathlines.LookupTable = internalStatusLUT

    SetActiveSource(domain)
    # reset view to fit data
    renderView1.ResetCamera()

    # reset view to fit data bounds
    renderView1.ResetCamera(282742.5, 770655.375, 5638386.0, 6006389.5, -77.9334335327, 10.0)
    animationScene1.GoToFirst()

    animationScene1.Play()

    # current camera placement for renderView1
    renderView1.CameraPosition = [471713.9270489631, 5820931.193431761, 683628.2383665707]
    renderView1.CameraFocalPoint = [471713.9270489631, 5820931.193431761, -33.96671676635]
    renderView1.CameraParallelScale = 305567.3237337823

    # save screenshot
    SaveScreenshot(outputFile, magnification=2, quality=100, view=renderView1)

    Delete(pathlines)
    Delete(temporalParticlesToPathlines1)
    Delete(tracksDisplay)
    Delete(tracks)
    Delete(domainDisplay)
    Delete(domain)
    Delete(depthLUT)
    Delete(internalStatusLUT)
    Delete(renderView1)
    Delete(animationScene1)


if __name__ == "__main__":
    import os
    import argparse

    cmdparser = argparse.ArgumentParser(description="Render Porpoise tracks using Paraview")
    cmdparser.add_argument("domain", help="VTU domain file")
    cmdparser.add_argument("pvdfile", help="PVD file referring to track data", nargs='+')
    cmdparser.add_argument("-k", "--keep-together", action="store_true", default=False, dest="keeptogether",
            help="Save plots adjacent to original PVD files.")

    args = cmdparser.parse_args()

    for track in args.pvdfile:
        if not args.keeptogether:
            outfile = os.path.normpath(track).replace(os.sep, "_")

        outfile = os.path.splitext(outfile)[0]
        print "Rendering tracks for %s..." % track
        renderPorpoiseTracks(args.domain, track, "%s.tracks.png" % outfile)
