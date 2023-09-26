// This script allows to detect the cells and write GeoJSON annotations
//
// Script author: Jianhong Ou
// Script should work with QuPath v0.3.0 or newer
// This source code is licensed under the MIT license

import qupath.ext.stardist.StarDist2D
import qupath.lib.gui.dialogs.Dialogs
import qupath.lib.scripting.QP
import qupath.lib.gui.scripting.QPEx

// set args by command line
//  conda activate /admin/apps/community/regeneromics/conda_envs/qupath
//  QuPath.sh script scripts/cellDetection.pipeline.groovy -i test.qptiff --save --args models/stardist/stardist_cell_seg_model.pb --args 'DAPI-01'
def modelPath = args[0]
def useChannel = args[1]
// tissueThreshold for tissue detection. There are issues for the Fluorescence image, it can only handle first 3 channels
def tissueThreshold  = 3//args[2]//3

def channelNames = getCurrentServer().getMetadata().getChannels().collect { c -> c.name }
println 'available channels: ' + channelNames

println 'modelPath is ' + modelPath
println 'select channel is ' + useChannel

def detect_tissue(imageData, threshold=3, requestedPixelSizeMicrons=100, minAreaMicrons=1000)
{
    if(imageData.isBrightfield()){
        throw new Exception("Cannot handle brightfield trained model!")
    }
    def imageType='Fluorescence'
    QP.setImageType(imageType)
    QP.logger.info("detecting tissue {} {} {}", imageData, imageType, threshold)

    QP.runPlugin('qupath.imagej.detect.tissue.SimpleTissueDetection2',
    imageData, '{"threshold": ' + threshold +
        ',  "requestedPixelSizeMicrons": ' +
        requestedPixelSizeMicrons +
        ',  "minAreaMicrons": '+ minAreaMicrons +
        ',  "maxHoleAreaMicrons": 100000000000000.0,  "darkBackground": true,  "smoothImage": true,  "medianCleanup": true,  "dilateBoundaries": false,  "smoothCoordinates": true,  "excludeOnBoundary": false,  "singleAnnotation": false}');

    def double areaMax = 0
    def hierarchy = imageData.getHierarchy()

    for (annotation in hierarchy.getAnnotationObjects()) {
        roi = annotation.getROI()
        area = roi.getArea()
        if (area > areaMax){
            areaMax = area
        }
    }

    //Minimum area is 10% of areaMax
    areaMin = Math.round(areaMax * .1)

    def int i = 1
    for (annotation in hierarchy.getAnnotationObjects()) {
        roi = annotation.getROI()
        area = roi.getArea()
        if (area >= areaMin){
            annotation.setName("Tissue"+i)
            annotation.setPathClass(PathClassFactory.getPathClass("Tissue"))
            i += 1
        }
        else {
            hierarchy.removeObject(annotation,true)
        }
    }

    hierarchy.getSelectionModel().setSelectedObject(null);
}

// select the whole image
QP.createFullImageAnnotation(true)
def pathObjects = QP.getSelectedObjects()

// Run detection for the selected objects
def imageData = QP.getCurrentImageData()

// start from new, clear all objects
clearAllObjects()
createSelectAllObject(true)

//detect tissue and cells
detect_tissue(imageData, tissueThreshold, 100)

// Customize how the StarDist detection should be applied
// Here some reasonable default options are specified
def stardist = StarDist2D
    .builder(modelPath)
    .channels(useChannel)            // Extract channel called 'DAPI'
    .normalizePercentiles(1, 99) // Percentile normalization
    .threshold(0.5)              // Probability (detection) threshold
    .pixelSize(0.37)              // Resolution for detection
    .cellExpansion(5)            // Expand nuclei to approximate cell boundaries
    .measureShape()              // Add shape measurements
    .measureIntensity()          // Add cell measurements (in all compartments)
    .build()

// Define which objects will be used as the 'parents' for detection
// Use QP.getAnnotationObjects() if you want to use all annotations, rather than selected objects
pathObjects = QP.getAnnotationObjects() //QP.getSelectedObjects()

if (pathObjects.isEmpty()) {
    QP.getLogger().error("No parent objects are selected!")
    return
}
stardist.detectObjects(imageData, pathObjects)
stardist.close() // This can help clean up & regain memory
println('Cell segmentation Done!')


def server = QP.getCurrentImageData().getServer()

// need to add annotations to hierarchy so qupath sees them
def hierarchy = QP.getCurrentHierarchy()

//*********Get GeoJSON automatically based on naming scheme
def path = GeneralTools.toPath(server.getURIs()[0]).toString()+".geojson"

// need to add annotations to hierarchy so qupath sees them
def cellObjects = getCellObjects()

// 'FEATURE_COLLECTION' is standard GeoJSON format for multiple objects
exportObjectsToGeoJson(cellObjects, path, "FEATURE_COLLECTION")

def csvpath = GeneralTools.toPath(server.getURIs()[0]).toString()+".tsv"
def header = cellObjects[0].getMeasurements().keySet()
def outcsv = new File(csvpath)
outcsv.withWriter {
   out ->
       out.writeLine cellObjects[0].getMeasurements().keySet().join('\t')
       cellObjects.forEach {
           out.writeLine it.getMeasurements().values().join('\t')
           }
}

println('Export cell information Done!')
