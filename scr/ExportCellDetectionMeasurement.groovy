// This script allows to detect the cells and write GeoJSON annotations
//
// Script author: Jianhong Ou
// Script should work with QuPath v0.3.0 or newer
// This source code is licensed under the MIT license
// Warning: this script only works for under QuPath GUI

import qupath.lib.gui.dialogs.Dialogs
import qupath.lib.scripting.QP
import qupath.lib.gui.scripting.QPEx

// define the output file prefix
// the measurements will be saved as prefix+geojson and prefix+tsv
// please make sure the parent folder is created!

var server = getCurrentServer()
def imagePath = GeneralTools.toPath(server.getURIs()[0])
def pPath = imagePath.getParent().toString()+"/" + File.createTempFile("measurementsExport", "").getName().toString() + "/"
File directory = new File(pPath)
if(! directory.exists()) {
   directory.mkdir()
}
def imageName = imagePath.getFileName().toString()
def outputPath = pPath + imageName
println('Measurements will export to '+pPath)

// need to add annotations to hierarchy so qupath sees them
def hierarchy = QP.getCurrentHierarchy()

//*********Get GeoJSON automatically based on naming scheme
def path = outputPath+".geojson"

// need to add annotations to hierarchy so qupath sees them
def cellObjects = getCellObjects()

// 'FEATURE_COLLECTION' is standard GeoJSON format for multiple objects
exportObjectsToGeoJson(cellObjects, path, "FEATURE_COLLECTION")

def csvpath = outputPath+".tsv"
def header_keys = cellObjects[0].getMeasurements().keySet()
def header = ['Cell.ID', 'Cell.classification', 'Cell.X', 'Cell.Y'] + header_keys
def outcsv = new File(csvpath)

outcsv.withWriter {
   out ->
       out.writeLine header.join('\t')
       cellObjects.forEach { cell ->
           def m = cell.getMeasurements()
           def mVal = []
           header_keys.forEach { key ->
               if(m[key]) {
                   mVal.add(m[key])
               }else {
                  mVal.add(0) 
               }
           }
           def val = [
               cell.getID(),
               cell.getClassifications()[0],
               cell.getROI().getCentroidX(),
               cell.getROI().getCentroidY()] + mVal
           out.writeLine val.join('\t')
       }
}

println('Export cell information Done!')
