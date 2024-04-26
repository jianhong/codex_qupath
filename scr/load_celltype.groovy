// This script allows to load the celltype from a csv file output from detectCellType.R
//
// Script author: Jianhong Ou
// Script should work with QuPath v0.3.0 or newer
// This source code is licensed under the MIT license
// Warning: this script only works for under QuPath GUI

import java.io.BufferedReader
import java.io.FileReader
import qupath.lib.io.*

def csvFilePath  = Dialogs.promptForFile(null)


// Read the CSV file
def csvData = []
try {
    def reader = new BufferedReader(new FileReader(csvFilePath))
    def line
    while ((line = reader.readLine()) != null) {
        csvData.add(line.split(','))
    }
    reader.close()
} catch (Exception e) {
    e.printStackTrace()
}

// Find the column indices for 'Object.ID' and 'CellType'
// default cell id is the first column, celltype is the second column
// otherwise, find it by 'Object.ID' and 'CellType'
def objectIdIndex = 0
def classIndex = 1

def headerRow = csvData[0]
headerRow.eachWithIndex { value, index ->
    if (value == 'Object.ID') {
        objectIdIndex = index
    } else if (value == 'Class') {
        classIndex = index
    }
}

//println objectIdIndex
//println classIndex
// Create the 'Object ID' to 'Class' map
def objectIdToClassMap = [:]

csvData.eachWithIndex { row, rowIndex ->
    if (rowIndex > 0) {
        def objectId = row[objectIdIndex].replaceAll("\"", "")
        def objectClass = row[classIndex].replaceAll("\"", "")
        objectIdToClassMap[objectId] = objectClass
    }
}
println 'Imported cell types'
println objectIdToClassMap
// Print the 'Object ID' to 'Class' map

def pathClasses = getQuPath().getAvailablePathClasses()

def cells = getCellObjects()
def detections = cells.collect {
    def newClass = it.getPathClass()
    if (!pathClasses.contains(newClass)){
        pathClasses.add(newClass)}
    PathObjects.createAnnotationObject(it.getROI(), newClass, it.getMeasurementList())
}

addObjects(detections)
def nucleusClass = PathClass.fromString('nucleus')
if (!pathClasses.contains(nucleusClass)){
        pathClasses.add(nucleusClass)}
def nucleus = cells.collect {
    PathObjects.createAnnotationObject(it.getNucleusROI(), nucleusClass, it.getMeasurementList())
}

addObjects(nucleus)

def celltype = cells.collect {
    def objID = it.getID().toString()
    def newClassID = objectIdToClassMap[objID]
    def newClass = PathClass.fromString(newClassID)
    if (!pathClasses.contains(newClass)){
        pathClasses.add(newClass)}

    PathObjects.createDetectionObject(it.getROI(), newClass, it.getMeasurementList())
}

addObjects(celltype)

fireHierarchyUpdate()

println "Done!"

