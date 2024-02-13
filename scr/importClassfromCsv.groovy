// This script allows to detect the cells and merge multinucleis annotations
//
// Script author: Jianhong Ou
// Script should work with QuPath v0.4.0 or newer
// This source code is licensed under the MIT license
// Warning: this script only works for under QuPath GUI

import qupath.lib.scripting.QP
import qupath.lib.objects.classes.PathClass
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

// Find the column indices for 'Object.ID' and 'Class'
def objectIdIndex = -1
def classIndex = -1

def headerRow = csvData[0]
headerRow.eachWithIndex { value, index ->
    if (value == 'Object.ID') {
        objectIdIndex = index
    } else if (value == 'Class') {
        classIndex = index
    }
}

println objectIdIndex
println classIndex
// Create the 'Object ID' to 'Class' map
def objectIdToClassMap = [:]

csvData.eachWithIndex { row, rowIndex ->
    if (rowIndex > 0) {
        def objectId = row[objectIdIndex]
        def objectClass = row[classIndex]
        objectIdToClassMap[objectId] = objectClass
    }
}
println objectIdToClassMap
// Print the 'Object ID' to 'Class' map

def pathClasses = getQuPath().getAvailablePathClasses()
print(pathClasses)
for (def detection : QP.getDetectionObjects()) {
    //def measurementList = detection.getMeasurementList()
    def allClass = detection.getClassifications()
    //println(allClass)
    objID = detection.getID().toString()
    
    def newClassID = objectIdToClassMap[objID]
    def newClass = PathClass.fromString(newClassID)
    if (!pathClasses.contains(newClass)){
        pathClasses.add(newClass)}
    //detection.setPathClass(newClass)
    // println(newClassID)
    if(!allClass.contains(newClassID)) {
        allClass = allClass + newClassID
        //println(allClass)
        detection.setClassifications(allClass)
    }
}

// Fire update event
fireHierarchyUpdate()

println("Done")