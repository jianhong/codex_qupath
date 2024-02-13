// This script allows to detect the cells and merge multinucleis annotations
//
// Script author: Jianhong Ou
// Script should work with QuPath v0.3.0 or newer
// This source code is licensed under the MIT license
// Warning: this script only works for under QuPath GUI

import qupath.ext.stardist.StarDist2D
import qupath.ext.biop.cellpose.Cellpose2D
import qupath.lib.objects.PathObjects
import qupath.lib.analysis.features.ObjectMeasurements
import qupath.lib.images.ImageData
import qupath.lib.images.servers.ConcatChannelsImageServer
import qupath.lib.images.servers.TransformedServerBuilder
import qupath.lib.scripting.QP
import qupath.lib.roi.RoiTools
import qupath.lib.roi.GeometryTools
import qupath.lib.regions.ImagePlane

// resample the image, do not change
var downsample = 1.0

def expandPixels= 200
// change the modelPath by the filepath of your stardist model
stardistModel = '../models/stardist/stardist_cell_seg_model.pb'
// cellpose model
def cellPoseModel = 'cyto'
// the cutoff distance of nucleis for multiple nucleis cells
def distanceCutoff = 2.0
// the measurementId to be used as a marker of multiple nucleis cells
// the cutoff percentage of measurement for the measurementId, 0.9 mean 90%, if it is no less than 1, then it will be absolute cutoff value
def measurementPercentileCutoff = ["cd68: Cytoplasm: Median":180, "hist3h2a: Membrane: Mean":50]
def cytoChannels = ["CD4", "cd68"]
def nucleiChannel = "DAPI-01"
// the cellCutoff will be used to filter the cells, leave it as [:] to skip the filter
// the cutoff percentage of measurement for the measurementId, 0.9 mean 90%, if it is no less than 1, then it will be absolute cutoff value
def cellCutoff = ["DAPI-01: Nucleus: Median":10]
def newDetection = true // true to remove old detections

println 'cytoChannels:'
println cytoChannels
println 'cell filters:'
println cellCutoff
println 'multiple nucleis cells measurement cutoff:'
println measurementPercentileCutoff
println 'cell distance cutoff:' + distanceCutoff
println 'cellPoseModel:' + cellPoseModel

if (newDetection) {
    clearDetections()
}

// Set some variables
var imageData = getCurrentImageData()
var server = getCurrentServer()
var pathObjects = QP.getSelectedObjects()
if (pathObjects.isEmpty()) {
    //createSelectAllObject(true)
    throw new Exception("can not get selected object.")
}
var cal = server.getPixelCalibration()
// Extract the cytoplasm channel, rescale, and mix to a new channel
//def ops = [
//    ImageOps.percentile(1, 99)
//    ]
//def op = ImageOps.buildImageDataOp().appendOps(*ops)
//def normServer = ImageOps.buildServer(imageData, op, server.getPixelCalibration())
def cytoplasmChannel = new TransformedServerBuilder( getCurrentServer() )
                .extractChannels(*cytoChannels)
                .averageChannelProject().build()

// Extract the one nuclear channel
def nuclearChannel = new TransformedServerBuilder( getCurrentServer() )
                 .extractChannels(nucleiChannel).build()

// Make a combined server. Notice the order here is DAPI first, then the average
def combined = new ConcatChannelsImageServer( getCurrentServer(), [nuclearChannel, cytoplasmChannel] )

//def channelNames = combined.getMetadata().getChannels().collect { c -> c.name }
//println 'available channels: ' + channelNames

//setChannelNames("DAPI", "cytoplasm")
//channelNames = combined.getMetadata().getChannels().collect { c -> c.name }
//println 'available channels: ' + channelNames

// Need to create a new in-place ImageData for cellpose later
def newImageData = new ImageData(combined)

// run Cellpose. Careful of the channels names
def cellpose = Cellpose2D.builder( cellPoseModel )
        .pixelSize( .37 )             // Resolution for detection in um
        .channels( "Average channels", nucleiChannel )	      // Select detection channel(s)
        .normalizePercentiles(1,99)
        .tileSize(512)                  // If your GPU can take it, make larger tiles to process fewer of them. 
        .cellposeChannels(1,2)         // Need these, otherwise it just sends the one channel. These will be sent directly to --chan and --chan2
        .diameter(0.0)                    // Median object diameter. Set to 0.0 for the `bact_omni` model or for automatic computation
//        .cellExpansion(5.0)              // Approximate cells based upon nucleus expansion
//        .measureShape()                // Add shape measurements
//        .measureIntensity()             // Add cell measurements (in all compartments)  
        .simplify(0)                   // Simplification 1.6 by default, set to 0 to get the cellpose masks as precisely as possible
        .build()

// Note here that it is the imageData we created above and not the result of getCurrentImageData()
cellpose.detectObjects(newImageData, pathObjects)

// Save cytoplasm
def cytoplasms = getDetectionObjects()
println("Cellpose Done!")

// Run stardist
def stardist = StarDist2D
        .builder(stardistModel)
        .channels(nucleiChannel)            // Extract channel called 'DAPI'
        .normalizePercentiles(1, 99) // Percentile normalization
        .threshold(0.5)              // Probability (detection) threshold
        .pixelSize(0.37)              // Resolution for detection
        .cellExpansion(5)            // Expand nuclei to approximate cell boundaries
//        .measureShape()              // Add shape measurements
//        .measureIntensity()          // Add cell measurements (in all compartments)
        .build()
stardist.detectObjects(imageData, pathObjects)

println("stardist Done!")

//clearSelectedObjects()

// Save stardist detections
def sdDetection = getDetectionObjects()
//println sdDetection
//println sdDetection.size()

// assign cell names
sdDetection.eachWithIndex{ cell, x ->
    cell.setName('stardistCell'+x.toString())
}

println("Start to merge the output of cellpose and stardist.")
// Create cells
cells = []
def toRemove = []
def plane = ImagePlane.getDefaultPlane()
def hierarchy = getCurrentHierarchy()

cytoplasms.each{ cytoplasm ->
    // expand from the current cell about defined pixel
    def cell1Geom = cytoplasm.getROI().getGeometry()
    def cell1Expansion = cell1Geom.buffer(expandPixels)
    def cell1ExpROI = GeometryTools.geometryToROI(cell1Expansion, plane)
    // get all cells in the expansion region
    cell2candidates = hierarchy.getObjectsForROI(null, cell1ExpROI).findAll{ it.isCell() && it.getName().startsWith('stardistCell')}
    cell2candidates.each{ nucleus ->      
        if ( cytoplasm.getROI().contains( nucleus.getROI().getCentroidX() , nucleus.getROI().getCentroidY())){
            newcell = PathObjects.createCellObject(
                RoiTools.combineROIs(cytoplasm.getROI(), nucleus.getNucleusROI(), RoiTools.CombineOp.ADD),
                nucleus.getNucleusROI(), getPathClass("SingleNucleiFromCellPose"), null )
            cells.add(newcell)
            sdDetection = sdDetection - nucleus
            toRemove.add(nucleus)
            }
        }
    }
 // add the undetected cells that can not be detected by Cellpose
 //println sdDetection.size()
sdDetection.each{ nucleus ->
    toRemove.add(nucleus)
    nucleus.setPathClass(getPathClass("SingleNucleiFormStarDist"))
    cells.add(nucleus)
 }

// remove the cell annotation
removeObjects(toRemove,true)
hierarchy.addObjects(cells)
// Add measurements
def measurements = ObjectMeasurements.Measurements.values() as List
def compartments = ObjectMeasurements.Compartments.values() as List
def shape = ObjectMeasurements.ShapeFeatures.values() as List
//def cells = getCellObjects()
for ( cell in cells ) {
    ObjectMeasurements.addIntensityMeasurements( server, cell, downsample, measurements, compartments )
    ObjectMeasurements.addCellShapeMeasurements( cell, cal,  shape )
    }
    
// def cells = getCellObjects()

// println cells[0].measurements.get(measurementId)
// reassign cell names
cells.eachWithIndex{ cell, x ->
    cell.setName('Cell'+x.toString())
}
// println measurement
// get the third quartile
def Quartiles(float pos, List<Double> val) {
    int length = val.size()
    double quartile
    int quartilePos = (int)(length * pos) - 1
    Arrays.sort(val)
    if (quartilePos % 1 == 0) {
        quartile = val[quartilePos]
    } else {
        quartile = (val[quartilePos] + val[quartilePos + 1]) / 2
    }
    return quartile
}
// check average measurement signals in cell cytoplam
def measurement_avgs = [:]
for(ele in measurementPercentileCutoff) {
    if(ele.value<1){
        def measurement = []
        for(cell in cells){
            def measurement_mean = cell.measurements.get(ele.key)
            if(!measurement_mean) {
               measurement_mean = 0 
            }
            measurement.add(measurement_mean)
        }
        //def measurement_avg = measurement.sum()/measurement.size()
        measurement_avg = Quartiles(ele.value, measurement)
        println "Cutoff "+ ele.key + "is: "+measurement_avg.toString()
    }else {
        measurement_avg = ele.value
    }
    measurement_avgs[ele.key] = measurement_avg
}

// two check: 1, if the measurement > measurement_avg; 2, the distance of the cells is 0
def multinucleis = [:]
def distancePixels_colection = []
def cell1candidates = []
def checkCellMeasurement(cell, measurement_avgs) {
    for(ele in measurement_avgs) {
        def cell_measurement = cell.measurements.get(ele.key)
        if (!cell_measurement) {
            return false 
        }
        if (cell_measurement<ele.value) {
            return false
        }
    }
    return true
}
for (cell in cells) {
    if(checkCellMeasurement(cell, measurement_avgs)) {
       cell1candidates.add(cell) 
    }
}
println "multinucleis cell candidates are "+cell1candidates.size()
//pSize=0
for (int i=0; i<cell1candidates.size(); i++) {
    def cell1 = cell1candidates[i]
    def cell1name = cell1.getName()
    def mnname = ''
    if(multinucleis[cell1name]) {
         mnname = multinucleis[cell1name]
    }
    // expand from the current cell about 1000 pixel
    def cell1Geom = cell1.getROI().getGeometry()
    def cell1Expansion = cell1Geom.buffer(expandPixels)//.difference(cell1Geom) this will not work for cellpose results
    def cell1ExpROI = GeometryTools.geometryToROI(cell1Expansion, plane)
    // get all cells in the expansion region
    cell2candidates = hierarchy.getObjectsForROI(null, cell1ExpROI).findAll{ it.isCell()}
    //if(pSize<5) {
    //    pSize++
    //    println 'cell2 candidate size is ' + cell2candidates.size()
    //}
    for (int j=0; j<cell2candidates.size(); j++) {
        def cell2 = cell2candidates[j]
        if(cell1.getName()==cell2.getName()) continue
        if(checkCellMeasurement(cell2, measurement_avgs)){
            def g1 = cell1.getNucleusROI().getGeometry()
            def g2 = cell2.getNucleusROI().getGeometry()
            double distancePixels = g1.distance(g2)
            distancePixels_colection.add(distancePixels)
            if(distancePixels<distanceCutoff) {
                def cell2name = cell2.getName()
                def mnname2 = ''
                //check if the cell2 contains multinucleis name
                if(multinucleis[cell2name]) {
                    mnname2 = multinucleis[cell2name]
                }
                if(!mnname) {
                    mnname = mnname2
                } else {
                    if(mnname2) {
                        if(mnname!=mnname2) {
                            // rename all mnname2 to mnname
                            multinucleis.findAll{it.value==mnname2}.each{multinucleis[it.key] = mnname}
                        }
                    }
                 }
                 if(!mnname) {
                     mnname = 'multinucleis' + (multinucleis.size()+1).toString()
                 }
                 multinucleis[cell1name] = mnname
                 multinucleis[cell2name] = mnname
            }
        }
    }
}
//println distancePixels_colection.min()
//println distancePixels_colection.max()
//println distancePixels_colection.sum()/distancePixels_colection.size()
//println multinucleis
def multinucleis_names = multinucleis.values().unique(false)
//println multinucleis_names
def multinucleiCells = []
toRemove = []
// println cells.size()
for(mnname in multinucleis_names) {
    def current_group_cells = []
    multinucleis.each { k, v ->
        if(v==mnname) {
            def this_cell = cells.find{it.getName()==k}
            current_group_cells.add(this_cell)
            // remove the cell
            cells = cells - this_cell
            toRemove.add(this_cell)
        }
    }
    // merge the ROI
    //println current_group_cells[0].dump()
    def cytoplasm = current_group_cells[0].getROI()
    def nucleus = current_group_cells[0].getNucleusROI()
    for(int i=1; i<current_group_cells.size(); i++) {
        try{
            cytoplasm = RoiTools.combineROIs(cytoplasm, current_group_cells[i].getROI(), RoiTools.CombineOp.ADD)
            nucleus   = RoiTools.combineROIs(nucleus, current_group_cells[i].getNucleusROI(), RoiTools.CombineOp.ADD)
        }catch(Exception e) {
            println "found non-noded intersection: "+e
        }
    }
    // create new cells by the new cytoplasm and nucleus
    def multinucleis_cell = PathObjects.createCellObject(cytoplasm, nucleus, getPathClass("MultiNucleis"), null )
    multinucleiCells.add(multinucleis_cell)
}
// println multinucleiCells
// println cells.size()

addObjects(multinucleiCells)

// remove the cell annotation
removeObjects(toRemove,true)

// Add measurements
cells = getCellObjects()
for ( cell in cells ) {
    ObjectMeasurements.addIntensityMeasurements( server, cell, downsample, measurements, compartments )
    ObjectMeasurements.addCellShapeMeasurements( cell, cal,  shape )
    }
println("Measurements Done!")

// filter cells
println "Filter cells by"
toRemove = []
measurement_avgs = [:]
for(ele in cellCutoff) {
    if(ele.value<1){
        def measurement = []
        for(cell in cells){
            def measurement_mean = cell.measurements.get(ele.key)
            if(!measurement_mean) {
               measurement_mean = 0 
            }
            measurement.add(measurement_mean)
        }
        //def measurement_avg = measurement.sum()/measurement.size()
        measurement_avg = Quartiles(ele.value, measurement)
        println "Cutoff "+ ele.key + "is: "+measurement_avg.toString()
    }else {
        measurement_avg = ele.value
    }
    measurement_avgs[ele.key] = measurement_avg
}
println measurement_avgs
for ( cell in cells ) {
    if(!checkCellMeasurement(cell, measurement_avgs)) {
        //println cell.measurements.get(measurement_avgs.keySet()[0])
            // remove the cell
            toRemove.add(cell)
       }
}
println("Removed "+toRemove.size()+" cells")
removeObjects(toRemove,true)

// Finished!
fireHierarchyUpdate()
println("Done!")
