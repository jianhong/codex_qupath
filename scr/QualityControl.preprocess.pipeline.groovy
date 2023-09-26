// This script allows to detect the cells and write GeoJSON annotations from a folder
//
// Script author: Jianhong Ou
// Script should work with QuPath v0.3.0 or newer
// This source code is licensed under the MIT license
// set args by command line
//  env=/admin/apps/community/regeneromics/conda_envs/qupath
//  conda activate $env
//  QuPath.sh script $env/scripts/QualityControl.preprocess.pipeline.groovy \
//      --save \
//      --args imagefolder=path/to/qptiff/folder \
//      --args model=stardist_cell_seg_model \
//      --args channel='DAPI-01'
//
//
// arguments: imagefolder:The qptiff file path. The project files will be
//                        created in the subfolder 'QuPathProject'.
//            model:      The path of cell segmentation model. The available 
//                        models are stardist_cell_segmodel, dsb2018_paper,
//                        dsb2018_heavy_augment, he_heavy_augment. Full path
//                        of model file are also acceptable, eg /path/to/model.pb
//                        Default is stardist_cell_seg_model
//            channel:    The DAPI channel name. This is used for StarDist2D 
//                        to define the nucleus. Default is 'DAPI-01'
//            tissueThreshold:
//                        The background signal threshold. Integer. Default is 3.
//                        There are issues for the Fluorescence image, 
//                        it can only handle first 3 channels.
//                        It may not work well
//
//  step 2
//  qptiff.tsv.boxplot.R path/to/qptiff/folder 50
//  here 50 is the y limitation of the plots.


import groovy.io.FileType
import java.awt.image.BufferedImage
import qupath.lib.images.servers.ImageServerProvider
import qupath.lib.gui.commands.ProjectCommands
import qupath.lib.gui.tools.GuiTools
import qupath.lib.gui.images.stores.DefaultImageRegionStore
import qupath.lib.gui.images.stores.ImageRegionStoreFactory
import qupath.lib.gui.dialogs.Dialogs
import qupath.lib.scripting.QP
import qupath.lib.gui.scripting.QPEx
import qupath.ext.stardist.StarDist2D

// define the default model path
def envs = '/admin/apps/community/regeneromics/conda_envs/qupath'
def starDist2D_model_path = envs + '/models/stardist'

params = [
    'model'           : 'stardist_cell_seg_model',
    'channel'         : 'DAPI-01',
    'tissueThreshold' : 3
]
for(arg in args){
    String[] val = arg.split("=")
    params[val[0]] = val[1]
}

File modelFile = new File(params.model)
if (!modelFile.exists()){
    modelFile = new File(starDist2D_model_path + File.separator + params.model + '.pb')
}
def modelPath = modelFile.toString()
if(!modelFile.exists()){
    throw new Exception("Can not locate the model at " + modelPath)
}

def useChannel = params.channel
// tissueThreshold for tissue detection.
def tissueThreshold  = params.tissueThreshold.toInteger()

if (!params.imagefolder) {
    throw new Exception("parameter 'imagefolder' is required!")
    
}

// Check if the image directory in there...
File imgDirectory = new File(params.imagefolder)
if(! imgDirectory.exists()){
    throw new Exception(params.imagefolder + ' does not exists!')
}

def projectName = "QuPathProject"
File projDirectory = new File(imgDirectory.toString() + File.separator + projectName)

if (!projDirectory.exists()){
    projDirectory.mkdirs()
}


println 'Project will be saved to ' + projDirectory.toString()
println 'the model path is ' + modelPath
println 'select channel is ' + useChannel


// Define the function to detect tissue.
// the tissue will be defined by imagej plugin 'SimpleTissueDetection2'
// The parameters are set to Fluorescence and darkBackground to true
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
        ',  "maxHoleAreaMicrons": 100000000000000.0' +
        ',  "darkBackground": true' +
        ',  "smoothImage": true' +
        ',  "medianCleanup": true' +
        ',  "dilateBoundaries": false' +
        ',  "smoothCoordinates": true' +
        ',  "excludeOnBoundary": false' +
        ',  "singleAnnotation": false}');
    
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
    return(imageData)
}

def cell_segmentation(imageData, modelPath, useChannel)
{
    // Check if the channel is in the images.
    def channelNames = imageData.getServer().getMetadata().getChannels().collect { c -> c.name }
    if (!channelNames.contains(useChannel)){
        throw new Exception('available channels: ' + channelNames.join(', ')) 
    }
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
    //pathObjects = QP.getAnnotationObjects() //QP.getSelectedObjects()
    def hierarchy=imageData.getHierarchy()
    def pathObjects=hierarchy.getAnnotationObjects()
    if (pathObjects.isEmpty()) {
        throw new Exception("No parent objects are selected!")
    }
    stardist.detectObjects(imageData, pathObjects)
    stardist.close() // This can help clean up & regain memory
    println('Cell segmentation Done!')
    return(imageData)
}

// export the signals
def export_signals(imageData){
    def server = imageData.getServer()

    // need to add annotations to hierarchy so qupath sees them
    def hierarchy = imageData.getHierarchy()
    
    //*********Get GeoJSON automatically based on naming scheme 
    def path = GeneralTools.toPath(server.getURIs()[0]).toString()+".geojson"
    
    // need to add annotations to hierarchy so qupath sees them
    def cellObjects = hierarchy.getCellObjects()
    
    // 'FEATURE_COLLECTION' is standard GeoJSON format for multiple objects
    exportObjectsToGeoJson(cellObjects, path, "FEATURE_COLLECTION")
    
    // define the tsv file path
    def csvpath = GeneralTools.toPath(server.getURIs()[0]).toString()+".tsv"
    // get table columns
    def header = cellObjects[0].getMeasurements().keySet()
    // save to tsv file
    def outcsv = new File(csvpath)
    outcsv.withWriter {
       out ->
           out.writeLine cellObjects[0].getMeasurements().keySet().join('\t')
           cellObjects.forEach {
               out.writeLine it.getMeasurements().values().join('\t')
               }
    }
    
    println('Export cell information Done!')
}


// Create project
def project = Projects.createProject(projDirectory , BufferedImage.class)

// Set up cache
def imageRegionStore = ImageRegionStoreFactory.createImageRegionStore(QuPathGUI.getTileCacheSizeBytes());

// Some filetypes are split between a name and a folder and we need to eliminate the folder from our recursive search.
// This is the case for vsi files for instance.
def skipList = []
imgDirectory.eachFileRecurse (FileType.FILES) { file ->
    if (file.name.endsWith(".vsi")) {
        print(file.name)
        f = new File(file.parent+File.separator+"_"+file.name.substring(0, file.name.length() - 4)+"_")
        skipList.add(f.toString()) //getCanonicalPath())
        return
    }
}

// Add files to the project
imgDirectory.eachFileRecurse (FileType.FILES) { file ->
    def imagePath = file.getCanonicalPath()
    skip = false
    for (p in skipList) {
        //print("--->"+p)
        if (imagePath.startsWith(p)) {
            skip = true
        }
        
    }
    if (skip == true) {
        //print("Skipping "+imagePath)
        return
    }
        
    // Skip a folder if there is a corresponding .vsi file.
    if (file.isDirectory()) {
        print(file.getParent())
        print(file.getName().startsWith('_') && file.getName().endsWith('_'))
        return
    }
    
    // Skip the project directory itself
    if (file.getCanonicalPath().startsWith(projDirectory.getCanonicalPath() + File.separator))
        return
        
    // I tend to add underscores to the end of filenames I want excluded
    // MacOSX seems to add hidden files that start with a dot (._), don't add those
    if (file.getName().endsWith("_") || file.getName().startsWith("."))
        return

    // Is it a file we know how to read?
    def support = ImageServerProvider.getPreferredUriImageSupport(BufferedImage.class, imagePath)
    if (support == null)
        return

    // iterate through the scenes contained in the image file
    support.builders.eachWithIndex { builder, i -> 
        sceneName = file.getName()
        
        if (sceneName.endsWith('.vsi')) {
            //This is specific to .vsi files, we do not add a scene name to a vsi file
            if (support.builders.size() >= 3 && i < 2) {
                return;
            }
        } else {
            if (support.builders.size() > 1)
                sceneName += " - Scene #" + (i+1)
        }
        // Add a new entry for the current builder and remove it if we weren't able to read the image.
        // I don't like it but I wasn't able to use PathIO.readImageData().
        entry = project.addImage(builder)
    
        try {
            imageData = entry.readImageData()
        } catch (Exception ex) {
            println sceneName +" -- Error reading image data " + ex
            project.removeImage(entry, true)
            return
        }
        
        println "Adding: " + sceneName
    
        // Set a particular image type automatically (based on /qupath/lib/gui/QuPathGUI.java#L2847)
        def imageType = GuiTools.estimateImageType(imageData.getServer(), imageRegionStore.getThumbnail(imageData.getServer(), 0, 0, true));
        imageData.setImageType(imageType)
        println "Image type estimated to be " + imageType
    
        // Write a thumbnail if we can
        var img = ProjectCommands.getThumbnailRGB(imageData.getServer());
        entry.setThumbnail(img)
        
        // Add an entry name (the filename)
        entry.setImageName(sceneName)
        
        // start from new, clear all objects
        clearAllObjects()
        createSelectAllObject(true)
        
        // detect tissue and cells
        imageData = detect_tissue(imageData, tissueThreshold, 100)
        
        // cell segmentation
        imageData = cell_segmentation(imageData, modelPath, useChannel)
        
        // export cell signals
        export_signals(imageData)

        // Adding image data to the project entry
        entry.saveImageData(imageData)
    }
}

// Changes should now be reflected in the project directory
project.syncChanges()

println 'Done!'

