// This script allows to read the GeoJSON annotations outputed by 
//   ExportCellDetectionMeasurement.groovy
//
// Script author: Jianhong Ou
// Script should work with QuPath v0.3.0 or newer
// This source code is licensed under the MIT license
// Warning: this script only works for under QuPath GUI

import qupath.lib.scripting.QP
import com.google.gson.Gson;
import com.google.gson.JsonElement;
import org.locationtech.jts.geom.Geometry;
import qupath.lib.analysis.features.ObjectMeasurements

def createROI(var jsonObject, Gson gson) {
    var geometry = gson.fromJson(jsonObject, Geometry.class);
    geometry = GeometryTools.homogenizeGeometryCollection(geometry);
    var roi = GeometryTools.geometryToROI(geometry, ImagePlane.getDefaultPlane());
    return roi;
}
def readCellObject(JsonElement element, Gson gson) {
    def cells = []
    for (feature in element.features) {
        // create a new cell object
        var jsonObject = feature.getAsJsonObject();
        def cytoplasm = createROI(jsonObject.get("geometry"), gson);
        def nucleus = createROI(jsonObject.get("nucleusGeometry"), gson);
        def prop = jsonObject.get("properties");
        // create new cells by the new cytoplasm and nucleus
        def cellclass = prop.classification.get("name").getAsString();
        def cell = PathObjects.createCellObject(cytoplasm, nucleus, getPathClass(cellclass), null )
        UUID ID = UUID.fromString(jsonObject.get("id").toString().replaceAll("\"", ""))
        cell.setID(ID)
        cells.add(cell)
    }
    return cells;
}

def read_geojson()
{
    def server = QP.getCurrentImageData().getServer()
    var cal = server.getPixelCalibration()
    var downsample = 1.0
    def measurements = ObjectMeasurements.Measurements.values() as List
    def compartments = ObjectMeasurements.Compartments.values() as List
    def shape = ObjectMeasurements.ShapeFeatures.values() as List

    // need to add annotations to hierarchy so qupath sees them
    def hierarchy = QP.getCurrentHierarchy()
        
    // locate the json file
    def JSONfile = Dialogs.promptForFile(null)
    
    def stream = new FileInputStream(JSONfile)
    var gson = GsonTools.getInstance();
    var reader = new InputStreamReader(new BufferedInputStream(stream))
    var element = gson.fromJson(reader, JsonElement.class);
    
    var cells = readCellObject(element, gson);
    // Add measurements
    for ( cell in cells ) {
        ObjectMeasurements.addIntensityMeasurements( server, cell, downsample, measurements, compartments )
        ObjectMeasurements.addCellShapeMeasurements( cell, cal,  shape )
        }
    println("Measurements Done!")
    hierarchy.addObjects(cells)
 
    fireHierarchyUpdate()
    println "Done!"
}

read_geojson()
