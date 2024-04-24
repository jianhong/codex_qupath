// This script allows to fill the detections
//
// Script author: Jianhong Ou
// Script should work with QuPath v0.3.0 or newer
// This source code is licensed under the MIT license
// Warning: this script only works for under QuPath GUI

// set cell type here, has priorities, first come first assign
// The key is the measurement key [marker: position: stats_method]
def classifier = [
    'cd68: Cytoplasm: Median':['cutoff':50, 'celltype':'M1 Macrophage'],
    'CD4: Cytoplasm: Median':['cutoff':50, 'celltype':'Helper T cells'],
    'sox9: Nucleus: Max':['cutoff':50, 'celltype':'Cartilage tumor cell'],
    'ki67: Nucleus: Max':['cutoff':50, 'celltype':'Proliferation Marker'],
    'lcp1: Cytoplasm: Max':['cutoff':50, 'celltype':'Tumor/Myeloid'],
    'cd206: Cytoplasm: Median':['cutoff':50, 'celltype':'M2 MACROPHAGE'],
    'CD163: Cytoplasm: Median':['cutoff':50, 'celltype':'M2 MACROPHAGE'],
    'cd45: Cytoplasm: Median':['cutoff':30, 'celltype':'HSC'],
    'CD14: Cytoplasm: Median':['cutoff':30, 'celltype':'MONOCYTE'],
    'cd8: Cytoplasm: Median':['cutoff':50, 'celltype':'Cytotoxic T Cells'],
    'CD11c: Cytoplasm: Median':['cutoff':50, 'celltype':'mDC/cDC - Dendritic cells'],
    'cd19: Cytoplasm: Median':['cutoff':50, 'celltype':'B CELLS'],
    'cd34: Cytoplasm: Median':['cutoff':30, 'celltype':'Endothelial'],
    'asma: Cytoplasm: Max':['cutoff':50, 'celltype':'Pericyte/Fibroblast']
]

def cells = getCellObjects()
def detections = cells.collect {
    PathObjects.createAnnotationObject(it.getROI(), it.getPathClass(), it.getMeasurementList())
}

addObjects(detections)

def nucleus = cells.collect {
    PathObjects.createAnnotationObject(it.getNucleusROI(), getPathClass('nucleus'), it.getMeasurementList())
}

addObjects(nucleus)

def celltype = cells.collect {
    def cell_classification = 'unknown'
    for(ele in classifier) {
        if(it.measurements.get(ele.key) > ele.value.cutoff) {
            cell_classification = ele.value.celltype
            continue
        }
    }
    PathObjects.createDetectionObject(it.getROI(), getPathClass(cell_classification), it.getMeasurementList())
}

addObjects(celltype)

println "Done!"
