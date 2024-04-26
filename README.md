# codex_qupath
The script collection for `QuPath`.

## Installation

You don't need to install the script.
While some scripts can be executed as commands,
most of them must be run within the `QuPath` GUI.
The cell segmentation requires the [cellpose](https://github.com/BIOP/qupath-extension-cellpose) and [stardist](https://github.com/qupath/qupath-extension-stardist) extension in `QuPath`.

```
git clone git@github.com:jianhong/codex_qupath.git
## list available scripts
ls codex_qupath/scr
## list available stardist models
ls codex_qupath/models/stardist
```
## Usage

### Step 1. Cell segmentation

This step involves merging the results from [Cellpose](https://www.cellpose.org/) and [Stardist](https://github.com/stardist/stardist) analyses.

First, you'll need to install the [Cellpose](https://github.com/BIOP/qupath-extension-cellpose) and [Stardist](https://github.com/qupath/qupath-extension-stardist) extensions. Don't forget to set the Python path for the Cellpose extension.

Next, navigate to line 25 of the `cellpose_stardist_multiNucleis_underGUI.groovy` script and update the path of the `stardistModel` variable to match the file path of your `cell_seg_model.pb`. Adjust the parameters from [line 26 to line 38](https://github.com/jianhong/codex_qupath/blob/main/scr/cellpose_stardist_multiNucleis_underGUI.groovy#L26-L38) accordingly.  You may want to change the `cytoChannels` and `measurementPercentileCutoff` according to your observations. The `distanceCutoff` is the maximal cutoff distance of nucleus for multiple nucleus cells.

Then, select a region in the opened TIFF view window in QuPath. Start with a small region when testing the script, and expand it once you achieve the desired results.

Finally, run the `cellpose_stardist_multiNucleis_underGUI.groovy` script by clicking the Run button in the Script Editor.

[![Youtube Video](https://img.youtube.com/vi/3SKKSDhlEkw/0.jpg)](https://youtu.be/3SKKSDhlEkw)

### Step 2. Export the cell segmentation

The script `ExportCellDetectionMeasurement.groovy` is designed to export cell segmentation data. Upon export, two files are generated and saved in a folder prefixed with "measurementsExport". 

- The `.tsv` file is compatible with `createSeuratObj.R`.
- The `.geojson` file is intended for reloading into `QuPath`. 

These exported files contain comprehensive information, including cell area, locations, marker signal statistics, and nucleus classification. The cell locations are particularly useful for neighborhood analysis.

### Step 3. Fill the cell with colors

The script `fill_detections.groovy` is employed to duplicate detected cells into `Detections` and `Annotations` objects within `QuPath`. Subsequently, user can establishe the classes in QuPath to assign colors.

Users may want to try different classifier for the cell type detection by reset the [line10-25](https://github.com/jianhong/codex_qupath/blob/main/scr/fill_detections.groovy#L10-L25). Here I set an example for the classifier.

[![Youtube Video](https://img.youtube.com/vi/4F6sQWKDFl4/0.jpg)](https://youtu.be/4F6sQWKDFl4)

## Credits

The scripts were originally a collection from published paper, [Image.sc Forum](https://forum.image.sc/), 
and [imagescientist.com](https://www.imagescientist.com/) or written by [Jianhong Ou](https://github.com/jianhong).

This work were supported by [Duke Regeneromics](https://sites.duke.edu/dukeregenerationcenter/regeneromics/) and [Visgauss Lab @Duke.edu](https://ortho.duke.edu/visgauss-laboratory).

## Contributions and Support

If you would like to contribute to this collection, please [create a pull request](https://github.com/jianhong/codex_qupath/compare).

For further information or help, don't hesitate to get in touch on [issue channel](https://github.com/jianhong/codex_qupath/issues/new/choose).

## Citations

If you use jianhong/codex_qupath for your analysis, please cite it using the following url: [https://github.com/jianhong/codex_qupath](https://github.com/jianhong/codex_qupath).

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.


