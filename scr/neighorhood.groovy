// This script allows to run neighborhood analysis
//
// Script author: Jianhong Ou
// Script should work with QuPath v0.3.0 or newer
// This source code is licensed under the MIT license
// Warning: this script only works for under QuPath GUI

import qupath.ext.biop.cellpose.CellposeSetup

// change the path to neighborhood_ana.py
def neighborhood_ana_script = 'neighborhood_ana.py'
// change the program waiting time, default is 6000s
def waitTime = 6000000
// print message or not
def silent = false

CellposeSetup cellposeSetup = CellposeSetup.getInstance();

// Make sure that cellposeSetup.getCellposePythonPath() is not empty
if (cellposeSetup.getCellposePytonPath().isEmpty()) {
    throw new IllegalStateException("Cellpose python path is empty. Please set it in Edit > Preferences");
}
// fill the command
neighborhood_ana_script = cellposeSetup.getCellposePytonPath() + neighborhood_ana_script

def task = neighborhood_ana_script.execute()
def result = [std: '', err: '']
def ready = false
Thread.start {
   def reader = new BufferedReader(new InputStreamReader(task.in))
   def line = ""
   while ((line = reader.readLine()) != null) {
       if (silent != false) {
           println "" + line
       }
       result.std += line + "\n"
   }
   ready = true
   reader.close()
}
task.waitForOrKill(waitTime)
def error = task.err.text

if (error.isEmpty()) {
   return result
} else {
   throw new RuntimeException("\n" + error)
}
