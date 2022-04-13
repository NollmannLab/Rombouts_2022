# batchDeconvolve.tcl
set tcl_files "/home/marcnol/Repositories/huygensDeconvolutionScripts/"

set dataDir "/mnt/grey/DATA/rawData_2020/Experiment_5/DAPI/005_merFISH_RAMM_DAPI"
set destDir "/home/marcnol/tmp/"

source "${tcl_files}merfish.tcl"

decon_segment "$tcl_files" "$dataDir" "$destDir"
