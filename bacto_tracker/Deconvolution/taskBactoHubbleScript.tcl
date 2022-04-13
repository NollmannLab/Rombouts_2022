# batchDeconvolve.tcl
set tcl_files "/mnt/PALM_dataserv/DATA/Sara/DATA/huygensDeconvolutionScripts/"
set dataDir "/mnt/grey/DATA/rawData_2020/Experiment_6_sara/tiff_raw_data/Ch_1/"
set destDir "/mnt/grey/DATA/rawData_2020/Experiment_6_sara/Ecoli_deconvolved/"

source "${tcl_files}/deconvBactoHubble.tcl"

decon_segment "$tcl_files" "$dataDir" "$destDir"

