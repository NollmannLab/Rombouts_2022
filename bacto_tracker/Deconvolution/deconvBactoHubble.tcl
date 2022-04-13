
proc decon_segment {tcl_files rootDir destDir} {

#source "${tcl_files}merfish_parameters.tcl"
source "${tcl_files}convert_to_Tiff.tcl"
source "${tcl_files}log_file.tcl"
source "${tcl_files}parse_analyse_output.tcl"
source "${tcl_files}do_TIFF_fileList_from_folder.tcl"
source "${tcl_files}Utilities.tcl"


# sets parameters
set dx 0.105#0.098
set dy $dx
set dz 0.25
set NA 1.2#1.3#1.45
# medium refractive index (water: 1.338)
set ri 1.334#1.406 
# lens refractive index (water immersion: 1.338, oil immersion: 1.515)
#set ril 1.515
set ril 1.334#1.406
set singleDir true

set wavelength_ex 561
set wavelength_em 600

#set wavelength_ex 488
#set wavelength_em 520
	
set srcType   "tiff"
set destType  "tiff16"

set TIME_start [clock clicks -milliseconds]
set log_filename "${destDir}output_${TIME_start}.log"
set writeMode "a"
set thresholds_filename "${destDir}thresholds_${TIME_start}.log"

# Retrieve the list of files stored in the source folder.
if {$singleDir == "true"} { 
	set fileList [glob "${rootDir}*.tif"]
} else {
	set fileList [do_TIFF_fileList_from_folder "$rootDir" "$log_filename" \
	"$destDir" "$tcl_files"]
	set fileAttributes0 [retrieve_file_attributes $fileList]
	array set fileAttributes $fileAttributes0
}

# Loop over the images located in the source folder.
set Nimages 0	
set totalImages [llength $fileList]

foreach imgInDisk $fileList {

	# Open this image.
	set rawImg [img open $imgInDisk -foreignTo float -cmode scale -series tool]

#	huPrint report "> File [expr $Nimages + 1] out of $totalImages loaded"	
	log_file $log_filename "---------------------------------------------- \
		\n > File [expr $Nimages + 1] out of $totalImages loaded" $writeMode


	# Sets parameters
	$rawImg setp -dx $dx -dy $dy -dz $dz -na $NA -ri $ri -ril $ril -ex $wavelength_ex -em $wavelength_em
#	for {set ch 0} {$ch < $Nchannels} {incr ch} {
#		$rawImg setp -chan $ch -ex $wavelength(ch$ch,ex) -chan $ch -em $wavelength(ch$ch,em)
#	}
	
	# Get the microscopic parameters of this raw image.
#	array set rawImgParams [$rawImg setp -tclReturn]

	# Create an empty image for the deconvolved version of the image.
	set destImg [$rawImg repl "deconvolved"]

	# Deconvolve the image with default settings. 
	# Add arguments to change the SNR, nr. of iterations, etc. 
	$rawImg cmle psf -> $destImg -bgMode wf -bgRadius 0.7 -blMode off 

	# Save the deconvolved image.
#	$rawImg save "${destDir}${rawImg}_converted" -type $destType -tiffMultiDir
	$destImg save "${destDir}${rawImg}_decon" -type $destType -tiffMultiDir
	
#	huPrint info "> Deconvolved Image saved to: ${destDir}${rawImg}_decon"	
	log_file $log_filename "> Deconvolved Image saved \
		to: ${destDir}${rawImg}_decon" $writeMode

	log_file $log_filename "> Segmented Images and objects saved " $writeMode

	set TIME_taken [expr [clock clicks -milliseconds] - $TIME_start]
	set TIME_human [hh:mm:ss [expr $TIME_taken / 1000 ] ]
	log_file $log_filename "> Cumulative time: $TIME_human" $writeMode

	# Continue to the next iteration of the loop: another raw image.
	incr Nimages

	#clears up memory
	$destImg del
	$rawImg del
	
}

# set total processing time
set TIME_taken [expr [clock clicks -milliseconds] - $TIME_start]
set TIME_human [hh:mm:ss [expr $TIME_taken / 1000 ] ]

#huPrint info ">number of images: $Nimages" 
#huPrint info ">time taken $TIME_human" 
	
# Reset verbosity level to old value
huOpt verb -mode $oldVerbosity

}





