
proc do_TIFF_fileList_from_folder { rootDir log_filename destDir tcl_files } {

source "${tcl_files}log_file.tcl"
set dirs1 [glob -type d "${rootDir}*"]
set fileList {}

log_file $log_filename "List of TIFF files to be processed in $rootDir" "w"
#huPrint info "root: $rootDir"

foreach ifile $dirs1 {
	set dirs2 [glob -type d "${ifile}/*"]	
	foreach ifile2 $dirs2 {
		set dirs3 [glob -type d "${ifile2}/*"]	
		foreach ifile3 $dirs3 {
				set dirs4 [glob -type d "${ifile3}/*"]	
				foreach ifile4 $dirs4 {
								
				set TIFF_files [glob "${ifile4}/*.tif"]
				foreach ifile5 $TIFF_files {
				#append ifile4 ".tif"				
				#huPrint info "TIFF: [file tail $ifile4] "
				log_file $log_filename "Full: [file \
					rootname $ifile5]" "a"
				lappend fileList $ifile5
						}
					}
				}
			}
	
}
	return $fileList

}

proc retrieve_file_attributes { fileList } {

set counter 0

foreach ifile $fileList {

#	huPrint info $ifile
	# set ifile [lindex $fileList 1]
	puts $ifile
	set filetail [file tail $ifile]
	set fileList [split $filetail "_"]
	set RT [string range [lindex $fileList 3] 2 end]
	set run [lindex $fileList 1]
	set roi [lindex $fileList 4]

#	set run [string range $filetail 5 7]
#	set posRT [string first "RT" $filetail]
#	set posLastUnderscore [string last "_" $filetail]
#	set substring [string range $filetail $posRT [expr $posLastUnderscore - 1] ]
#	set posUnderscore [string first "_" $substring]
#
#	set RT [string range $substring 2 [expr $posUnderscore -1 ] ]
#	set roi [string range $substring [ expr $posUnderscore + 1] [string len substring] ]

	set fileAttributes($counter.RT) $RT
	set fileAttributes($counter.roi) $roi
	set fileAttributes($counter.run) $run
	
	incr counter
	}

return [array get fileAttributes]

}

#set destDir   "/home/marcnol/SVI/Images/deconvolved/"
#set log_filename "${destDir}dir_parse_output.txt"
#set rootDir /mnt/PALM_dataserv/DATA/Andres/MERFISH/20170623_MERFISH_2_Readouts/005_merFISH_merFISH_slide4_RT7_8/
#set tcl_files "/home/marcnol/SVI/Batch/"

#do_TIFF_fileList_from_folder $rootDir $log_filename $destDir $tcl_files
