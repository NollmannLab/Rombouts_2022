
proc do_TIFF_fileList_from_folder { rootDir log_filename destDir tcl_files } {

source "${tcl_files}log_file.tcl"
set dirs1 [glob -type d "${rootDir}*"]
set fileList {}

log_file $log_filename "List of TIFF files to be processed in $rootDir" "w"

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

	puts $ifile
	set filetail [file tail $ifile]
	set fileList [split $filetail "_"]
	set RT [string range [lindex $fileList 3] 2 end]
	set run [lindex $fileList 1]
	set roi [lindex $fileList 4]
	set fileAttributes($counter.RT) $RT
	set fileAttributes($counter.roi) $roi
	set fileAttributes($counter.run) $run
	
	incr counter
	}

return [array get fileAttributes]

}
