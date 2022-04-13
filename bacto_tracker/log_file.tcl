#sends input string to file

proc log_file {filename text writeMode} {



		# open the filename for writing	
		set fileId [open $filename $writeMode]

		# send the data to the file -
		#  failure to add '-nonewline' will result in an extra newline
		# at the end of the file
		puts $fileId "$text"

		# close the file, ensuring the data is written out before you continue
		#  with processing.
		close $fileId


}

