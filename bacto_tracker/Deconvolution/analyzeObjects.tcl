# This function analyzes objects in an image
# The input is the imageName and threshold value used for segmentation
# The function shows the center-of-mass for each object in the Task Report

proc analyzeObjects {imageName threshold} {
	# adjust baseline (correcting for intensity offset): 
	$imageName adjbl -> a

	# label the image:
	set nrObjects [a labelDirect -> b $threshold $threshold ]

	# analyze the image, using the labeled image
	a analyze b -garb 100 -style tcl

	# create a variable, and store the analysis results as a tcl list
	set objectStats [a analyze b -garb 100 -style tcl]

	# each object is an array in the list. Loop through each object:
	for {set i 1} {$i < $nrObjects} {incr i} {
		array set object [lindex $objectStats $i]

		# store the center-of-mass for each object in a new variable:
		set cm $object(CM)

		# Write the results to the Task report window:
		huOpt report "CM object $i: $cm"
	}
}

