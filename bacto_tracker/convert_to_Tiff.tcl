# This Huygens script converts RAMM images to stacked tiff

proc convert2Tiff {srcImg destImg chan} {
	set dims [$srcImg getdims]
	# z dimensions of the input image
	set zdim [lindex $dims 2]
        # new Channel dimension
        set newZ [expr $zdim / $chan]

	set type [$srcImg getconfig -mode type]
	set span $dims 

	# for-loop that creates the new channels
		for { set ch 0 } { $ch < $chan } { incr ch } {
		    lset span 2 $newZ
		    # create temporary image for current channel
		    set destImg:Ch$ch [img create $destImg:Ch$ch \
		                        -dim $span -type $type]
		    # for-loop that builds the new z-stack in current channel
		     for { set z 0 } { $z < $newZ } { incr z } {
		        lset span 2 1
		        $srcImg cp -> $destImg:Ch$ch -span $span \
		            -srco [list 0 0 [expr $ch + $chan * $z]] \
		            -desto [list 0 0 $z]
		    }
		    # join the current channel with previous channels (only needed for ch>0)
		    if { $ch > 0 } {
		        $destImg:Ch0 join $destImg:Ch$ch -> $destImg:Ch0
		        # delete current channel, as it is not needed anymore
		        $destImg:Ch$ch del
		    }
		}

	# copy the (joined) channel image into the the destination image
	$destImg:Ch0 -> $destImg

	# allow this if you want a thumb image to be shown... it may slow down conversion.	

	# delete the (joined) channel image, as it is not needed anymore
	$destImg:Ch0 del

	
}



