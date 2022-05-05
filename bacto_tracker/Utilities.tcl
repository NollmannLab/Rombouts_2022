#utilities

proc hh:mm:ss {secs} {
        set h [expr {$secs/3600}]
        incr secs [expr {$h*-3600}]
        set m [expr {$secs/60}]
        set s [expr {$secs%60}]
        format "%02.2d:%02.2d:%02.2d" $h $m $s
}

############################################################################
# from a List of objectes N_obj, it finds the diff(N_obj) then finds what
# vector index has the abs(diff(N_obj))<diff_thres
# and returns this index


proc find_best_threshold { N_obj diff_thres} { 

 set i 0
 set r ""
 
 foreach iList [lrange $N_obj 1 end] {
 	append r "[expr abs($iList - [lindex $N_obj $i]) ] "
 	incr i
 	}
 
 set ithresh_max 0	
 foreach iList $r {
 	if { $iList < $diff_thres } {
 		break
 		} else {
 			incr ithresh_max
 			set threshold_index [lindex $N_obj $ithresh_max]
 		}
 }
 		
 
 return $ithresh_max
 
}

############################################################################

proc thresholds_labels { destImg thresholds garbage dims border_fill type} {
	
	set thresholded_image [$destImg repl "thresholded"] 

	$destImg thres -> $thresholded_image $thresholds -mode all
	set newImg [img create newImg -dim $dims -type $type]

	for {set ch 0} {$ch < [llength $thresholds]} {incr ch} {
		set thresImage_0 [$thresholded_image split -mode one -chan $ch]
		$thresImage_0 borderFill $border_fill 

		$thresImage_0 label -> $newImg

		set segm_output [ $thresImage_0 analyze $newImg \
			 -garb [lindex $garbage $ch] -style tcl ] 
			 
		set N [llength $segm_output]
		
		append Nobj " { thresh [lindex $thresholds $ch] ch $ch Nobj $N }"
		$thresImage_0 del		
	}

	$newImg del
	$thresholded_image del
	
	return $Nobj
}

############################################################################

##
proc otsu { destImg ch } {

	
	set hist_data [$destImg hist -chan $ch -exportData]
	set data [lindex $hist_data 11 ]
	array set data_array $data
	set Nbins [array size data_array]
	set Npixels [expr [lindex [$destImg getdims] 0] * [lindex [$destImg getdims] 1 ] * [lindex [$destImg getdims] 2]]

	set sumB 0
	set WB 0
	set maximum 0.0
	set sum1 0
	set total 0
	
	for {set iBin 0} {$iBin < $Nbins} {incr iBin} {
		incr sum1 [expr $iBin * $data_array($iBin,0,0,0)]
		incr total $data_array($iBin,0,0,0)

	}

	for {set iBin 0} {$iBin < $Nbins} {incr iBin} {
		incr WB $data_array($iBin,0,0,0)
		set WF [expr $total - $WB]
		
		if {$WB == 0 | $WF ==0} {
			continue
		}
	
		incr sumB [expr ($iBin ) * $data_array($iBin,0,0,0)]
		set mF [expr ($sum1 - $sumB) / $WF ]
		set between [expr $WF * $WB * (($sumB / $WB) - $mF) * (($sumB / $WB) - $mF) ]
		
		if {$between >= $maximum} {
			set otsu_level $iBin
			set maximum $between
		}
	}
	
	return $otsu_level
	
}


















