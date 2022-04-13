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


# set l { 1000 200 15 11 9 8 8 8}
# set diff_thres 1
 
 set i 0
 set r ""
 
 foreach iList [lrange $N_obj 1 end] {
 	append r "[expr abs($iList - [lindex $N_obj $i]) ] "
 	incr i
 	}
 
 set ithresh_max 0	
 foreach iList $r {
 #	echo ">>$iList, $ithresh_max"
 	if { $iList < $diff_thres } {
 #		echo "out of here: $iList < $diff_thres"
 		break
 		} else {
 			incr ithresh_max
 			set threshold_index [lindex $N_obj $ithresh_max]
 		}
 }
 		
 
 return $ithresh_max
# echo $threshold_index
 
}

# huPrint info "Utilities.tcl read"


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
#		set N [lindex [split [$newImg range] " "] 1]
#		huPrint info "channel $ch, garbage: [lindex $garbage $ch]"
		
		append Nobj " { thresh [lindex $thresholds $ch] ch $ch Nobj $N }"
		$thresImage_0 del		
	}

	$newImg del
	$thresholded_image del
	
	return $Nobj
}

############################################################################

#function level = otsu(histogramCounts)
#total = sum(histogramCounts); % '''total''' is the number of pixels in the given image. 
#%% OTSU automatic thresholding method
#sumB = 0;
#wB = 0;
#maximum = 0.0;
#sum1 = dot( (0:255), histogramCounts); 

#for ii=1:256
#    wB = wB + histogramCounts(ii);
#    wF = total - wB;
#    if (wB == 0 || wF == 0)
#        continue;
#    end
#    sumB = sumB +  (ii-1) * histogramCounts(ii);
#    mF = (sum1 - sumB) / wF;
#    between = wB * wF * ((sumB / wB) - mF) * ((sumB / wB) - mF);
#    if ( between >= maximum )
#        level = ii;
##        maximum = between;
  #  end
#end
#end

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

	#	set ptr [open /home/marcnol/Dropbox/projects/methodological/SVI/Images/test/otsu.dat "w"]
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
#		puts $ptr "ii=$iBin hist=$data_array($iBin,0,0,0) otsu=$otsu_level ; max=$maximum; bet=$between WF=$WF WB=$WB sumB=$sumB mf=$mF"			
	}
	
#	close $ptr
		
	
	
	return $otsu_level
	
}



#proc histogram { destImg ch Nbins} {
	
	
	
	
#}
















