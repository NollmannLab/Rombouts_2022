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

huPrint info "Utilities.tcl read"

