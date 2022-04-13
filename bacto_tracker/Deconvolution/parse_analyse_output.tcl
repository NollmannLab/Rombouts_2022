

proc parse_analyse_output { segm_output filename } {
#set segm_output [ $tresholded_image analyze $labeled_image \
#			 -garb $garbage_value	-chan $channel -style tcl ] 

set fileId [open $filename "w"]

set val_length [llength $segm_output]

set file_header "label, vol, surf, sVol, min, max, sum, mean, var, len, wiLat, \
	wiAx, CMx, CMy, CMz, bboxOffsetx, bboxOffsety, bboxOffsetz, bboxSpanx, bboxSpany, bboxSpanz"
puts $fileId $file_header

for {set i 0} {$i < $val_length} {incr i} {
	set ioutput [lindex $segm_output $i]
	array set one_object $ioutput

	set label $one_object(label)
	set vol $one_object(vol)
	set surf $one_object(surf)
	set sVol $one_object(sVol)
	set min $one_object(min)
	set max $one_object(max)
	set sum $one_object(sum)
	set mean $one_object(mean)
	set var $one_object(var)
	set len $one_object(len)
	set wiLat $one_object(wiLat)
	set wiAx $one_object(wiAx)
	set CM $one_object(CM)
	set bboxOffset $one_object(O)
	set bboxSpan $one_object(S)

#	huPrint info "label $label, CM $CM"
	set output_string "$label, $vol, $surf, $sVol, $min, $max, $sum, $mean,\
			$var, $len, $wiLat, $wiAx"
 	foreach iCM $CM {
		append output_string ", $iCM"
	}
	foreach ibboxOffset $bboxOffset {
		append output_string ", $ibboxOffset"
	}
	foreach ibboxSpan $bboxSpan {
		append output_string ", $ibboxSpan"
	}
	puts $fileId $output_string
}

close $fileId

}
