

proc manual_thresholding1 { img treshold_value0 garbage_value0} { 

#set treshold_value0 100

set tresholded_image [$img repl "thresholded"]
set labeled_image [$img repl "labeled"]


$img thres -> $tresholded_image {$treshold_value0}
$tresholded_image label -> $labeled_image -garb $garbage_value0

$tresholded_image show
$labeled_image show


}


proc manual_thresholding2 { } { 
set imgInDisk "/mnt/PALM_dataserv/DATA/Andres-Delphine SIM/Squeeze_Embryo/20160802_Microarray_early+late_new_acquisitions/20160802_R508_squeeze_early_embryo_MicroArray_004_SIR_ALX.dv"
	
set treshold_value0 10000
set treshold_value1 2500

set garbage_value0 10
set garbage_value1 10

set img0 [img open $imgInDisk -foreignTo float -cmode scale]

set tresholded_image [$img0 repl "thresholded" ]
set labeled_image [$img0 repl "labeled"]

$img0 thres -> $tresholded_image "$treshold_value0 $treshold_value1" -mode all

#set img0_ch1_gauss [$img0 repl "gauss"]
#$img0_ch1 kuwahara -> $img0_ch1_gauss -width 1 1 1 3 0

set channel 1
$tresholded_image label -> $labeled_image -chan $channel
$tresholded_image analyze $labeled_image -chan $channel -style tcl

#-garb $garbage_value0
#$labeled_image show
#$tresholded_image show
#$img0 show


}


