########################################################>
# Do not change
########################################################<
# Declaration of paths and file types.
#set sourceDir "/home/marcnol/PALM/marcnol/software/huygens_script/"
#set merfishDir "/mnt/grey/DATA/rawData_2020/Experiment_5/DAPI/005_merFISH_RAMM_DAPI/"
set srcType   "tiff"
set destType  "tiff16"

########################################################>
# MAIN PARAMETERS TO SET
########################################################<
#set destDir   "/home/marcnol/SVI/Images/test/"

#set singleDir true
set singleDir false

# used to define image parameters
set experiment_type "DAPI"
#set experiment_type "RT"
#set experiment_type "3channels"

# order of channels
#set channel_order($experiment_type) "0"
#set channel_order(RT) "2"
#set channel_order(beads) "1"
#set channel_order(other) "2"

# set global versus specific 
set globalSegmPars true

# method used for thresholding: fixed search or otsu. None does no segmentation!
#set thresholding_method "fixed"
#set thresholding_method "search"
#set thresholding_method "otsu"
set thresholding_method "none"

set Nchannels 3

set border_fill "50 50 1"

########################################################<
# Settings level 2
########################################################<

# pixel size, um
# emCCD front illuminated
#set dx 0.08
# sCMOS
set dx 0.105#0.098
set dy $dx
set dz 0.25
set NA 1.2#1.3#1.45
# medium refractive index (water: 1.338)
set ri 1.334#1.406 
# lens refractive index (water immersion: 1.338, oil immersion: 1.515)
#set ril 1.515
set ril 1.334#1.406

########################################################<


# set image parameters
switch $Nchannels {
	"1" {
	set wavelength(ch0,ex) 405
	set wavelength(ch0,em) 440
	}
	"2" {
		switch $experiment_type {
		"DAPI" {
			set wavelength(ch0,ex) 405
			set wavelength(ch1,ex) 561
			set wavelength(ch0,em) 440
			set wavelength(ch1,em) 600			
		} 
		"RT" {
			
			set wavelength(ch0,ex) 561
			set wavelength(ch1,ex) 640
			set wavelength(ch0,em) 600
			set wavelength(ch1,em) 670

		}
		}
	}
	"3" {
	set wavelength(ch0,ex) 405
	set wavelength(ch1,ex) 488
	set wavelength(ch2,ex) 561
	set wavelength(ch0,em) 440
	set wavelength(ch1,em) 520
	set wavelength(ch2,em) 600
	set globalSegmPars false
	}
}	


# clear parameters if they already exist
if {[info exists garbage_value] == 1} {
	unset garbage_value
	}

if {[info exists channel] == 1} {
		unset channel
	}
if {[info exists threshold_value] == 1} {
	unset threshold_value
	}
	
	# thresholding_method fixed search or otsu
	
# sets segmentation parameters
switch $thresholding_method {
	"otsu" {
			# run using single directory. Files with 2 channels
			set garbage_value(ch0) 1000
			set garbage_value(ch1) 100
			set garbage_value(ch2) 100
	}
	"search" {
			# it is using channel 0 for the time being...
			set threshold_value(ch0) "1000,20000,10"
			set garbage_value(ch0) 25000

			set threshold_value(ch1) "1000,20000,10"
			set garbage_value(ch1) 100

			set threshold_value(ch2) "1000,20000,10"
			set garbage_value(ch2) 100
			
			# minimum difference in the number of objects for automatic thresholding
			set minDifference 3
		
	}
	"fixed" {

			if {$singleDir == "true"} { 
				# run using single directory. Files with 2 channels
				set threshold_value(ch0) 2000
				set garbage_value(ch0) 10000
				set threshold_value(ch1) 2000
				set garbage_value(ch1) 100
				set threshold_value(ch2) 1000
				set garbage_value(ch2) 100		
		
			} else {	
				# run using specific segmentation parameters
				# DAPI
				set threshold_value(RTPI) 2500
				set garbage_value(RTPI) 10000
				# beads
				set threshold_value(beads) 2500
				set garbage_value(beads) 10000		
				# RTs
				set threshold_value(RT1) 2500
				set garbage_value(RT1) 10000
				set threshold_value(RT2) 2500
				set garbage_value(RT2) 10000
				set threshold_value(RT3) 2500
				set garbage_value(RT3) 10000
				set threshold_value(RT4) 2500
				set garbage_value(RT4) 10000
				set threshold_value(RT5) 2500
				set garbage_value(RT5) 10000
				set threshold_value(RT6) 2500
				set garbage_value(RT6) 10000
				set threshold_value(RT7) 2500
				set garbage_value(RT7) 10000
				set threshold_value(RT8) 2500
				set garbage_value(RT8) 10000
				set threshold_value(RT9) 2500
				set garbage_value(RT9) 10000
				set threshold_value(RT10) 2500
				set garbage_value(RT10) 10000
				set threshold_value(RT11) 2500
				set garbage_value(RT11) 10000
				set threshold_value(RT12) 2500
				set garbage_value(RT12) 10000
				set threshold_value(RT13) 2500
				set garbage_value(RT13) 10000
				set threshold_value(RT14) 2500
				set garbage_value(RT14) 10000
				set threshold_value(RT15) 2500
				set garbage_value(RT15) 10000
				set threshold_value(RT16) 2500
				set garbage_value(RT16) 10000
				set threshold_value(RT17) 2500
				set garbage_value(RT17) 10000
				set threshold_value(RT18) 2500
				set garbage_value(RT18) 10000
				set threshold_value(RT19) 2500
				set garbage_value(RT19) 10000
				}		
	}
}
