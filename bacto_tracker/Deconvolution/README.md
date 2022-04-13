# huygensDeconvolutionScripts
Set of scripts to deconvolve images in huygens in batch mode





## Installing

To install, run this command in the folder where you want to install  ```huygensDeconvolutionScripts```

```bash
$ git clone https://github.com/marcnol/huygensDeconvolutionScripts
```

or your preferred method to clone a github repository.

If you don't have access, let #marcnol know!



You should now have a folder called ```huygensDeconvolutionScripts```



## Running 

Run the deconvolution script is pretty simple.



1. First login into ``lopevi`` or in 192.168.8.7



2. Modify ```batchDeconvolve.tcl``` with your input (i.e. ```dataDir```) and output (i.e. ```destDir```)  folders as follows:

```tcl
# batchDeconvolve.tcl
set tcl_files "/home/marcnol/Repositories/huygensDeconvolutionScripts/"

set dataDir "/mnt/grey/DATA/rawData_2020/Experiment_5/DAPI/005_merFISH_RAMM_DAPI"
set destDir "/home/marcnol/tmp/"

source "${tcl_files}merfish.tcl"

decon_segment "$tcl_files" "$dataDir" "$destDir"
```

```tcl_files``` only has to be modified the first time you run the script to indicate the folder where you store your scripts.



3. Change parameters in `merfish_parameters.tcl` if needed (i.e. `$experiment_type, $Nchannels`)



4. Then run the ```batchDeconvolve.tcl``` as follows

```bash
$ hucore -task /home/marcnol/Repositories/huygensDeconvolutionScripts/batchDeconvolve.tcl
```

adapt the directory above to the place where you stored your scripts.



## Disconnect a terminal without losing your run

This is easily done by resorting to ```screen```, a command line Linux program.



Start a screen session called ```hucore``` by running

```bash
screen -S hucore
```



Then run hucore

```bash
$ hucore -task /home/marcnol/Repositories/huygensDeconvolutionScripts/batchDeconvolve.tcl
```



Then click on ```ctrl+a```  then ```d``` . This command will "detach" your terminal.



You can now logout.



If you want to see where the script is in deconvolving images, you can log back into ```lopevi```, then run:

```bash
$ screen -x hucore
```

You should now see again hucore output.



If you want to see all the open screen sessions, run:

```bash
$ screen -ls
```

you should see something like:

```bash
There are screens on:
	16863.hucore	(07/27/2020 02:17:53 PM)	(Detached)
```



