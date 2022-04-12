# bacto_tracker

This repository contains the code implemented to segment and track bacterial cells in crowded environments.

Software packages were written by Sara Rombouts, Jean-Bernard Fiche and Marcelo Nollmann at the Center for Structural Biology, a joint Institute from the Centre National de la Recherche Scientifique (CNRS), the Institut National de la Sante et la Recherche Medicale (INSERM), and the University of Montpellier, France.

Code used for the acquisition and analysis of HiM data are included in the `/src` folder.



## Documentation

A full documentation can be found [here]().



## Data

Example datasets can be obtained from Zenodo [here](). 

After running Steps 1-4 in your raw datasets, you will obtain a directory structure as follows:

```
.
├── Tiling_Drift_PostProcess
├── Tracking

```

To convert the MATLAB output files into CSVs containing single cell tracks, you will need to run `Convert_output_to_csv.m` followed by `loadDataMATLAB.py`.

`Convert_output_to_csv.m` will create a folder `/Track_IDs` with CSV files where the tracks are encoded by two parameters: `frame number, cell ID number`. In addition, `Convert_output_to_csv.m` will generate a `/masks` folder that contains CSV files encoding `cell ID, x-centroid (px), y-centroid (px)` for each frame.

These two collections of CSV files will be used by `loadDataMATLAB.py` to generate a single CSV file where each single cell trajectory is encoded by `ID, frame number, x-centroid (px), y-centroid (px)`. 

After this conversion is made, you will be able to calculate track properties, as in the example in this [Jupyter lab]().



## Contact

Go to the [NollmannLab](https://www.nollmannlab.org/contact.html) group web-page.

