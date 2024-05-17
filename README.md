# IsoCorrectoRBoost
This is a faster tool for Natural Isotopes Abundance (NIA) correction based on [**IsoCorrectoR**](https://genomics.ur.de/files/IsoCorrectoR/.). We have optimized some details based on IsoCorrector, resulting in a significant improvement in correction speed, particularly for large datasets.

## Data Preparation
We need 3 files for NIA stripping: Measurement File, Element File and Molecule File.
### Measurement File
This file contains the measured data that needs to be corrected. See example: [MeasurementFile](data_to_strip/tracer_1_Normal_subset_NIA.csv)
### Element File
All the elements that occur in the molecules to be corrected. See example: [ElementFile](data_source/ElementFile.csv)  
If any question, please check 6.1.3 in IsoCorrectoR vignettes [here](https://www.bioconductor.org/packages/devel/bioc/vignettes/IsoCorrectoRGUI/inst/doc/IsoCorrectoRGUI.html#input-files-and-parameters).
### Molecule File
All the molecules to be corrected. See example: [MoleculeFile](data_source/MoleculeFile.csv)
If any question, please check 6.1.1 in IsoCorrectoR vignettes [here](https://www.bioconductor.org/packages/devel/bioc/vignettes/IsoCorrectoRGUI/inst/doc/IsoCorrectoRGUI.html#input-files-and-parameters).
