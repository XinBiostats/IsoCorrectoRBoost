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
**Note**:  
- **a** The molecule name must be consistent with that in dataset.  
- **b**	If the molecule only contains one specific element atom, you need to add number “1” after that element. For example, L.Glutamic.acid, C5H9N"1"O4LabC5.  
- **c**	For each molecule(-fragment) to be corrected, the number of atoms of all elements relevant for correction needs to be given in a sum formula, for example: C6H12O2N1LabC2. The prefix Lab marks the tracer element. In the example, C6 indicates that there are in total 6 atoms of carbon in the molecule or fragment considered. Then, LabC2 provides the information that of those 6 carbons, 2 positions may actually be labeled due to incorporation from the tracer substrate. The other 4 positions cannot contain tracer from the tracer substrate e.g. because they stem from derivatization.  

If any question, please check 6.1.1 in IsoCorrectoR vignettes [here](https://www.bioconductor.org/packages/devel/bioc/vignettes/IsoCorrectoRGUI/inst/doc/IsoCorrectoRGUI.html#input-files-and-parameters).


