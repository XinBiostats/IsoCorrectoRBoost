wd <- '/Users/xin.ma/Desktop/sunlab/NIA/dataset1/'
setwd(wd)

source('./code/source.R')
source('./code/util.R')


# Format Data
DataFormatted <- FormatData(RawMeasurementFile = './data_to_strip/tracer_1_Normal_subset_NIA.csv',
                            FirstCompound = 'L.Lactic.acid_.M.Hm0',
                            RawMoleculeFile = './data_source/MoleculeFile.csv')

# Delete low-quality pixels
DataQC <- rmNullPixel(MeasurementFileDir = './data_to_strip/',
                      pattern = "formatted\\.csv",
                      csvReturn = TRUE,
                      OnlyDeletePixelsWOIsotopologs = FALSE,
                      verbose = FALSE,
                      verboseFeature = FALSE,
                      rmDataStore = "NewDir",
                      outdir = "rmOutput")

# NIA Stripping
CorrectedData <- NIAStripping(MeasurementFile = './rmOutput/tracer_1_Normal_subset_NIA_formatted_rm0.csv',
                              ElementFile = './data_source/ElementFile.csv',
                              MoleculeFile = './data_source/MoleculeFile.csv',
                              ObservationInfo = DataFormatted$ObservationInfo,
                              CorrectTracerImpurity = FALSE,
                              CorrectTracerElementCore = TRUE,
                              CalculateMeanEnrichment = TRUE,
                              UltraHighRes = FALSE,
                              DirOut = "./OutputIsoCorrectoRBoost",
                              FileOut = "result",
                              FileOutFormat = "csv",
                              ReturnResultsObject = TRUE,
                              CorrectAlsoMonoisotopic = FALSE,
                              CalculationThreshold = 10^-8,
                              CalculationThreshold_UHR = 8,
                              parallel = TRUE,
                              verbose=TRUE,
                              Testmode = FALSE)


