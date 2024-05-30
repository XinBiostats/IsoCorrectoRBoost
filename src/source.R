######### source ###########
library(dplyr)
library(data.table)
############################
### Data Preprocessing
############################

list2df <- function(x){
  MAX.LEN <- max(sapply(x, length), na.rm = TRUE)
  DF <- data.frame(lapply(x, function(x) c(x, rep(NA, MAX.LEN - length(x)))))
  colnames(DF) <- paste("V", seq(ncol(DF)), sep = "")
  return(DF)
}

order_molecule <- function(x){
  split_molecule <- strsplit(x,"_")
  molecule <- sapply(split_molecule, '[', 1)
  suffix <- as.numeric(sapply(split_molecule, '[', 2))
  sorted_indices <- order(molecule,suffix)
  return(sorted_indices)
}

FormatData <- function(RawMeasurementFile=NA,FirstCompound,RawMoleculeFile=NA){
  FormattedDataList <- list()
  
  # format MeasurementFile
  if (!is.na(RawMeasurementFile)){
    data <- fread(RawMeasurementFile,drop=1)
    
    selected_cols <- names(data)[which(names(data)==FirstCompound):ncol(data)]
    
    subdata <- data %>%
      select(all_of(selected_cols)) #drop redundant columns
    
    subdata_t <- t(subdata) #transpose
    
    rownames(subdata_t) <- gsub("^(.*?)_.*?(\\d+(\\.\\d+)?)$", "\\1_\\2", rownames(subdata_t)) # reformat molecule name
    rownames(subdata_t) <- gsub('\\.','-',rownames(subdata_t))
    
    unique_molecules <- unique(sapply(rownames(subdata_t), function(x) strsplit(x, "_")[[1]][1], USE.NAMES = FALSE)) # check unique molecules
    
    molecule_order <- order_molecule(rownames(subdata_t)) # reorder molecules
    subdata_t <- subdata_t[molecule_order,]
    colnames(subdata_t) <- data$spotId #set colnames
    
    # remove those molecules without any isotope
    molecules <- apply(X = list2df(strsplit(x = rownames(subdata_t),split = "_"))[1,], MARGIN = 2, FUN = as.character)
    molecules_count <- table(molecules)
    molecules_filtered <- names(molecules_count[molecules_count != 1])
    subdata_t = subdata_t[which(molecules %in% molecules_filtered),]
    
    newfile <- sub("\\.csv$", "_formatted.csv", RawMeasurementFile)
    
    write.csv(subdata_t, file=newfile,row.names=TRUE)
    
    # extract observations information
    ObsInfo <- data %>%
      select(colnames(data)[0:(which(colnames(data)==FirstCompound)-1)])
    
    FormattedDataList$ObservationInfo <- ObsInfo
    FormattedDataList$FormattedMeasurementData <- subdata_t
  }
  
  # format MoleculeFile
  if (!is.na(RawMoleculeFile)){
    MoleculeData <- fread(RawMoleculeFile)
    MoleculeData$Molecule <- gsub('\\.','-',MoleculeData$Molecule)
    write.csv(MoleculeData, file=RawMoleculeFile,row.names=FALSE)
    FormattedDataList$FormattedMoleculeData <- MoleculeData
  }
  
  return(FormattedDataList)
}


rmNullPixel <- function(MeasurementFileDir = NA,
                            pattern = "formatted\\.csv",
                            csvReturn = TRUE,
                            OnlyDeletePixelsWOIsotopologs = FALSE,
                            verbose = FALSE,
                            verboseFeature = FALSE,
                            rmDataStore = c("NewDir", "InputDir"),
                            outdir = "rmOutput"){
  
  
  ### functions needed
  
  #### list2df
  
  list2df <- function(x)
  {
    MAX.LEN <- max(sapply(x, length), na.rm = TRUE)
    DF <- data.frame(lapply(x, function(x) c(x, rep(NA, MAX.LEN - length(x)))))
    colnames(DF) <- paste("V", seq(ncol(DF)), sep = "")
    DF
  }
  
  #### Null Pixel Deplete
  
  NullPixel_deplet <- function(MeasurementFileDir,
                               verbose = F,
                               OnlyDeleteIsotopologPeak) {
    
    MeasurementFile <- read.csv(file = MeasurementFileDir, header = T, row.names = 1)
    
    lipids <- unique(apply(X = list2df(strsplit(x = rownames(MeasurementFile),
                                                split = "_"))[1,],
                           MARGIN = 2,
                           FUN = as.character))
    
    list_isotopologs <- list()
    
    count = 0
    
    for (i in 1:length(lipids)) {
      
      count = count + 1
      
      if (verbose == T) {
        
        cat("...", "\n")
        cat("Processing Feature No. :", i, "i.e.,", lipids[i])
        cat("...", "\n")
        
      }
      
      isotopolog_envelope <- MeasurementFile[grep(pattern = lipids[i], x = rownames(MeasurementFile)),]
      
      count_NA_0 = 0
      
      setNA <- function(col,OnlyDeleteIsotopologPeak){
        
        sum_isotopolog <- sum(col[-1]) 
        
        if (sum_isotopolog == 0){
          return (NA)
        } else if (sum_isotopolog != 0 & col[1] == 0 & OnlyDeleteIsotopologPeak == F){
          return (NA)
        } else {
          return (col)
        }
      }
      
      isotopolog_envelope[] <- lapply(isotopolog_envelope,setNA,OnlyDeletePixelsWOIsotopologs)
      list_isotopologs[[i]] <- isotopolog_envelope
    }
    
    
    ### assembling a matrix from the different length isotopologs
    
    out_df <- matrix(NA, nrow = 0, ncol = dim(MeasurementFile)[2])
    
    for (i in 1:length(list_isotopologs)) {
      
      runner <- list_isotopologs[[i]]
      
      out_df <-  rbind(out_df, runner)
      
    }
    return(out_df)
  }
  
  
  ### Main
  
  reps_test <- list.files(path = MeasurementFileDir, pattern = pattern,
                          all.files = FALSE, full.names = TRUE, recursive = TRUE,
                          ignore.case = FALSE, include.dirs = TRUE, no.. = FALSE)
  
  if (length(reps_test) > 0) {
    
    
    reps <- list.files(path = MeasurementFileDir, pattern = pattern,
                       all.files = FALSE, full.names = TRUE, recursive = TRUE,
                       ignore.case = FALSE, include.dirs = TRUE, no.. = FALSE)
    
    
    count = 0
    
    zero_rm_MatList <- list()
    
    for (i in 1:length(reps)) {
      
      count = count + 1
      
      if (verbose == T) {
        
        cat("...", "\n")
        cat("Processing File No. :", i, "i.e.,", reps[i])
        cat("...", "\n")
        
      }
      
      zero_rm_MatList[[i]] <- NullPixel_deplet(MeasurementFileDir = reps[i],
                                               verbose = verboseFeature,
                                               OnlyDeleteIsotopologPeak = OnlyDeletePixelsWOIsotopologs)
      
    }
    
  } else {
    
    stop("ERROR: no files found under the given pattern parameter in the function call")
    
  }
  
  
  if (csvReturn == T) {
    
    for (i in 1:length(zero_rm_MatList)) {
      
      test_print <- cbind(rownames(zero_rm_MatList[[i]]), zero_rm_MatList[[i]])
      colnames(test_print) <- c("Measurements/Samples", colnames(zero_rm_MatList[[i]]))
      
      file_name_spl <- strsplit(reps[i], split = "/")[[1]]
      
      rep_name <- strsplit(file_name_spl[[length(strsplit(reps[i], split = "/")[[1]])]], "\\.")[[1]][1]
      
      if (rmDataStore == "NewDir") {
        
        dir.create(path = outdir, showWarnings = F)
        
        file_name <- paste0(outdir, "/", rep_name,  "_rm0.csv")
        
        
      } else if (rmDataStore == "InputDir"){
        
        file_name <- paste0(paste0(file_name_spl[-length(file_name_spl)], collapse = "/"), "/", rep_name,  "_rm0.csv")
        
      }
      
      write.csv(x = test_print, file = file_name, row.names = F, quote = F)
      
    }
    
  } else {
    
  }
  
  return(zero_rm_MatList)
  
}
############################
### ISOCorrectoR
############################
# Element Parameter Extraction

#' @importFrom stringr str_split
#' @importFrom tibble as_tibble
#' @importFrom magrittr '%>%'
#' 
ElementInfoExtraction <- function(ElementData, UltraHighRes, logEnvironment, verbose) {
  if(verbose){message(date(), " :: processing element file ...")}
  
  # Using two nested loops, the isotope information of each element in column 2 is extracted by splitting it using the delimiters '/'
  # (isotopes1/isotope2...) and '_' (abundance-mass shift of given isotope).
  
  ElementList <- list()
  
  for (Element in seq_len(nrow(ElementData))) {
    
    Isotopes <- stringr::str_split(ElementData[Element, 2], pattern = "/", simplify = TRUE)
    
    ElementArray <- matrix(NA, nrow = length(Isotopes), ncol = 2)
    colnames(ElementArray) <- c("IsotopeAbundance", "MassShift")
    for (Isotope in seq_len(length(Isotopes))) {
      
      IsoAbundanceMassShift <- stringr::str_split(Isotopes[Isotope], pattern = "_", simplify = TRUE)
      
      for (i in seq_len(ncol(ElementArray))) {
        ElementArray[Isotope, i] <- as.numeric(IsoAbundanceMassShift[i])
      }
      
    }  #Isotope
    
    # Sort Isotopes of each element in ElementArray ascending by probability. This is required to be able to later calculate probabilities in descending
    # order, making it possible to stop the calculation at a defined threshold without losing higher probability values.
    
    ElementArray_df <- as.data.frame(ElementArray)
    
    ElementArray_df <- ElementArray_df[order(ElementArray_df$IsotopeAbundance, decreasing = TRUE),] 
    
    ElementArray_tbl <- tibble::as_tibble(ElementArray_df)
    
    tmp.list <- list()
    tmp2.list <- list()
    tmp.list[[1]] <- ElementArray_tbl
    tmp2.list[[1]] <- tmp.list
    
    # now add rest of ElementData
    for (i in 3:ncol(ElementData)) {
      tmp <- ElementData[Element, i]
      
      # Making NAs of ElementArray 0
      
      if (is.na(tmp)) {
        tmp <- 0
      }
      tmp2.list[[i - 1]] <- tmp
    }  #i
    
    names(tmp2.list) <- c("Isotopes", colnames(ElementData)[3:ncol(ElementData)])
    ElementList[[Element]] <- tmp2.list
    rm(tmp.list)
    rm(tmp2.list)
    
  }  #Element
  
  names(ElementList) <- ElementData[[1]]
  
  checkElementDataLogic(ElementList, UltraHighRes, logEnvironment, verbose=verbose)
  
  if(verbose){message(date(), " :: processing element file [OK]\n")}
  return(ElementList)
}  # end of function



# Load input files, check their structure, extract information and check input logic

ExtractInputFileInformation <- function(ElementFile, MoleculeFile,
                                        MeasurementFile, UltraHighRes,
                                        CorrectTracerImpurity, logEnvironment,
                                        verbose) {
  ElementData <-
    checkElementDataStructure(
      data = readFileByExt(ElementFile, verbose = verbose), logEnvironment =
        logEnvironment, verbose = verbose
    )
  MoleculeData <-
    checkMoleculeDataStructure(
      data = readFileByExt(MoleculeFile, verbose = verbose),
      logEnvironment = logEnvironment, verbose = verbose
    )
  RawData <-
    checkRawData(
      inputFile = MeasurementFile, logEnvironment = logEnvironment,
      verbose = verbose
    )
  
  # Extract user input from the element information file
  
  ElementInfo <-
    ElementInfoExtraction(
      ElementData = ElementData,
      UltraHighRes = UltraHighRes,
      logEnvironment = logEnvironment,
      verbose = verbose
    )
  
  # Extract user input from the molecule information file
  
  MoleculeInfo <- MoleculeInfoExtraction(
    MoleculeData = MoleculeData,
    ElementInfo = ElementInfo,
    UltraHighRes = UltraHighRes,
    CorrectTracerImpurity = CorrectTracerImpurity,
    MoleculesFound = RawData$MoleculesFound,
    logEnvironment = logEnvironment,
    verbose = verbose
  )
  
  # Automatic determination of the transitions/measured isotopologues expected
  # according to the molecule parameters entered in the molecule information file
  
  if (!UltraHighRes) {
    MoleculeInfo <- CalculateTransitionsNR(
      MoleculeInfo = MoleculeInfo,
      ElementInfo = ElementInfo,
      verbose = verbose
    )
  } else {
    MoleculeInfo <- CalculateTransitionsUHR(
      MoleculeInfo = MoleculeInfo,
      ElementInfo = ElementInfo,
      verbose = verbose
    )
  }
  
  # Extraction of measurement information from the measurement file provided
  
  resultRawDataExtractionList <-
    RawDataExtraction(
      data = RawData[[1]],
      MoleculeArray = MoleculeInfo,
      logEnvironment = logEnvironment,
      verbose = verbose
    )
  
  MoleculeInfo <- resultRawDataExtractionList[["MoleculeInfo"]]
  
  MeasurementInfo <- resultRawDataExtractionList[["dataRaw"]]
  
  return(list(
    "MoleculeInfo" = MoleculeInfo,
    "ElementInfo" = ElementInfo,
    "MeasurementInfo" = MeasurementInfo
  ))
}

# Extraction of molecule parameters

#' @importFrom stringr str_c
#' @importFrom stringr str_detect
#' @importFrom stringr str_replace
#' @importFrom stringr str_extract
#' @importFrom stringr str_extract_all
#' @importFrom magrittr '%>%'
#'  
MoleculeInfoExtraction <- function(MoleculeData, ElementInfo, UltraHighRes, CorrectTracerImpurity, MoleculesFound, logEnvironment, verbose) {
  
  if(verbose){message(date(), " :: processing molecule file ...")}
  
  #Find molecules from measurement file in molecule data
  
  rownames(MoleculeData) <- as.character(MoleculeData[, 1])
  
  MoleculeLocation.vec <- vector()
  MoleculeLocationLabel.vec <- vector()
  
  for (MoleculeFound in MoleculesFound) {
    
    loc <- grep(paste0("^", MoleculeFound, "$"), rownames(MoleculeData))
    
    if (length(loc) == 1) {
      MoleculeLocation.vec <- c(MoleculeLocation.vec, loc)  # simply add current MoleculeName to vector of all valid MoleculeNames
      MoleculeLocationLabel.vec <- c(MoleculeLocationLabel.vec, MoleculeFound)
    } else {
      notification <- stringr::str_c("Molecule '", MoleculeFound, "' from the measurement data was not found in the molecule file.")
      errorHandler(notification, logEnvironment, "error", verbose=verbose)
    }
  }  #MoleculeFound
  
  names(MoleculeLocation.vec) <- MoleculeLocationLabel.vec
  
  MoleculeList <- list()  # big list containing all necessary information on a molecule
  
  # MoleculeLocation.vec contains location indices of found molecules !
  
  #Analyze imported molecule formulas with regex
  
  for (MoleculeNo in seq_len(length(MoleculeLocation.vec))) {
    
    MoleculeLocation <- MoleculeLocation.vec[MoleculeNo]
    
    MoleculeInformation <- MoleculeData[MoleculeLocation,]
    
    if(verbose){message(date(), " :: :: found molecule: #", MoleculeNo, " [", names(MoleculeLocation.vec)[MoleculeNo], "]")}
    # Check if MoleculeFile contains entries in both the product ion/MS ion 
    # and the neutral loss column or only in the product ion/MS ion column.  
    # This yields NumberFragments for a given molecule and determines whether 
    # MS or MS/MS correction will be applied in the following.
    
    FragmentList <- ParseMoleculeInformation(MoleculeInformation=MoleculeInformation,
                                             ElementInfo=ElementInfo,UltraHighRes=UltraHighRes,
                                             verbose=verbose)
    
    MoleculeList[[MoleculeNo]] <- FragmentList
    
  }  # MoleculeNo
  
  names(MoleculeList) <- names(MoleculeLocation.vec)
  
  #Check logic of extracted molecule information
  
  checkMoleculeDataLogic(MoleculeList, UltraHighRes, CorrectTracerImpurity, ElementInfo, logEnvironment, verbose=verbose)
  
  #Provide log information on molecules and tracers
  
  logEnvironment$param$molecules <- names(MoleculeList)
  
  tracerLogInfo <- list()
  tracerNames <- character()
  
  #Get tracerNames, a vector of all tracer elements used
  
  for (molecule in names(MoleculeList)) {
    
    for (fragmentNo in seq_len(sum(stringr::str_detect(names(MoleculeList[[molecule]]), "Fragment")))) {
      
      fragment <- paste0("Fragment_", fragmentNo)
      
      tracers <- MoleculeList[[molecule]][[fragment]]$Tracer
      
      if (length(tracers) > 0) {
        
        tracerNames <- c(tracerNames, names(tracers))
        
      }
      
    }
    
  }
  
  tracerNames <- unique(tracerNames)
  
  for (tracer in tracerNames) {
    
    if (CorrectTracerImpurity == TRUE) {
      
      tracerLogInfo[[tracer]] <- ElementInfo[[tracer]]$Tracer.purity
      
    } else {
      
      tracerLogInfo[[tracer]] <- NA
      
    }
    
  }
  
  logEnvironment$param$tracers <- tracerLogInfo
  
  if(verbose){message(date(), " :: processing molecule file [OK]\n")}
  
  return(MoleculeInfo = MoleculeList)
}  # function()


# Analysis of molecule formula with regex

# The strings imported from MoleculeFile containing information on the number 
# of atoms per element are analyzed using regular expressions.  The regular
# expression Characters extracts all characters, the regular expression 
# Numbers extracts all numbers.  TracerElement looks for the tracer element
# identifier.

ParseMoleculeInformation <- function(MoleculeInformation, ElementInfo, UltraHighRes, verbose) {
  
  Characters <- "([A-Za-z]*)"
  Numbers <- "([0-9]+[0-9]|[0-9])"
  TracerElement <- "(Lab)"
  MoleculeFormulaRegex <- "(|Lab)[A-Z][a-z]?[1-9][0-9]*"
  
  if (is.na(MoleculeInformation[1,3])) {
    NumberFragments <- 1
  } else {
    NumberFragments <- 2
  }
  
  # collect features for each individual fragment
  tmpFragmentList <- list()
  
  # FragmentList contains one or two tmpFragmentList()s
  FragmentList <- list()
  
  # The following is done for each Molecule and each fragment 
  # of a given molecule
  for (Fragment in seq_len(NumberFragments)) {
    
    # For each molecule(-fragment), the string containing the element count 
    # information is split into the individual elements using a regex
    
    Elements <- stringr::str_extract_all(MoleculeInformation[1,Fragment + 1], MoleculeFormulaRegex)[[1]]
    
    # Now each element in Elements, consisting of the element ID and the 
    # element count, is again separated into ID and element count using 
    # the regular expressions constructed above. 
    # Additionally, a check is made to identify the 
    # tracer element information ('Lab').
    
    ElementCount.vec <- vector()
    Element.vec <- vector()
    TracerCount.vec <- vector()
    Tracer.vec <- vector()
    
    NumberElementsandTracer <- length(Elements)
    
    for (Element in seq_len(NumberElementsandTracer)) {
      
      # TRUE => current Element is Tracer // FALSE => current Element is no Tracer
      if (Elements[Element] %>% stringr::str_detect(TracerElement)) {
        
        TracerCount.vec <- c(TracerCount.vec, Elements[Element] %>% stringr::str_replace(TracerElement, "") %>% stringr::str_extract(Numbers) %>% 
                               as.numeric)
        Tracer.vec <- c(Tracer.vec, Elements[Element] %>% stringr::str_replace(TracerElement, "") %>% stringr::str_extract(Characters))
      } else {
        
        ElementCount.vec <- c(ElementCount.vec, Elements[Element] %>% stringr::str_extract(Numbers) %>% as.numeric)
        Element.vec <- c(Element.vec, Elements[Element] %>% stringr::str_extract(Characters))
      }  #else
    }  #Element    
    names(ElementCount.vec) <- Element.vec
    names(TracerCount.vec) <- Tracer.vec
    
    tmpFragmentList <- list(Element = ElementCount.vec, Tracer = TracerCount.vec)
    
    NumberElements <- length(ElementCount.vec)
    NumberTracers <- length(TracerCount.vec)
    
    # 'NumberElements' - 'NumberTracers' gives the number of non-tracer elements NumberElementsNonTracer.
    NumberElementsNonTracer <- NumberElements - NumberTracers
    
    if (!UltraHighRes) {
      
      # If a tracer is present in the molecule(-fragment) considered, 
      # the tracer parameters are extracted.  'MaxLabel' is the maximum amount
      # of tracer isotope that is expected to be found in the 
      # molecule(-fragment) due to metabolism while 'nTracerMax' is 
      # the maximum amount of tracer element (labelled or
      # unlabelled) in that same species. 'IDTracer' is the tracer elements ID.
      
      if (NumberTracers > 0) {
        
        MaxLabel <- max(TracerCount.vec)
        IDTracer <- TracerCount.vec %>% which.max %>% names
        
        # This part deals with the possibility of having a tracer element 
        # that shows isotopes with a negative mass shift.  
        # In this case IsoCombinationsMaster()
        # has to calculate `PlacesToAssign` differently for the tracer.
        
        if (sum(ElementInfo[[IDTracer]][["Isotopes"]][[1]][["MassShift"]] < 0) > 0) {
          NegIsoTracer <- 1  # => TRUE
        } else {
          NegIsoTracer <- 2  # => FALSE
        }
        
        names(NegIsoTracer) <- IDTracer
        
        tmp.which <- which(names(ElementCount.vec) %in% names(TracerCount.vec))
        nTracerMax <- as.numeric(ElementCount.vec[tmp.which])
        
        if (length(nTracerMax) == NumberTracers) {
          
          names(nTracerMax) <- names(TracerCount.vec)
          
        }
        
        names(MaxLabel) <- IDTracer
        
        tmpFragmentList[["MaxLabel"]] <- MaxLabel
        tmpFragmentList[["IDTracer"]] <- IDTracer
        tmpFragmentList[["nTracerMax"]] <- nTracerMax
        tmpFragmentList[["NegIsoTracer"]] <- NegIsoTracer
      } else {
        tmpFragmentList[["MaxLabel"]] <- NA
        tmpFragmentList[["IDTracer"]] <- NA
        tmpFragmentList[["nTracerMax"]] <- NA
        tmpFragmentList[["NegIsoTracer"]] <- NA
      }  #if(NumberTracers>0)
      
      # In this section, the non-tracer element parameters are extracted in 
      # the same way as for the tracer element.
      
      ElementsNonTracerList <- list()
      
      if (NumberElementsNonTracer > 0) {
        
        NonTracer <- 1
        
        tmpElementNonTracer.vec <- vector()
        tmpElementNonTracerCount.vec <- vector()
        tmpElementZeroTracer.vec <- vector()
        tmpElementZeroTracerCount.vec <- vector()
        
        for (Element in seq_len(NumberElements)) {
          # tracers exist
          if (NumberTracers > 0) {
            # current element is no tracer
            if (names(ElementCount.vec)[Element] != IDTracer) {
              tmpElementNonTracer.vec <- c(tmpElementNonTracer.vec, names(ElementCount.vec)[Element])
              tmpElementNonTracerCount.vec <- c(tmpElementNonTracerCount.vec, as.numeric(ElementCount.vec[Element]))
            }
          } else {
            tmpElementZeroTracer.vec <- c(tmpElementZeroTracer.vec, names(ElementCount.vec)[Element])
            tmpElementZeroTracerCount.vec <- c(tmpElementZeroTracerCount.vec, as.numeric(ElementCount.vec[[Element]]))
            
          }  #NumberTracers
        }  #Element
        
        names(tmpElementNonTracerCount.vec) <- tmpElementNonTracer.vec
        names(tmpElementZeroTracerCount.vec) <- tmpElementZeroTracer.vec
        
        tmpFragmentList[["NonTracer"]] <- tmpElementNonTracerCount.vec
        tmpFragmentList[["ZeroTracer"]] <- tmpElementZeroTracerCount.vec
        
      } else {
        # NumberElementsNonTracer>0
        
        tmpFragmentList[["NonTracer"]] <- c()
        tmpFragmentList[["ZeroTracer"]] <- c()
        
      }
      
    } else if (UltraHighRes) 
    {
      
      if (NumberTracers > 0) {
        
        # For multiple tracer correction, the gathering of tracer parameters 
        # has to loop through the number of tracers present in a 
        # given molecule (-fragment).
        
        MaxLabel.vec <- vector()
        IDTracer.vec <- vector()
        nTracerMax.vec <- vector()
        
        for (TracerNo in seq_len(NumberTracers)) {
          MaxLabel <- as.numeric(TracerCount.vec[TracerNo])
          IDTracer <- names(TracerCount.vec)[TracerNo]
          
          MaxLabel.vec <- c(MaxLabel.vec, MaxLabel)
          IDTracer.vec <- c(IDTracer.vec, IDTracer)
          
          for (Element in seq_len(NumberElements)) {
            # if current element is a tracer
            if (names(tmpFragmentList[["Element"]][Element]) == IDTracer) 
            {
              nTracerMax <- tmpFragmentList[["Element"]][Element] %>% as.numeric
              nTracerMax.vec <- c(nTracerMax.vec, nTracerMax)
            }  #if
          }  #Element
        }  #TracerNo
        names(MaxLabel.vec) <- IDTracer.vec
        
        if (length(nTracerMax.vec) == NumberTracers) {
          
          names(nTracerMax.vec) <- IDTracer.vec
          
        }
        
        tmpFragmentList[["MaxLabel"]] <- MaxLabel.vec
        tmpFragmentList[["IDTracer"]] <- IDTracer.vec
        tmpFragmentList[["nTracerMax"]] <- nTracerMax.vec
        
        # Gaining natural abundance information associated with the 
        # tracer isotopes for probability calculations
        
        # Number of isotopes per Element
        NumberIso <- unlist(lapply(ElementInfo, function(x) nrow(data.frame(x[["Isotopes"]]))))
        
        ### store information for each individual Tracer Element
        NatAbuTracerList <- list()
        NatAbuBaseList <- list()
        MassShiftTracerList <- list()
        
        MassShiftTracer.vec <- vector()
        
        for (TracerNo in names(TracerCount.vec)) {
          MassShiftTracer <- ElementInfo[[TracerNo]][[2]]
          
          # for each Tracer isotope
          NatAbuTracer.vec <- vector()
          NatAbuBase.vec <- vector()
          
          for (IsotopeNo in seq_len(NumberIso[TracerNo])) {
            
            if (data.frame(ElementInfo[[TracerNo]][[1]])[IsotopeNo, 2] == MassShiftTracer) {
              
              NatAbuTracer.vec <- c(NatAbuTracer.vec, data.frame(ElementInfo[[TracerNo]][[1]])[IsotopeNo, 1])  # 1=>IsotopeAbundance
              
            } else if (data.frame(ElementInfo[[TracerNo]][[1]])[IsotopeNo, 2] == 0) {
              
              NatAbuBase.vec <- c(NatAbuBase.vec, data.frame(ElementInfo[[TracerNo]][[1]])[IsotopeNo, 1])
              
            } else {
              
              if(verbose){message(date(), " :: skipping IsotopeNo ", IsotopeNo, " for Tracer ", TracerNo)}
              
            }
          }  #IsotopeNo
          
          NatAbuTracerList[[TracerNo]] <- NatAbuTracer.vec
          NatAbuBaseList[[TracerNo]] <- NatAbuBase.vec
          MassShiftTracerList[[TracerNo]] <- MassShiftTracer
          
        }  #TracerNo
        
        tmpFragmentList[["NatAbuTracer"]] <- unlist(NatAbuTracerList)
        tmpFragmentList[["NatAbuBase"]] <- unlist(NatAbuBaseList)
        tmpFragmentList[["MassShiftTracer"]] <- unlist(MassShiftTracerList)
      } else {
        if(verbose){message(date(), " :: NO TRACER FOR MOLECULE ", MoleculeInformation[1,1], " AND FRAGMENT #", Fragment)}
      }  #NumberTracers>0
    }  #if(UltraHighRes==1)
    
    FragmentList[[Fragment]] <- tmpFragmentList
  }  # Fragment
  
  names(FragmentList) <- stringr::str_c("Fragment_", rep(seq_len(NumberFragments)))
  
  return(FragmentList)
  
}


CalculateTransitionsNR <- function(MoleculeInfo, ElementInfo, verbose) {
  
  if(verbose){message(date(), " :: calculating transitions ...")}
  MoleculesTotal <- length(MoleculeInfo)
  MoleculesName <- names(MoleculeInfo)
  
  # For all molecules, determine the maximum amount of tracer isotopes expected
  # in product ion or neutral loss.
  for (MoleculeNo in seq_len(MoleculesTotal)) {
    
    ### assemble required data
    MoleculeData <- MoleculeInfo[[MoleculeNo]]  # entire 'MoleculeInfo' data for current molecule
    
    NumberFragments <- sum(stringr::str_detect(names(MoleculeData), "Fragment"))  # %>%sum
    NumberTracers <- unlist(lapply(MoleculeData, function(x) length(x[["Tracer"]])))
    MaxLabel <- unlist(lapply(MoleculeData, function(x) x[["MaxLabel"]]))
    IDTracer <- unlist(lapply(MoleculeData, function(x) x[["IDTracer"]]))
    
    if (NumberFragments == 2) {
      
      #if(verbose){message(date(), " :: [", MoleculesName[MoleculeNo], "] // NumberFragments==", NumberFragments)}
      if(verbose){message(date(), " :: :: molecule [", MoleculesName[MoleculeNo], "] with ",NumberFragments," fragments")}
      if (NumberTracers[[1]] > 0) {
        MaxLabelProduct <- MaxLabel[[1]]
      } else {
        MaxLabelProduct <- 0
      }
      
      if (NumberTracers[[2]] > 0) {
        MaxLabelNeutralLoss <- MaxLabel[[2]]
      } else {
        MaxLabelNeutralLoss <- 0
      }
      if (NumberTracers[[1]] > 0) {
        Tracer <- IDTracer[[1]]
      } else if (NumberTracers[[2]] > 0) {
        Tracer <- IDTracer[[2]]
      } else {
        stop(date(), " :: [CalculateTransitions()] No tracer element specified for molecule ", MoleculesName[MoleculeNo], "!")
      }
      
      MaxLabelPrecursor <- MaxLabelProduct + MaxLabelNeutralLoss
      
      # CALCULATION OF MS/MS TRANSITIONS
      
      TransitionNo <- 0
      LabelNeutralLoss <- 0
      ProductIon_SiteSaturation <- 0
      SiteSaturationCorrection <- 0
      
      # This part of the function computes how label in the precursor ion can be distributed among product ion and neutral loss, within the constraints given by
      # MaxLabelProduct and MaxLabelNeutralLoss. This yields the expected MS/MS transitions.
      
      tmpTransitionsExpectedList <- list()
      tmpTransitionsExpectedListNames.vec <- vector()
      
      for (LabelPrecursor in 0:MaxLabelPrecursor) {
        # Until the label in the precursor has exceeded the amount of MaxLabelNeutralLoss, all label is considered to be in the neutral loss at this stage.
        
        if (LabelNeutralLoss < MaxLabelNeutralLoss) {
          LabelNeutralLoss <- LabelPrecursor
        }
        
        # ProductIon_SiteSaturation is a means to determine by how much the label in the precursor exceeds the maximum capacity of the product ion given in
        # MaxLabelProduct.
        
        ProductIon_SiteSaturation <- (MaxLabelProduct - LabelPrecursor) * (-1)
        if (ProductIon_SiteSaturation > 0) {
          SiteSaturationCorrection <- ProductIon_SiteSaturation
        }
        
        # At this stage, the MS/MS transitions are calculated and written into the molecule information list.  The vector named Precursor contains the mass shift
        # (in relation to M+0) associated with the precursor (Number of label in the precursor multiplied with the tracer isotope mass shift given in
        # ElementList).  The vector named NeutralLoss contains the mass shift of the neutral loss. This is determined in a loop which starts at the maximum
        # possible amount of current LabelPrecursor from the parent loop in the neutral loss. It then decreases the amount of LabelNeutralLoss by 'shifting' label
        # from the neutral loss to the product ion (NLLabel_in_ProductIon). Consequently, the product ion mass shift in the vector named ProductIon is given as
        # the difference between precursor mass shift and neutral loss mass shift of the current loop iteration.  This way all combinations of product ion and
        # neutral loss labelling for a given precursor labelling state are derived. The shifting just described has limits given by MaxLabelProduct.  This is
        # where SiteSaturationCorrection becomes relevant.  If the amount of label in the precursor exceeds the maximum amount of labelling possible in the
        # product ion, it will reduce the number of iterations of the shifting loop by just this difference. In addition, the loop generates specific name tags
        # for each transition of a given molecule by combining the molecule name from the molecules information file with '_T' for transition, followed by the
        # precuror label x and product ion label y as x.y (Name_x.y).
        
        # TransitionNo<-0
        for (NLLabel_in_ProductIon in 0:(LabelNeutralLoss - SiteSaturationCorrection)) {
          
          TransitionNo <- TransitionNo + 1
          
          tmpTrans1 <- LabelPrecursor * ElementInfo[[Tracer]][[2]]
          tmpTrans2 <- (LabelNeutralLoss - NLLabel_in_ProductIon) * ElementInfo[[Tracer]][[2]]
          tmpTrans3 <- tmpTrans1 - tmpTrans2
          
          tmpTransitionsExpectedList[[TransitionNo]] <- list(ProductIon = tmpTrans3, NeutralLoss = tmpTrans2, Precursor = tmpTrans1)
          
          tmpTransitionsExpectedListNames.vec <- c(tmpTransitionsExpectedListNames.vec, stringr::str_c(MoleculesName[MoleculeNo], "_", tmpTrans1, ".", tmpTrans3))
          
        }  # NLLabel
        
      }  # LabelPrecursor
      
      names(tmpTransitionsExpectedList) <- tmpTransitionsExpectedListNames.vec
      
      TransitionsExpected.df <- data.frame(ProductIon = as.numeric(unlist(lapply(tmpTransitionsExpectedList, function(x) x["ProductIon"]))), NeutralLoss = as.numeric(unlist(lapply(tmpTransitionsExpectedList, 
                                                                                                                                                                                    function(x) x["NeutralLoss"]))), Precursor = as.numeric(unlist(lapply(tmpTransitionsExpectedList, function(x) x["Precursor"]))))
      
      rownames(TransitionsExpected.df) <- tmpTransitionsExpectedListNames.vec
      
      TransitionsExpected <- TransitionsExpected.df
      
      # In this section the expected MS (not MS/MS) measurements and their name tags are generated for each molecule An MS name tag has the structure 'Name_x'
      
    } else if (NumberFragments == 1) 
    {
      
      if (NumberTracers[[1]] > 0) {
        MaxLabelProduct <- MaxLabel[[1]]
        Tracer <- IDTracer[[1]]
      } else {
        MaxLabelProduct <- 0
        warning(date(), " [WARNING] No tracer element specified for Molecule #", MoleculeNo, ".")
      }
      
      TransitionNo <- 0
      tmpTransitionsExpectedList <- list()
      tmpTransitionsExpectedListNames.vec <- vector()
      for (LabelIon in 0:MaxLabelProduct) {
        TransitionNo <- TransitionNo + 1
        tmpTrans1 <- LabelIon * ElementInfo[[Tracer]][[2]]
        
        tmpTransitionsExpectedList[[TransitionNo]] <- list(ProductIon = tmpTrans1, Precursor = tmpTrans1)
        tmpTransitionsExpectedListNames.vec <- c(tmpTransitionsExpectedListNames.vec, stringr::str_c(MoleculesName[MoleculeNo], "_", tmpTrans1))
        
      }  #LabelIon
      names(tmpTransitionsExpectedList) <- tmpTransitionsExpectedListNames.vec
      
      TransitionsExpected.df <- data.frame(ProductIon = as.numeric(unlist(lapply(tmpTransitionsExpectedList, function(x) x["ProductIon"]))), Precursor = as.numeric(unlist(lapply(tmpTransitionsExpectedList, 
                                                                                                                                                                                  function(x) x["Precursor"]))))
      
      
      rownames(TransitionsExpected.df) <- tmpTransitionsExpectedListNames.vec
      
      TransitionsExpected <- TransitionsExpected.df
      
    }  # NumberFragments==1
    MoleculeInfo[[MoleculeNo]][["TransitionsExpected"]] <- TransitionsExpected
  }  #  MoleculeNo
  
  if(verbose){message(date(), " :: calculating transitions [OK]\n")}
  return(MoleculeInfo)
  
}


#' @importFrom magrittr '%>%'

RawDataExtraction <- function(data, MoleculeArray, logEnvironment, verbose) {
  
  if(verbose){message(date(), " :: processing raw data ...")}
  
  for (i in 2:ncol(data)) {
    data[, i] <- as.numeric(data[, i])
  }
  
  # remove rows that are ALL NA (=> don't remove individual NAs !)
  data <- data[!apply(apply(data, 1, is.na), 2, sum) == ncol(data), ]
  # remove columns that are all NA
  data <- data[, !apply(apply(data, 2, is.na), 2, sum) == nrow(data)]
  
  rownames(data) <- as.character(data[, "Measurements/Samples"])
  data <- as.matrix(data[, -which(colnames(data) %in% "Measurements/Samples"), drop = FALSE])
  txt_extracted <- rownames(data)
  
  # Extraction of measurement locations from the imported data to check for presence.
  
  MoleculesTotal <- length(MoleculeArray)
  MoleculeName <- names(MoleculeArray)
  
  MoleculeLocationList <- list()
  NumberTransitionsList <- list()
  TransitionLocationExpectedList <- list()
  TransitionLocationList <- list()
  MissingMoleculesList <- list()
  MissingTransitionsList <- list()
  
  for (MoleculeNo in seq_len(MoleculesTotal)) {
    
    MoleculeData <- MoleculeArray[[MoleculeNo]]  # entire 'MoleculeInfo' data for current molecule
    
    TransitionsExpected <- MoleculeData[["TransitionsExpected"]]  # from CalculateTransitions()
    NumberTransitions <- nrow(TransitionsExpected)
    
    # Check: Do the measurement labels in the first column of the data match the theoretically expected measurement labels?
    
    MissingTransitions.vec <- vector()
    TransitionLocationExpected.vec <- vector()
    
    for (TransitionNo in seq_len(NumberTransitions)) {
      
      tmpIdx <- stringr::str_detect(txt_extracted, paste0("^", rownames(TransitionsExpected)[TransitionNo], "$")) %>% which
      
      MissingTransitions.vec[TransitionNo] <- 0
      
      if (length(tmpIdx) == 0) {
        TransitionLocationExpected.vec[TransitionNo] <- 0
        MissingTransitions.vec[TransitionNo] <- 1
        
        notification <- stringr::str_c("In measurement data file: The expected measurement ID '", rownames(TransitionsExpected)[TransitionNo], "' could not be found in the measurement data file.\nCorrection will be performed, however the results may be less accurate. Please check your input file for typos.", 
                                       "\nBe especially careful when considering fraction and mean enrichment values from molecules with missing measurements.")
        errorHandler(notification, logEnvironment, "warning", verbose=verbose)
        
      } else {
        
        TransitionLocationExpected.vec[TransitionNo] <- tmpIdx
        
      }
      
    }  #TransitionNo
    
    if (sum(TransitionLocationExpected.vec) == 0) {
      
      notification <- stringr::str_c("In measurement data file: None of the expected measurement IDs found for molecule '", names(MoleculeArray[MoleculeNo]), 
                                     "'.")
      errorHandler(notification, logEnvironment, "error", verbose=verbose)
      
    }
    
    names(TransitionLocationExpected.vec) <- rownames(TransitionsExpected)
    names(MissingTransitions.vec) <- rownames(TransitionsExpected)
    
    # If measurement IDs are missing, a new list of Measurement IDs is produced containing only those that were actually found in the input data.  This list
    # is then used for further computations.
    
    Transitions <- TransitionsExpected[MissingTransitions.vec == 0, , drop = FALSE]
    TransitionLocation <- TransitionLocationExpected.vec[MissingTransitions.vec == 0]
    
    TransitionLocationExpectedList[[MoleculeNo]] <- TransitionLocationExpected.vec
    TransitionLocationList[[MoleculeNo]] <- TransitionLocation
    
    MoleculeArray[[MoleculeNo]][["Transitions"]] <- Transitions
    MoleculeArray[[MoleculeNo]][["MissingTransitions"]] <- MissingTransitions.vec
    MoleculeArray[[MoleculeNo]][["TransitionLocationExpected"]] <- TransitionLocationExpected.vec
    MoleculeArray[[MoleculeNo]][["TransitionLocation"]] <- TransitionLocation
    
  }  #MoleculeNo
  
  if(verbose){message(date(), " :: processing raw data [OK]\n")}
  return(list(MoleculeInfo = MoleculeArray, dataRaw = data))
  
}  #RawDataExtraction


# Computation of isotope combinations, associated probabilities and mass shifts 
# for an element

#' @importFrom stringr str_c

IsoCombinationsMaster <- function(MoleculeData, MoleculeName, MoleculeNo, Fragment, ElementInfo, CalculationThreshold, CorrectTracerImpurity, CorrectTracerElementCore) {
  
  MoleculeFragmentData <- MoleculeData[[stringr::str_c("Fragment_", Fragment)]]  # get individual fragment data
  
  Transitions <- MoleculeData[["Transitions"]]
  
  NumberTransitions <- nrow(Transitions)
  
  NumberTracers <- length(MoleculeFragmentData[["Tracer"]])
  
  # Isotope combinations for tracer element
  
  IsoCombinationsTracerResult <- IsoCombinationsTracer(MoleculeFragmentData=MoleculeFragmentData,
                                                       Fragment=Fragment,
                                                       ElementInfo = ElementInfo,
                                                       Transitions=Transitions,
                                                       NumberTracers=NumberTracers,
                                                       CorrectTracerImpurity=CorrectTracerImpurity, 
                                                       CorrectTracerElementCore=CorrectTracerElementCore,
                                                       CalculationThreshold=CalculationThreshold)
  
  ProbTracerList <- IsoCombinationsTracerResult$ProbTracerList
  MassShiftTracerList <- IsoCombinationsTracerResult$MassShiftTracerList
  NumberIsoCombTracerEff <- IsoCombinationsTracerResult$NumberIsoCombTracerEff
  IsotopesTracer <- IsoCombinationsTracerResult$IsotopesTracer
  LabelPresent <- IsoCombinationsTracerResult$LabelPresent
  
  # Isotope combinations for non-tracer elements
  
  if (NumberTracers > 0) {
    
    IDTracer <- MoleculeFragmentData[["IDTracer"]] 
    
  } else {
    
    IDTracer <- NA
    
  }
  
  IsoCombinationsNonTracerResult <- IsoCombinationsNonTracer(MoleculeFragmentData=MoleculeFragmentData, 
                                                             Fragment=Fragment,
                                                             ElementInfo=ElementInfo,
                                                             NumberTracers=NumberTracers,
                                                             IDTracer=IDTracer,
                                                             CalculationThreshold=CalculationThreshold)
  
  ProbElemList <- IsoCombinationsNonTracerResult$ProbElemList
  MassShiftElemList <- IsoCombinationsNonTracerResult$MassShiftElemList
  NumberIsoCombEff <- IsoCombinationsNonTracerResult$NumberIsoCombEff
  IsotopesElem <- IsoCombinationsNonTracerResult$IsotopesElem
  
  # Compute tracer impurity combinations
  
  ProbTracerImpurityList <- list()
  MassShiftTracerImpurityList <- list()
  NumberImpurityCombEff <- list()
  
  if (CorrectTracerImpurity) {
    
    if (NumberTracers > 0) {
      
      IsoCombinationsImpurityResult <- IsoCombinationsTracerImpurity(MoleculeFragmentData=MoleculeFragmentData,
                                                                     ElementInfo=ElementInfo,
                                                                     Transitions=Transitions,
                                                                     NumberTransitions=NumberTransitions,
                                                                     IDTracer=IDTracer,
                                                                     LabelPresent=LabelPresent,
                                                                     CalculationThreshold=CalculationThreshold)
      
      ProbTracerImpurityList <- IsoCombinationsImpurityResult$ProbTracerImpurityList
      MassShiftTracerImpurityList <- IsoCombinationsImpurityResult$MassShiftTracerImpurityList
      NumberImpurityCombEff <- IsoCombinationsImpurityResult$NumberImpurityCombEff
      
    } 
  } else {
    ProbTracerImpurityList <- "No Correction for Tracer Impurity performed"
    MassShiftTracerImpurityList <- "No Correction for Tracer Impurity performed"
    NumberImpurityCombEff <- "No Correction for Tracer Impurity performed"
  } #CorrectTracerImpurity=TRUE
  
  IsoCombinationsResult <- list(ProbTracerList = ProbTracerList, ProbElemList = ProbElemList, ProbTracerImpurityList = ProbTracerImpurityList, MassShiftTracerList = MassShiftTracerList, 
                                MassShiftElemList = MassShiftElemList, MassShiftTracerImpurityList = MassShiftTracerImpurityList, NumberIsoCombEff = unlist(NumberIsoCombEff), 
                                NumberImpurityCombEff = NumberImpurityCombEff, NumberIsoCombTracerEff = unlist(NumberIsoCombTracerEff), IsotopesTracer = IsotopesTracer, IsotopesElem = IsotopesElem, 
                                LabelPresent = LabelPresent)
  
  return(IsoCombinationsResult)
  
}  #IsoCombinationsMaster()

# COMPUTE TRACER IMPURITY COMBINATIONS

# Tracer purity is the probability that a tracer element atom in the 
# tracer substrate that should be labelled actually is labelled. 
# E.g.  a 99.0% isotopic purity 1,2-13C-Glucose has a 99.0% chance for 
# each of its carbon atom positions 1 and 2 of containing a 13C. 
# Thus, there is a 1.0% chance for each of those carbons that they are not 13C.
# Consequently, molecules that contain the tracer
# isotope due to metabolic transfer from the tracer substrate and not due to 
# natural abundance have a certain chance of contributing to labelling states
# with less tracer incorporated. This is due to the decrease in mass shift 
# associated with tracer impurity (e.g. 12C instead of 13C at a carbon position).
# In the following, for each labelling state the probability of having a 
# certain number of 'impure' tracer positions as well as the associated 
# (negative) mass shift are is calculated. 
# The maximum amount of impure tracer positions possible depends on the number
# of label present in a given transition.

IsoCombinationsTracerImpurity <- function(MoleculeFragmentData, ElementInfo, Transitions, NumberTransitions,
                                          IDTracer, LabelPresent, CalculationThreshold) { 
  
  ProbTracerImpurityList <- list()
  MassShiftTracerImpurityList <- list()
  NumberImpurityCombEff <- list()
  
  for (TransitionNo in seq_len(NumberTransitions)) {
    
    tmpProbTracerImpurity_Array <- vector()
    tmpMassShiftTracerImpurity_Array <- vector()
    tmpNumberImpurityCombEff <- vector()
    
    for (ImpureTracer in 0:LabelPresent[TransitionNo]) {
      
      ProbImpurity <- choose(LabelPresent[TransitionNo], ImpureTracer) * (ElementInfo[[IDTracer]][["Tracer.purity"]]^(LabelPresent[TransitionNo] - 
                                                                                                                        ImpureTracer)) * ((1 - ElementInfo[[IDTracer]][["Tracer.purity"]])^ImpureTracer)
      
      if (ProbImpurity > CalculationThreshold) {
        tmpProbTracerImpurity_Array[ImpureTracer + 1] <- ProbImpurity
        tmpMassShiftTracerImpurity_Array[ImpureTracer + 1] <- (-1) * ImpureTracer * ElementInfo[[IDTracer]][["Tracer.isotope.mass.shift"]]
        NumberImpurityCombEff[TransitionNo] <- ImpureTracer
      } else {
        ImpureTracer <- LabelPresent + 1
      }  #ProbImpurity > CalculationThreshold
    }  #ImpureTracer
    
    ProbTracerImpurityList[[TransitionNo]] <- tmpProbTracerImpurity_Array
    MassShiftTracerImpurityList[[TransitionNo]] <- tmpMassShiftTracerImpurity_Array
  }  #TransitionNo
  
  names(ProbTracerImpurityList) <- rownames(Transitions)
  names(MassShiftTracerImpurityList) <- rownames(Transitions)
  names(NumberImpurityCombEff) <- rownames(Transitions)
  
  return(list("ProbTracerImpurityList"=ProbTracerImpurityList,
              "MassShiftTracerImpurityList"=MassShiftTracerImpurityList,
              "NumberImpurityCombEff"=NumberImpurityCombEff))
  
}


# COMPUTE ISOTOPE COMBINATIONS OF TRACER ELEMENT

IsoCombinationsTracer <- function(MoleculeFragmentData, Fragment, ElementInfo, Transitions, NumberTracers, CorrectTracerImpurity, CorrectTracerElementCore, CalculationThreshold) {
  ProbTracerList <- list()
  MassShiftTracerList <- list()
  NumberIsoCombTracerEff <- list()
  IsotopesTracer <- list()
  
  NumberTransitions <- nrow(Transitions)
  NegIsoTracer <- MoleculeFragmentData[["NegIsoTracer"]]
  nTracerMax <- MoleculeFragmentData[["nTracerMax"]]
  MaxLabel <- MoleculeFragmentData[["MaxLabel"]]
  
  LabelPresent <- rep(0, nrow(Transitions))
  
  if (NumberTracers > 0) {
    IDTracer <- MoleculeFragmentData[["IDTracer"]]
    Element <- IDTracer
    LabelPresent <- vector()
    nTracer <- vector()
    
    for (TransitionNo in seq_len(NumberTransitions)) {
      
      # x is passed to the daughter function 'IsoCombinations'
      
      x <- TransitionNo
      
      # ElementInfo[[Element]][[2]] is the tracer isotope mass shift.
      # Dividing by this transforms the mass shift of a transition into
      # the number of label present.
      
      LabelPresent[x] <- Transitions[x, Fragment] / ElementInfo[[Element]][[2]]
      
      nTracer[x] <- as.numeric(nTracerMax - LabelPresent[x])
      
      # Do not constrain PlacesToAssign if tracer purity correction is active,
      # as tracer purity produces negative mass shifts which leads to similar
      # effects as NegIsoTracer.
      
      if (CorrectTracerImpurity || NegIsoTracer == 1) {
        PlacesToAssign <- as.numeric(ElementInfo[[Element]][[2]]) * MaxLabel
      } else {
        PlacesToAssign <- as.numeric(ElementInfo[[Element]][[2]] * (MaxLabel - LabelPresent[x]))
      }
      
      # Correct molecule core or not
      
      if (CorrectTracerElementCore) {
        AvailablePlacesTotal <- nTracer[x]
      } else if (!CorrectTracerElementCore) {
        AvailablePlacesTotal <- as.numeric(nTracerMax) - as.numeric(MaxLabel)
      } else {
        stop(date(), " :: Invalid value for `CorrectTracerElementCore`: ", CorrectTracerElementCore)
      }
      
      # Definition of the IsoCluster variable that is transferred to the
      # daughter function IsoCombinations() that computes the actual
      # isotope combinations
      
      if (PlacesToAssign <= AvailablePlacesTotal) {
        IsoCluster <- PlacesToAssign
      } else {
        IsoCluster <- AvailablePlacesTotal
      }
      
      IsoCombinationsResultList <- IsoCombinations(Fragment = Fragment, x = x, Element = Element, ElementInfo = ElementInfo, AvailablePlacesTotal = AvailablePlacesTotal, IsoCluster = IsoCluster, CalculationThreshold = CalculationThreshold)
      NumberIsoComb <- IsoCombinationsResultList[["NumberIsoComb"]]
      tmpProbTracer_Array <- vector()
      tmpMassShiftTracer_Array <- vector()
      tmpIsotopesTracer <- list()
      for (i in seq_len(NumberIsoComb)) {
        tmpProbTracer_Array[i] <- IsoCombinationsResultList[["ProbArray"]][[i]]
        tmpMassShiftTracer_Array[i] <- IsoCombinationsResultList[["MassShiftArray"]][[i]]
      } # i
      
      ProbTracerList[[x]] <- tmpProbTracer_Array
      MassShiftTracerList[[x]] <- tmpMassShiftTracer_Array
      NumberIsoCombTracerEff[[x]] <- NumberIsoComb
      IsotopesTracer[[x]] <- IsoCombinationsResultList[["Isotopes"]]
    } # TransitionsNo (represented by x)
    
    names(ProbTracerList) <- rownames(Transitions)
    names(MassShiftTracerList) <- rownames(Transitions)
    names(NumberIsoCombTracerEff) <- rownames(Transitions)
    names(IsotopesTracer) <- rownames(Transitions)
  } else {
    LabelPresent <- rep(0, NumberTransitions)
  }
  
  return(list(
    "ProbTracerList" = ProbTracerList,
    "MassShiftTracerList" = MassShiftTracerList,
    "NumberIsoCombTracerEff" = NumberIsoCombTracerEff,
    "IsotopesTracer" = IsotopesTracer,
    "LabelPresent" = LabelPresent
  ))
}



#' @importFrom dplyr last
#' @importFrom dplyr arrange
#' @importFrom stringr str_c
#' @importFrom tibble as_tibble
#' @importFrom magrittr '%>%'

IsoCombinations <- function(Fragment, x, Element, ElementInfo, AvailablePlacesTotal, IsoCluster, CalculationThreshold) {
  
  # NumberIso is the number of isotopes per element.
  # The number of nested loops used depends on the number of
  # isotopes of the current element.
  
  NumberIso <- nrow(ElementInfo[[Element]][["Isotopes"]][[1]])
  
  Combinations <- vector()
  Combinations[1] <- 1
  
  # Calculation of Combinations[IsotopeNo+1], the number of isotope combinations
  # when considering IsotopeNo different isotopes (also including combinations
  # where the sum of isotopes exceeds the number of available places IsoCluster).
  
  for (IsotopeNo in seq_len(NumberIso)) {
    Combinations[IsotopeNo + 1] <- Combinations[IsotopeNo] * ((IsoCluster) + 1)
  } # IsotopeNo
  
  TotalCombinations <- dplyr::last(Combinations) # get last element
  
  # Computation of the isotope combinations list `CombinationsArrayPreList`.
  # At this point, combinations where the sum of the different isotopes exceeds
  # IsoCluster are still considered, leading to invalid combinations that are
  # later removed. In its columns, CombinationsArrayPreList contains the number
  # of occurences of a given isotope in a given combination.
  # Each row corresponds to a combination.
  
  CombinationsArrayPreList <- list()
  
  for (IsotopeNo in seq_len(NumberIso)) {
    p <- 1
    tmpCombinationsArrayPre <- vector()
    
    # first, for each column (isotope) the repetitive column element is computed
    
    for (IsoCount in 0:IsoCluster) {
      tmpCombinationsArrayPre[p:(p + Combinations[IsotopeNo] - 1)] <- IsoCount
      p <- p + Combinations[IsotopeNo]
    } # IsoCount
    
    # Then, the repetitive column element is repeatedly inserted until the
    # number of rows matches `TotalCombinations`
    
    for (Repeat in seq_len((TotalCombinations / Combinations[IsotopeNo + 1])) - 1) {
      if (Repeat > 0) {
        tmpCombinationsArrayPre[p:(p + Combinations[IsotopeNo + 1] - 1)] <- tmpCombinationsArrayPre[seq_len(Combinations[IsotopeNo + 1])]
        p <- p + Combinations[IsotopeNo + 1]
      }
    } # Repeat
    
    CombinationsArrayPreList[[IsotopeNo]] <- tmpCombinationsArrayPre
  } # IsotopeNo
  names(CombinationsArrayPreList) <- stringr::str_c("Isotope_", seq_len(NumberIso))
  
  # CombinationsArrayPreList is sorted to get the most probable combinations on top
  
  CombinationsArrayPreList2 <- CombinationsArrayPreList %>% as_tibble() %>% dplyr::arrange(!!!syms(names(CombinationsArrayPreList))) %>% as.list()
  
  CombinationsArrayPre.df <- as.data.frame(lapply(CombinationsArrayPreList2, function(x) rev(x)))
  
  RowSum <- apply(CombinationsArrayPre.df, 1, sum)
  j <- 1
  
  # CombinationsArray is derived from CombinationsArrayPreList by removing
  # all combinations that do not have a row sum of IsoCluster
  
  CombinationsArray <- NULL
  
  for (i in seq_len(TotalCombinations)) {
    if (RowSum[i] == IsoCluster) {
      tmpA <- vector()
      
      tmpA <- as.numeric(CombinationsArrayPre.df[i, seq_len(NumberIso)])
      tmpA[1] <- tmpA[1] + (AvailablePlacesTotal - IsoCluster)
      CombinationsArray <- rbind(CombinationsArray, tmpA)
      
      j <- j + 1
    } # RowSum[i]==IsoCluster
  } # i
  
  # The number of combinations in CombinationsArray is derived,
  # CombinationsArray is stored in Isotopes
  
  NumberIsoCombMax <- nrow(CombinationsArray)
  
  Isotopes <- CombinationsArray
  
  # The probabilities and mass shifts of the isotope combinations in
  # CombinationsArray are computed using IsoCombinationsProbMassShift()
  
  IsoComb <- 1
  ProbArray.vec <- vector()
  MassShiftArray.vec <- vector()
  NumberIsoComb.vec <- vector()
  
  for (IsoCombPre in seq_len(NumberIsoCombMax)) {
    icpms <- IsoCombinationsProbMassShift(
      Element = Element, NumberIso = NumberIso, ElementInfo = ElementInfo, CombinationsArray = CombinationsArray,
      IsoCombPre = IsoCombPre, AvailablePlacesTotal = AvailablePlacesTotal
    )
    
    # Due to CalculationThreshold which limits calculations to probability
    # values large enough to be relevant, probability and mass shift are
    # usually not computed for all isotope combinations.
    
    if (icpms[["ProbIsoComb"]] > CalculationThreshold) {
      ProbArray.vec[IsoComb] <- icpms[["ProbIsoComb"]]
      MassShiftArray.vec[IsoComb] <- icpms[["MassShiftIsoComb"]]
      NumberIsoComb <- IsoComb
      IsoComb <- IsoComb + 1
    } else {
      if (IsoComb == 1) {
        NumberIsoComb <- 0
      }
      
      IsoCombPre <- NumberIsoCombMax + 1
    }
  } # IsoCombPre
  
  return(list(Isotopes = Isotopes, ProbArray = ProbArray.vec, MassShiftArray = MassShiftArray.vec, NumberIsoComb = NumberIsoComb))
} # IsoCombinations

# Calculation of probability and mass shift associated with an isotope combination

IsoCombinationsProbMassShift <- function(Element, NumberIso, ElementInfo, CombinationsArray, IsoCombPre, AvailablePlacesTotal) {
  
  ProbFactorialTermLog <- lfactorial(AvailablePlacesTotal)
  ProbExpTerm <- 1
  MassShiftIsoComb <- 0
  ElementTable <- data.frame(ElementInfo[[Element]][[1]])
  
  for (IsotopeNo in seq_len(NumberIso)) {
    
    IsotopeCount <- as.numeric(CombinationsArray[IsoCombPre, IsotopeNo])
    IsotopeAbundance <- ElementTable[IsotopeNo, 1]
    IsotopeMassShift <- ElementTable[IsotopeNo, 2]
    
    # Probability factorial term (log factorial as factorial can lead to numbers that are too big)
    ProbFactorialTermLog <- ProbFactorialTermLog - lfactorial(IsotopeCount)
    
    # Probability exponential term
    ProbExpTerm <- ProbExpTerm * (IsotopeAbundance^IsotopeCount)
    
    # mass shift
    MassShiftIsoComb <- MassShiftIsoComb + IsotopeMassShift * IsotopeCount
  }  #IsotopeNo
  
  ProbIsoComb <- exp(ProbFactorialTermLog) * ProbExpTerm
  
  return(list(ProbIsoComb = ProbIsoComb, MassShiftIsoComb = MassShiftIsoComb))
  
}


# COMPUTE ISOTOPE COMBINATIONS OF NON-TRACER ELEMENT

# Combinations, probabilites and mass shifts of the non-tracer elements are
# computed, given that not only the tracer element is provided in
# the molcule file. As with the tracer element,
# IsoCombinationsMaster() feeds into IsoCombinations() also for the other elements.

IsoCombinationsNonTracer <- function(MoleculeFragmentData, Fragment, ElementInfo,
                                     NumberTracers, IDTracer, CalculationThreshold) {
  ProbElemList <- list()
  MassShiftElemList <- list()
  NumberIsoCombEff <- list()
  IsotopesElem <- list()
  
  ElementNonTracer_label <- "NonTracer"
  ElementsNonTracer <- MoleculeFragmentData[[ElementNonTracer_label]]
  
  if (length(ElementsNonTracer) == 0) {
    ElementsNonTracer <- MoleculeFragmentData[["ZeroTracer"]]
    ElementNonTracer_label <- "ZeroTracer"
  }
  
  NumberElementsNonTracer <- length(ElementsNonTracer)
  
  if (NumberElementsNonTracer > 0) {
    for (NonTracer in seq_len(NumberElementsNonTracer)) {
      
      # Definition of parameters for the non-tracer elements
      x <- NonTracer
      
      Element <- names(ElementsNonTracer)[x]
      
      if (NumberTracers > 0) {
        PlacesToAssign <- as.numeric(ElementInfo[[IDTracer]][[2]] * MoleculeFragmentData[["MaxLabel"]])
      } else {
        PlacesToAssign <- 0
      } # NumberTracers > 0
      
      AvailablePlacesTotal <- as.numeric(MoleculeFragmentData[[ElementNonTracer_label]][Element])
      
      if (PlacesToAssign <= AvailablePlacesTotal) {
        IsoCluster <- PlacesToAssign
      } else {
        IsoCluster <- AvailablePlacesTotal
      }
      
      # Calling of IsoCombinations()
      IsoCombinationsResultList <- IsoCombinations(
        Fragment = Fragment, x = x, Element = Element, ElementInfo = ElementInfo,
        AvailablePlacesTotal = AvailablePlacesTotal, IsoCluster = IsoCluster, CalculationThreshold = CalculationThreshold
      )
      NumberIsoComb <- IsoCombinationsResultList[["NumberIsoComb"]]
      tmpProbElem_Array <- vector()
      tmpMassShiftElem_Array <- vector()
      
      # tmpIsotopesElem <- list()
      
      for (i in seq_len(NumberIsoComb)) {
        tmpProbElem_Array[i] <- IsoCombinationsResultList[["ProbArray"]][i]
        tmpMassShiftElem_Array[i] <- IsoCombinationsResultList[["MassShiftArray"]][i]
      } # i
      
      ProbElemList[[x]] <- tmpProbElem_Array
      MassShiftElemList[[x]] <- tmpMassShiftElem_Array
      NumberIsoCombEff[[x]] <- NumberIsoComb
      IsotopesElem[[x]] <- IsoCombinationsResultList[["Isotopes"]]
    } # NonTracer
    
    names(ProbElemList) <- names(ElementsNonTracer)
    names(MassShiftElemList) <- names(ElementsNonTracer)
    names(NumberIsoCombEff) <- names(ElementsNonTracer)
    names(IsotopesElem) <- names(ElementsNonTracer)
  } # NumberElementsNonTracer > 0
  
  return(list(
    "ProbElemList" = ProbElemList,
    "MassShiftElemList" = MassShiftElemList,
    "NumberIsoCombEff" = NumberIsoCombEff,
    "IsotopesElem" = IsotopesElem
  ))
}

# Calculate probabilites and mass shifts of element combinations

# With the probabilities and mass shifts of the different
# isotope combinations calculated, the next step is to combine
# the isotope combination probabilities and mass shifts of
# the different elements (and the tracer impurity, if checked) using the function
# 'ElementCombinations'. For a given molecule(-fragment) the function combinatorically
# multiplies all isotope combination probabilities of all non-tracer elements,
# the tracer element and the tracer impurity with each other
# In the same way, the isotope combination mass shifts are added to each
# other. In the end, the probability that a molecule(-fragment) contains a combination
# of specific element isotope combinations for its various elements
# (and a certain number of 'impure' tracer atoms) is yielded,
# together with the mass shift associated with
# such a combination. Because combinatorically multiplying all
# isotope combinations is resource intensive, the
# calculation stops at certain probability thresholds defined by the
# parameter CalculationThreshold. In those cases the probability of a combination is so low
# that it does not affect the correction.

ElementCombinations <- function(MoleculeInfo, MoleculeNo, Fragment, IsoCombinationsResult, CalculationThreshold,
                                CorrectTracerImpurity) {
  
  MoleculeFragmentData <- MoleculeInfo[[MoleculeNo]][[stringr::str_c("Fragment_", Fragment)]]
  ProbTracerList <- IsoCombinationsResult[["ProbTracerList"]]
  MassShiftTracerList <- IsoCombinationsResult[["MassShiftTracerList"]]
  NumberIsoCombTracerEff <- IsoCombinationsResult[["NumberIsoCombTracerEff"]]
  ProbElemList <- IsoCombinationsResult[["ProbElemList"]]
  MassShiftElemList <- IsoCombinationsResult[["MassShiftElemList"]]
  NumberIsoCombEff <- IsoCombinationsResult[["NumberIsoCombEff"]]
  ProbTracerImpurityList <- IsoCombinationsResult[["ProbTracerImpurityList"]]
  MassShiftTracerImpurityList <- IsoCombinationsResult[["MassShiftTracerImpurityList"]]
  NumberImpurityCombEff <- IsoCombinationsResult[["NumberImpurityCombEff"]]
  ElementsNonTracer <- MoleculeFragmentData[["NonTracer"]]
  
  if (length(ElementsNonTracer) == 0) {
    ElementsNonTracer <- MoleculeFragmentData[["ZeroTracer"]]
  }
  
  NumberElementsNonTracer <- length(ElementsNonTracer)
  
  # CALCULATION OF THE NON-TRACER ELEMENT COMBINATIONS
  
  ElementCombinationsNonTracerResult <- ElementCombinationsNonTracer(ProbElemList=ProbElemList,
                                                                     MassShiftElemList=MassShiftElemList,
                                                                     NumberIsoCombEff=NumberIsoCombEff,
                                                                     ElementsNonTracer=ElementsNonTracer,
                                                                     NumberElementsNonTracer = NumberElementsNonTracer,
                                                                     CalculationThreshold=CalculationThreshold)
  
  CombElemProbList <- ElementCombinationsNonTracerResult$CombElemProbList
  CombElemMassShiftList <- ElementCombinationsNonTracerResult$CombElemMassShiftList
  NumberEleComb <- ElementCombinationsNonTracerResult$NumberEleComb
  
  # CALCULATION OF THE TRACER ELEMENT - NON-TRACER ELEMENT COMBINATIONS
  
  NumberTracers <- length(MoleculeFragmentData[["Tracer"]])
  Transitions <- MoleculeInfo[[MoleculeNo]][["Transitions"]]
  NumberTransitions <- nrow(Transitions)
  
  ElementCombinationsTracerResult <- ElementCombinationsTracer(ProbTracerList=ProbTracerList,
                                                               MassShiftTracerList=MassShiftTracerList, 
                                                               NumberIsoCombTracerEff=NumberIsoCombTracerEff, 
                                                               CombElemProbList=CombElemProbList, 
                                                               CombElemMassShiftList=CombElemMassShiftList,
                                                               NumberEleComb=NumberEleComb,
                                                               Fragment=Fragment,
                                                               NumberElementsNonTracer=NumberElementsNonTracer,
                                                               NumberTracers=NumberTracers, 
                                                               Transitions=Transitions,
                                                               NumberTransitions=NumberTransitions, 
                                                               CalculationThreshold=CalculationThreshold)
  
  TracerElemProbList <- ElementCombinationsTracerResult$TracerElemProbList
  TracerElemMassShiftList <- ElementCombinationsTracerResult$TracerElemMassShiftList
  NumberTracerElemComb <- ElementCombinationsTracerResult$NumberTracerElemComb
  
  # COMBINATION WITH TRACER IMPURITY STATES
  
  if (CorrectTracerImpurity && (NumberTracers > 0)) {
    
    ElementCombinationsTracerImpurityResult <- ElementCombinationsTracerImpurity(ProbTracerImpurityList=ProbTracerImpurityList,
                                                                                 MassShiftTracerImpurityList=MassShiftTracerImpurityList,
                                                                                 NumberImpurityCombEff=NumberImpurityCombEff,
                                                                                 TracerElemProbList=TracerElemProbList,
                                                                                 TracerElemMassShiftList=TracerElemMassShiftList,
                                                                                 NumberTracerElemComb=NumberTracerElemComb,
                                                                                 Transitions=Transitions,
                                                                                 NumberTransitions=NumberTransitions,
                                                                                 CalculationThreshold=CalculationThreshold) 
    
    TracerImpurityCombProbList <- ElementCombinationsTracerImpurityResult$TracerImpurityCombProbList
    TracerImpurityCombMassShiftList <- ElementCombinationsTracerImpurityResult$TracerImpurityCombMassShiftList
    NumberTracerImpurityComb <- ElementCombinationsTracerImpurityResult$NumberTracerImpurityComb
    
  } else {
    
    TracerImpurityCombProbList <- list()
    TracerImpurityCombMassShiftList <- list()
    NumberTracerImpurityComb <- vector()
    
  } 
  
  return(list(TracerElemProbList = TracerElemProbList, TracerElemMassShiftList = TracerElemMassShiftList, TracerImpurityCombProbList = TracerImpurityCombProbList, 
              TracerImpurityCombMassShiftList = TracerImpurityCombMassShiftList, NumberTracerImpurityComb = NumberTracerImpurityComb, NumberTracerElemComb = NumberTracerElemComb))
  
}  #ElementCombinations()

# CALCULATION OF THE NON-TRACER ELEMENT COMBINATIONS 

ElementCombinationsNonTracer <- function(ProbElemList, MassShiftElemList, NumberIsoCombEff, 
                                         ElementsNonTracer, NumberElementsNonTracer, CalculationThreshold) {
  
  # If the number of non-tracer elements is > 1, all isotope combination probabilities of non-tracer
  # elements 1 and 2 of a molecule(-fragment) are multiplied with each other and the corresponding mass shifts are added.  Like in 'IsoCombinations' there
  # is a probability threshold which resets to the outer loop if the probability calculated is below.
  
  CombElemProbList <- list()
  tmpCombElemProb_Array <- vector()
  
  CombElemMassShiftList <- list()
  tmpCombElemMassShift_Array <- vector()
  
  NumberEleComb <- list()
  tmpNumberEleComb <- vector()
  
  if (NumberElementsNonTracer > 1) {
    
    EleComb <- 1
    
    for (IsoCombinationA in seq_len(as.numeric(NumberIsoCombEff[1]))) {
      
      for (IsoCombinationB in seq_len(as.numeric(NumberIsoCombEff[[2]]))) {
        
        tmpCombElemProb <- ProbElemList[[1]][IsoCombinationA] * ProbElemList[[2]][IsoCombinationB]
        
        if (tmpCombElemProb > CalculationThreshold) {
          
          tmpCombElemProb_Array[EleComb] <- tmpCombElemProb
          tmpCombElemMassShift_Array[EleComb] <- MassShiftElemList[[1]][IsoCombinationA] + MassShiftElemList[[2]][IsoCombinationB]
          tmpNumberEleComb <- EleComb
          
          EleComb <- EleComb + 1
        } else {
          IsoCombinationB <- NumberIsoCombEff[[2]] + 1
        }  #tmpCombElemProb > CalculationThreshold
      }  #IsoCombinationB
    }  #IsoCombinationA
    
    CombElemProbList[[2]] <- tmpCombElemProb_Array
    CombElemMassShiftList[[2]] <- tmpCombElemMassShift_Array
    NumberEleComb[[2]] <- tmpNumberEleComb
    
    rm(tmpCombElemProb_Array)
    rm(tmpCombElemMassShift_Array)
    rm(tmpNumberEleComb)
    
  } else if (NumberElementsNonTracer > 0) {
    
    # If the number of non-tracer elements is 1, the probability and mass shift arrays are equal to those generated in 'IsoCombinationsResult'.
    
    tmpCombElemProb_Array[seq_len(NumberIsoCombEff[[1]])] <- ProbElemList[[1]][seq_len(NumberIsoCombEff[[1]])]
    tmpCombElemMassShift_Array[seq_len(NumberIsoCombEff[[1]])] <- MassShiftElemList[[1]][seq_len(NumberIsoCombEff[[1]])]
    tmpNumberEleComb <- NumberIsoCombEff[[1]]
    
    CombElemProbList[[1]] <- tmpCombElemProb_Array
    CombElemMassShiftList[[1]] <- tmpCombElemMassShift_Array
    NumberEleComb[[1]] <- tmpNumberEleComb
    
    rm(tmpCombElemProb_Array)
    rm(tmpCombElemMassShift_Array)
    rm(tmpNumberEleComb)
    
  }  #NumberElementsNonTracer>1
  
  # If there are more than 2 non-tracer elements, the element combinations are calculated in iterations of the following outer loop until all elements have
  # been covered. The algorithm uses the CombElemProbList covering the previously combined element probabilities and multiplies these with the isotope
  # combination probabilities of the next element (ProbElemList).  This is fed back into CombElemProbList at the current element index of the loop). The
  # same is done additively for the mass shift array CombElemMassShiftList.
  
  if (NumberElementsNonTracer > 2) {
    
    tmpCombElemProb_Array <- vector()
    tmpCombElemMassShift_Array <- vector()
    tmpNumberEleComb <- vector()
    
    for (NonTracer in 3:NumberElementsNonTracer) {
      EleComb <- 1
      for (ElementCombination in seq_len(NumberEleComb[[NonTracer - 1]])) {
        for (IsoCombination in seq_len(NumberIsoCombEff[[NonTracer]])) {
          
          tmpCombElemProb <- CombElemProbList[[NonTracer - 1]][[ElementCombination]] * ProbElemList[[NonTracer]][[IsoCombination]]
          if (tmpCombElemProb > CalculationThreshold) {
            
            tmpCombElemProb_Array[EleComb] <- tmpCombElemProb
            tmpCombElemMassShift_Array[EleComb] <- CombElemMassShiftList[[NonTracer - 1]][[ElementCombination]] + MassShiftElemList[[NonTracer]][[IsoCombination]]
            tmpNumberEleComb <- EleComb
            EleComb <- EleComb + 1
            
          } else {
            
            IsoCombination <- NumberIsoCombEff[[NonTracer]] + 1
            
          }  #tmpCombElemProb > CalculationThreshold
        }  #IsoCombination
        
      }  #ElementCombination
      CombElemProbList[[NonTracer]] <- tmpCombElemProb_Array
      CombElemMassShiftList[[NonTracer]] <- tmpCombElemMassShift_Array
      NumberEleComb[[NonTracer]] <- tmpNumberEleComb
      
    }  #NonTracer
  }  #NumberElementsNonTracer>2
  
  return(list("CombElemProbList"=CombElemProbList,
              "CombElemMassShiftList"=CombElemMassShiftList,
              "NumberEleComb"=NumberEleComb))
  
}

# CALCULATION OF THE TRACER ELEMENT - NON-TRACER ELEMENT COMBINATIONS

# Given that there are more elements than just the tracer element, the element combinations from CombElemProbList are now multiplied with the isotope
# combinations of the tracer element, for all labelling states. This is analogous to the previous algorithms from this function.  However, in the
# computation of the mass shift array, the intrinsic mass shift (due to tracer incorporation) associated with the respective labelling state is added in
# addtion.  In the end, the probability array TracerElemProbList and the mass shift array TracerElemMassShiftList are derived. They contain the
# probabilities and mass shifts of all relevant elemental combinations of isotope combinations for all labelling states of a given molecule(-fragment).

ElementCombinationsTracer <- function(ProbTracerList, MassShiftTracerList, NumberIsoCombTracerEff, 
                                      CombElemProbList, CombElemMassShiftList, NumberEleComb,
                                      Fragment, NumberElementsNonTracer, NumberTracers, Transitions,
                                      NumberTransitions, CalculationThreshold) {
  
  TracerElemProbList <- list()
  TracerElemMassShiftList <- list()
  NumberTracerElemComb <- vector()
  
  if (NumberTracers > 0) {
    
    if (NumberElementsNonTracer > 0) {
      
      for (TransitionNo in seq_len(NumberTransitions)) {
        
        tmpTracerElemProb_Array <- vector()
        tmpTracerElemMassShift_Array <- vector()
        EleComb <- 1
        
        for (ElementCombination in seq_len(NumberEleComb[[NumberElementsNonTracer]])) {
          
          for (IsoCombination in seq_len(NumberIsoCombTracerEff[[TransitionNo]])) {
            
            TracerElemProb <- CombElemProbList[[NumberElementsNonTracer]][[ElementCombination]] * ProbTracerList[[TransitionNo]][[IsoCombination]]
            if (TracerElemProb > CalculationThreshold) {
              tmpTracerElemProb_Array[EleComb] <- TracerElemProb
              
              tmpTracerElemMassShift_Array[EleComb] <- CombElemMassShiftList[[NumberElementsNonTracer]][ElementCombination] + MassShiftTracerList[[TransitionNo]][IsoCombination] + 
                Transitions[TransitionNo, Fragment]
              
              NumberTracerElemComb[TransitionNo] <- EleComb
              
              EleComb <- EleComb + 1
              
            } else {
              IsoCombination <- NumberIsoCombTracerEff[[TransitionNo]] + 1
            }  #TracerElemProb > CalculationThreshold
          }  #IsoCombination
        }  #ElementCombination
        
        TracerElemProbList[[TransitionNo]] <- tmpTracerElemProb_Array
        TracerElemMassShiftList[[TransitionNo]] <- tmpTracerElemMassShift_Array
        
      }  #TransitionNo
    } else {
      
      # If the only element considered is the tracer element, the probability- and mass shift arrays from 'IsoCombinationsResult' are directly fed into the
      # probability array TracerElemProbList and the mass shift array TracerElemMassShiftList.  However, for each transition its intrinsic mass shift is
      # added to the mass shift array in addition.
      
      for (i in seq_len(NumberTransitions)) {
        NumberTracerElemComb[i] <- NumberIsoCombTracerEff[[i]]
      }  #i
      
      for (TransitionNo in seq_len(NumberTransitions)) {
        tmpTracerElemProb_Array <- vector()
        tmpTracerElemMassShift_Array <- vector()
        
        for (IsoCombination in seq_len(NumberIsoCombTracerEff[[TransitionNo]])) {
          
          tmpTracerElemProb_Array[IsoCombination] <- ProbTracerList[[TransitionNo]][IsoCombination]
          
          tmpTracerElemMassShift_Array[IsoCombination] <- MassShiftTracerList[[TransitionNo]][IsoCombination] + Transitions[TransitionNo, Fragment]
          
        }  #IsoCombination
        TracerElemProbList[[TransitionNo]] <- tmpTracerElemProb_Array
        TracerElemMassShiftList[[TransitionNo]] <- tmpTracerElemMassShift_Array
      }  #TransitionNo
    }  #NumberElementsNonTracer>0
  } else {
    
    # If the molecule(-fragment) in question contains no tracer element, the probability array TracerElemProbList and the mass shift array
    # TracerElemMassShiftList are equal to CombElemProbList and CombElemMassShiftList for all transitions of this molecule(-fragment).
    
    for (TransitionNo in seq_len(NumberTransitions)) {
      
      tmpTracerElemProb_Array <- vector()
      tmpTracerElemMassShift_Array <- vector()
      
      for (i in seq_len(NumberEleComb[[NumberElementsNonTracer]])) {
        
        tmpTracerElemProb_Array[i] <- CombElemProbList[[NumberElementsNonTracer]][i]
        tmpTracerElemMassShift_Array[i] <- CombElemMassShiftList[[NumberElementsNonTracer]][i]
      }  #i
      
      TracerElemProbList[[TransitionNo]] <- tmpTracerElemProb_Array
      TracerElemMassShiftList[[TransitionNo]] <- tmpTracerElemMassShift_Array
      NumberTracerElemComb[TransitionNo] <- NumberEleComb[[NumberElementsNonTracer]]
    }  #TransitionNo
    
  }  #NumberTracers>0
  
  names(TracerElemProbList) <- rownames(Transitions)
  names(TracerElemMassShiftList) <- rownames(Transitions)
  names(NumberTracerElemComb) <- rownames(Transitions)
  
  return(list("TracerElemProbList"=TracerElemProbList,
              "TracerElemMassShiftList"=TracerElemMassShiftList,
              "NumberTracerElemComb"=NumberTracerElemComb))
  
}

# COMBINATION WITH TRACER IMPURITY STATES

# If correction for tracer impurity is switched on, the arrays TracerElemProbList and TracerElemMassShiftList are additonally combined with the
# probabilites and mass shifts associated with different numbers of 'impure' tracer atoms in the molecule(-fragment). This yields the arrays
# TracerImpurityCombProbList and TracerImpurityCombMassShiftList. They give the probability and mass shift associated with defined elemental isotope
# combinations and a defined number of 'impure' tracer atoms.

ElementCombinationsTracerImpurity <- function(ProbTracerImpurityList,MassShiftTracerImpurityList,NumberImpurityCombEff,
                                              TracerElemProbList,TracerElemMassShiftList,NumberTracerElemComb,
                                              Transitions,NumberTransitions,CalculationThreshold) {
  
  TracerImpurityCombProbList <- list()
  TracerImpurityCombMassShiftList <- list()
  NumberTracerImpurityComb <- vector()
  
  for (TransitionNo in seq_len(NumberTransitions)) {
    
    tmpTracerImpurityCombProb_Array <- vector()
    tmpTracerImpurityCombMassShift_Array <- vector()
    ImpurityComb <- 1
    for (ImpureTracer in 0:NumberImpurityCombEff[[TransitionNo]]) {
      
      for (EleComb in seq_len(NumberTracerElemComb[TransitionNo])) {
        
        TracerImpurityCombProb <- ProbTracerImpurityList[[TransitionNo]][ImpureTracer + 1] * TracerElemProbList[[TransitionNo]][EleComb]
        
        if (TracerImpurityCombProb > CalculationThreshold) {
          
          tmpTracerImpurityCombProb_Array[ImpurityComb] <- TracerImpurityCombProb
          tmpTracerImpurityCombMassShift_Array[ImpurityComb] <- MassShiftTracerImpurityList[[TransitionNo]][ImpureTracer + 1] + TracerElemMassShiftList[[TransitionNo]][EleComb]
          
          NumberTracerImpurityComb[TransitionNo] <- ImpurityComb
          ImpurityComb <- ImpurityComb + 1
        } else {
          EleComb <- NumberTracerElemComb[TransitionNo] + 1
        }  #TracerImpurityCombProb>CalculationThreshold
      }  #EleComb
    }  #ImpureTracer
    TracerImpurityCombProbList[[TransitionNo]] <- tmpTracerImpurityCombProb_Array
    TracerImpurityCombMassShiftList[[TransitionNo]] <- tmpTracerImpurityCombMassShift_Array
  }  #TransitionNo
  
  names(TracerImpurityCombProbList) <- rownames(Transitions)
  names(TracerImpurityCombMassShiftList) <- rownames(Transitions)
  names(NumberTracerImpurityComb) <- rownames(Transitions)
  
  return(list("TracerImpurityCombProbList"=TracerImpurityCombProbList,
              "TracerImpurityCombMassShiftList"=TracerImpurityCombMassShiftList,
              "NumberTracerImpurityComb"=NumberTracerImpurityComb))
  
}

# Calculates the sum of probabilities with a common mass shift

# The function 'MassShiftProbabilities' yields the
# probability that a given transition of a
# molecule(-fragment) (through naturally occuring isotopes)
# produces a certain mass shift in relation to the same molecule(-fragment)
# containing no isotopes of higher mass. This probability is gained
# by summing all entries of a molecule(-fragment) of a given
# transition in TracerElemProbList that have the same mass shift in
# TracerElemMassShiftList. If correction for tracer impurity is turned
# on and if the molecule(-fragment) in question can contain tracer,
# 'MassShiftProbabilities' uses TracerImpurityCombProbList
# and TracerImpurityCombMassShiftList instead.

#' @importFrom stringr str_c
#' 
MassShiftProbabilities <- function(MoleculeInfo, MoleculeNo, Fragment, MaxMassShift, ElementCombinations, CorrectTracerImpurity) {
  
  MoleculeData <- MoleculeInfo[[MoleculeNo]][[stringr::str_c("Fragment_", Fragment)]]
  
  Transitions <- MoleculeInfo[[MoleculeNo]][["Transitions"]]
  NumberTransitions <- nrow(Transitions)
  NumberTracers <- length(MoleculeData[["Tracer"]])
  
  if (CorrectTracerImpurity) {
    if (NumberTracers > 0) {
      
      ProbabilityList <- ElementCombinations[["TracerImpurityCombProbList"]]
      MassShiftList <- ElementCombinations[["TracerImpurityCombMassShiftList"]]
      
    } else {
      ProbabilityList <- ElementCombinations[["TracerElemProbList"]]
      MassShiftList <- ElementCombinations[["TracerElemMassShiftList"]]
      
    }
    
  } else {
    
    ProbabilityList <- ElementCombinations[["TracerElemProbList"]]
    MassShiftList <- ElementCombinations[["TracerElemMassShiftList"]]
    
  }
  
  # This loop goes through all mass shifts until it reaches the maximum possible
  # mass shift of the current molecule(-fragment).  For each labelling state,
  # it finds indexes in MassShiftList associated with the given mass shift 
  # and uses them to find the corresponding probabilites in ProbabilityList.
  # Those values are eventually summed up to yield the probability that a 
  # certain labelling state produces a certain mass shift in relation to the
  # unlabelled species due to natural abundance and possibly tracer impurity.
  
  CumProbList <- list()
  
  for (MassShift in 0:MaxMassShift) {
    
    tmpCumProb.vec <- vector()
    for (TransitionNo in seq_len(NumberTransitions)) {
      
      tmpIdx <- which(MassShiftList[[TransitionNo]] == MassShift)
      
      tmpCumProb.vec[TransitionNo] <- sum(ProbabilityList[[TransitionNo]][tmpIdx])
      
    }
    CumProbList[[MassShift + 1]] <- tmpCumProb.vec
  }  #MassShift
  
  names(CumProbList) <- stringr::str_c("MassShift_", 0:MaxMassShift)
  
  return(CumProbList)
}


# Generation of a mass shift vector of labelling states

# The purpose of the function 'MassShiftTransitions' is
# to generate a vector (FragmentMassShift) which contains
# the mass shift of a molecule(-fragment) of a given
# transition in relation to the same molecule(-fragment)
# containing no isotopes of higher mass. This is required to
# be able to yield the probability that a certain transition
# or full MS ion x is derived from another transition/full MS ion
# y due to natural stable isotope abundance. This probability is yielded
# by matching the mass shift k of transition x in FragmentMassShift
# with the mass shift k of transition y in CumProbList
# (this is the mass shift 'produced' by y through natural abundance).
# The probability is found associated with mass shift k in CumProbList.

MassShiftTransitions <- function(MoleculeInfo, MoleculeNo, Fragment) {
  
  Transitions <- MoleculeInfo[[MoleculeNo]][["Transitions"]]
  NumberTransitions <- nrow(Transitions)
  
  FragmentMassShift <- vector()
  
  for (TransitionNo in seq_len(NumberTransitions)) {
    FragmentMassShift[TransitionNo] <- Transitions[TransitionNo, Fragment]
  }  #TransitionNo
  
  return(FragmentMassShift)
  
}  #MassShiftTransitions


# Calculation of probability that a given labelling state produces the mass shift of another labelling state

# In the function 'FragmentsCombinedProbability' the final
# probability matrix CombFragmentProb is produced. In this matrix, for each
# transition y of a molecule the probability that it produces a transition x
# of higher mass through natural isotope abundance is given.
# This is done by matching the mass shifts from
# FragmentMassShift and CumProbList, as described before for the
# function 'MassShiftTransitions'. For MS/MS data, the
# probability that the product ion y produces the product ion
# mass shift of transition x and and the probability that neutral
# loss y produces the neutral loss mass shift of transition x
# are multiplied. This way the overall probability that
# transition y produces transition x through natural isotope abundance
# is yielded.

#' @importFrom stringr str_detect
#' @importFrom magrittr '%>%'
#' 
FragmentsCombinedProbability <- function(MoleculeInfo, MoleculeNo, CumProbList, FragmentMassShift) {
  
  NumberFragments <- stringr::str_detect(names(MoleculeInfo[[MoleculeNo]]), "Fragment") %>% sum
  Transitions <- MoleculeInfo[[MoleculeNo]][["Transitions"]]
  NumberTransitions <- nrow(Transitions)
  
  CombinedProb <- matrix(NA, nrow = NumberTransitions, ncol = NumberTransitions)
  colnames(CombinedProb) <- rownames(Transitions)
  rownames(CombinedProb) <- rownames(Transitions)
  
  if (NumberFragments == 2) {
    for (x in seq_len(NumberTransitions)) {
      for (y in seq_len(NumberTransitions)) {
        CombinedProb[x, y] <- as.data.frame(CumProbList[[1]])[y, FragmentMassShift[[1]][x] + 1] * as.data.frame(CumProbList[[2]])[y, FragmentMassShift[[2]][x] + 1]
      }  #y
    }  #x
  } else if (NumberFragments == 1) {
    for (x in seq_len(NumberTransitions)) {
      for (y in seq_len(NumberTransitions)) {
        CombinedProb[x, y] <- as.data.frame(CumProbList[[1]])[y, FragmentMassShift[[1]][x] + 1]
      }  #y
    }  #x
  }
  
  return(CombinedProb)
  
}


# CALCULATION OF PROBABILITY MATRIX
#
# In this section, the probability matrix used for data correction is produced.
# To this end, for each molecule to be corrected, molecule/fragment
# information from the molecule information input file is used
# to calculate the probability that a certain transition or full MS ion x
# is derived from another transition/full MS ion y due to
# natural stable isotope abundance/tracer purity.

ProbNormalRes <- function(MoleculeInfo,MoleculesTotal,ElementInfo,CorrectTracerElementCore,
                          CorrectTracerImpurity,CalculationThreshold, verbose) {
  
  if(verbose){message(date(), " :: calculating probability matrix with parameters: ")
    message(date(), " :: :: CorrectTracerElementCore: ", CorrectTracerElementCore)
    message(date(), " :: :: CorrectTracerImpurity: ", CorrectTracerImpurity)
    message(date(), " :: :: CalculationThreshold: ", CalculationThreshold)}
  
  MaxMassShiftList <- list()
  icmResultList <- list()
  ecResultList <- list()
  msprobResultList <- list()
  mstResultList <- list()
  
  CombinedProbList <- list()
  
  for (MoleculeNo in seq_len(MoleculesTotal)) {
    
    MoleculeData <- MoleculeInfo[[MoleculeNo]]
    MoleculeName <- names(MoleculeInfo[MoleculeNo])
    
    NumberFragments <-
      sum(stringr::str_detect(names(MoleculeData), "Fragment"))
    Transitions <- MoleculeData[["Transitions"]]
    
    tmpIcmResultList <- list()
    tmpEcResultList <- list()
    tmpMsprobResultList <- list()
    tmpMstResultList <- list()
    tmpMaxMassShiftList <- list()
    
    # The following functions are run for both fragments of a molecule
    # in the case of MS/MS measurements or just for one molecule(-fragment) for
    # full-MS measurements.
    
    for (Fragment in seq_len(NumberFragments)) {
      # For each molecule, the expected transitions are checked and MaxMassShift
      # gives the mass shifts of the labelled molecule(-fragment(s)) with the highest mass.
      
      MaxMassShift <- apply(Transitions, 2, max)[Fragment]
      tmpMaxMassShiftList[[Fragment]] <- MaxMassShift
      
      # CALCULATE ISOTOPE COMBINATION PROBABILITIES AND MASS SHIFTS
      
      IsoCombinationsResult <-
        IsoCombinationsMaster(
          MoleculeData = MoleculeData,
          MoleculeName = MoleculeName,
          MoleculeNo = MoleculeNo,
          Fragment = Fragment,
          ElementInfo = ElementInfo,
          CorrectTracerElementCore = CorrectTracerElementCore,
          CorrectTracerImpurity = CorrectTracerImpurity,
          CalculationThreshold = CalculationThreshold
        )
      
      tmpIcmResultList[[Fragment]] <- IsoCombinationsResult
      
      # CALCULATE ELEMENT COMBINATION PROBABILITIES AND MASS SHIFTS
      
      ElementCombinationsResult <-
        ElementCombinations(
          MoleculeInfo = MoleculeInfo,
          MoleculeNo = MoleculeNo,
          Fragment = Fragment,
          IsoCombinationsResult = IsoCombinationsResult,
          CalculationThreshold = CalculationThreshold,
          CorrectTracerImpurity = CorrectTracerImpurity
        )
      
      tmpEcResultList[[Fragment]] <- ElementCombinationsResult
      
      #CALCULATE CUMULATIVE PROBABILITIES OF MASS SHIFTS
      
      CumProbList <-
        MassShiftProbabilities(
          MoleculeInfo = MoleculeInfo,
          MoleculeNo = MoleculeNo,
          Fragment = Fragment,
          ElementCombinations = ElementCombinationsResult,
          MaxMassShift = MaxMassShift,
          CorrectTracerImpurity = CorrectTracerImpurity
        )
      
      tmpMsprobResultList[[Fragment]] <- CumProbList
      
      #GET MASS SHIFTS OF MEASUREMENTS/TRANSITIONS
      
      mstResult <-
        MassShiftTransitions(
          MoleculeInfo = MoleculeInfo,
          MoleculeNo = MoleculeNo,
          Fragment = Fragment
        )
      
      tmpMstResultList[[Fragment]] <- mstResult
      
    }#Fragment
    names(tmpIcmResultList) <-
      paste0("Fragment_", seq(seq_len(NumberFragments)))
    names(tmpEcResultList) <-
      paste0("Fragment_", seq(seq_len(NumberFragments)))
    names(tmpMsprobResultList) <-
      paste0("Fragment_", seq(seq_len(NumberFragments)))
    names(tmpMstResultList) <-
      paste0("Fragment_", seq(seq_len(NumberFragments)))
    names(tmpMaxMassShiftList) <-
      paste0("Fragment_", seq(seq_len(NumberFragments)))
    
    # CALCULATE FINAL PROBABILTIY MATRIX CombFragmentProb
    
    CombFragmentProb <-
      FragmentsCombinedProbability(
        MoleculeInfo = MoleculeInfo,
        MoleculeNo = MoleculeNo,
        CumProbList = tmpMsprobResultList,
        FragmentMassShift =
          tmpMstResultList
      )
    
    CombinedProbList[[MoleculeNo]] <-
      CombFragmentProb
    
    MaxMassShiftList[[MoleculeNo]] <- tmpMaxMassShiftList
    icmResultList[[MoleculeNo]] <- tmpIcmResultList
    ecResultList[[MoleculeNo]] <- tmpEcResultList
    msprobResultList[[MoleculeNo]] <- tmpMsprobResultList
    mstResultList[[MoleculeNo]] <- tmpMstResultList
    
  }#MoleculeNo  
  
  names(MaxMassShiftList) <- names(MoleculeInfo)
  names(icmResultList) <- names(MoleculeInfo)
  names(ecResultList) <- names(MoleculeInfo)
  names(msprobResultList) <- names(MoleculeInfo) # CumProbList
  names(mstResultList) <- names(MoleculeInfo) # FragmentMassShift
  names(CombinedProbList) <-
    names(MoleculeInfo) # CombFragmentProb
  
  if(verbose){message(date(), " :: calculating probability matrix [OK]")}
  
  return(CombinedProbList)
  
}


# Data correction using the probability matrix

# The function 'RawDataCorrection' performs the actual
# data correction. It uses the values extracted from the measurement file
# and corrects molecule by molecule for each sample. This is done by
# numerically solving a linear equation system. Here,
# the uncorrected value of a given transition is a linear combination
# of the corrected transition values with their
# probability of contributing to the uncorrected value
# as their coefficients (derived from the probability matrix,
# ProbMatrixComplete). The algorithm works with the
# constraint that a solution of the linear equation
# system cannot be < 0.

#' @importFrom pracma lsqlincon
RawDataCorrection <- function(UncorrectedData, MoleculeData, MoleculeName,
                              ProbMatrix, MoleculeNo, SampleNo, SampleName,
                              roundDigits, logEnvironment, verbose) {
  NumberTransitions <- nrow(MoleculeData[["Transitions"]])
  
  ProbMatrixComplete <- ProbMatrix
  
  # In the following, the linear equation system is solved with linear
  # inequality constraints: ConstraintMatrix*SolutionVector <= ConstraintVector.
  # The ConstraintVector is 0, the ConstraintMatrix is -1 at each
  # ConstraintMatrix(TransitionNo, TransitionNo) position.
  # It is thereby assured that the solutions cannot be < 0.
  # If a solution is < 0, ConstraintMatrix*SolutionVector becomes positive and
  # the constraint ConstraintMatrix*SolutionVector <=
  # ConstraintVector (which is 0) is not fulfilled.
  
  ConstraintVector <- rep(0, NumberTransitions)
  ConstraintMatrix <- matrix(0, nrow = NumberTransitions, ncol = NumberTransitions)
  
  for (TransitionNo in seq_len(NumberTransitions)) {
    ConstraintMatrix[TransitionNo, TransitionNo] <- -1
  }
  
  # Check the vector of uncorrected values for values that are not a number
  # (e.g. NA) and thus missing.
  CheckForNaN <- is.na(UncorrectedData)
  
  # If there are no missing values, perform the correction for the
  # current vector of uncorrected values UncorrectedData.
  # This is done by numerically solving a linear equation system where each
  # uncorrected value is a linear combination of unknown corrected values with
  # their probabilites of contributing to the uncorrected value as the constants
  # (from ProbMatrix).  ||ProbMatrix*CorrTransitionsPlus - UncorrectedData||^2 is minimized in a linear least squares
  # approach with the constraint that the corrected values (CorrTransitionsPlus)
  # must not be < 0. In addition to CorrTransitionsPlus, also the residuals of
  # the solving process are given (CorrResiduals).
  
  if (sum(CheckForNaN) == 0) {
    CorrTransitionsPlus <- pracma::lsqlincon(C = ProbMatrix, d = as.numeric(UncorrectedData), A = ConstraintMatrix, b = ConstraintVector)
    ManResiduals <- ProbMatrix %*% CorrTransitionsPlus - UncorrectedData
  } else {
    notification <- paste0(
      "In measurement data file: Measurement data of molecule ", MoleculeName, " in sample ", SampleName, " contains NA values. The correction performed may be less acurate.",
      "\nBe especially careful when considering fraction and mean enrichment values from samples with missing values."
    )
    errorHandler(notification, logEnvironment, "warning", verbose = verbose)
    
    CorrTransitionsPlus <- vector()
    ManResiduals <- vector()
    
    # Get indices of missing values. The rows with these indices removed from
    # UncorrectedData and the ConstraintVector. Row and column with these
    # indices are removed from ProbMatrix and ConstraintMatrix.
    # This way, the missing values are completely removed from the
    # linear equation system.
    
    # get index of NaN values
    nan.index <- as.numeric(which(CheckForNaN))
    
    # remove respective rows/columns
    UncorrectedData.nan <- UncorrectedData[-nan.index]
    ProbMatrix.nan <- ProbMatrix[-nan.index, -nan.index]
    ConstraintMatrix.nan <- ConstraintMatrix[-nan.index, -nan.index]
    ConstraintVector.nan <- ConstraintVector[-nan.index]
    
    # Here, the linear equation system is solved with a set of equations that is reduced because of the missing values. If the missing values are expected to
    # be high, the effect on the accuracy of the correction is more pronounced than if they are expected to be close to 0.
    
    if (length(nan.index) < NumberTransitions) {
      CorrTransitionsPlusNaN <- pracma::lsqlincon(C = ProbMatrix.nan, d = as.numeric(UncorrectedData.nan), A = ConstraintMatrix.nan, b = ConstraintVector.nan)
      ManResidualsNaN <- ProbMatrix.nan %*% CorrTransitionsPlusNaN - UncorrectedData.nan
    } else {
      notification <- paste0(
        "In measurement data file: All measurements are NA for molecule ", MoleculeName, " in sample ", SampleName ,
        "."
      )
      errorHandler(notification, logEnvironment, "warning", verbose)
    }
    
    # To be able to correctly assign corrected values to their measurement tag when missing values are present, a full size corrected value vector
    # CorrTransitionsPlus is made from the reduced size corrected values vector CorrTransitionsPlusNaN.  This is done by adding the missing values as NA at
    # the associated vector index.  In the same way a full size CorrResiduals vector is produced. A warning is given in the WarningsLog for each missing
    # value.
    
    i <- 1
    
    # NumberTransitions<-ncol(ProbMatrix)
    
    for (TransitionNo in seq_len(NumberTransitions)) {
      if (!CheckForNaN[TransitionNo]) {
        # 2017-08-28
        
        CorrTransitionsPlus[TransitionNo] <- CorrTransitionsPlusNaN[i]
        ManResiduals[TransitionNo] <- ManResidualsNaN[i]
        
        i <- i + 1
      } else {
        CorrTransitionsPlus[TransitionNo] <- NA
        ManResiduals[TransitionNo] <- NA
      } # if
    } # TransitionNo
  } # sum(CheckForNaN)==0
  
  CorrTransitionsPlus <- round(CorrTransitionsPlus, digits = roundDigits)
  ManResiduals <- round(ManResiduals, digits = roundDigits)
  
  # The corrected values in CorrTransitionsPlus are values that correspond to
  # the full integral of the isotopologue abundance distribution of a given
  # labeled species. To get the value (CorrTransitions) that corresponds to the
  # species in the distribution that contains no natural isotopes of higher
  # mass, the corresponding CorrTransitionsPlus value has to be multiplied
  # with the probability that the species contains no isotopes of higher mass
  # due to natural abundance
  
  CorrTransitions <- vector()
  CorrTransitionsFractions <- vector()
  CorrTransitionsPlusFractions <- vector()
  RelativeResiduals <- vector()
  
  # Calculation of residuals relative to corrected data
  
  RelativeResiduals <- round(ManResiduals / CorrTransitionsPlus, digits = roundDigits)
  
  for (TransitionNo in seq_len(NumberTransitions)) {
    if (is.na(CorrTransitionsPlus[TransitionNo]) == FALSE) {
      CorrTransitions[TransitionNo] <- CorrTransitionsPlus[TransitionNo] * ProbMatrixComplete[TransitionNo, TransitionNo]
    } else {
      CorrTransitions[TransitionNo] <- NA
    }
  } # TransitionNo
  
  round(CorrTransitions, digits = roundDigits)
  
  # To be able to calculate fractions, the vectors of the corrected values are summed.
  CorrTransitionsPlusSum <- sum(CorrTransitionsPlus[which(CheckForNaN == FALSE)])
  CorrTransitionsSum <- sum(CorrTransitions[which(CheckForNaN == FALSE)])
  
  # Calculation of fractions
  
  CorrTransitionsFractions <- vector()
  CorrTransitionsPlusFractions <- vector()
  
  for (TransitionNo in seq_len(NumberTransitions)) {
    CorrTransitionsFractions[TransitionNo] <- CorrTransitions[TransitionNo] / CorrTransitionsSum
    CorrTransitionsPlusFractions[TransitionNo] <- CorrTransitionsPlus[TransitionNo] / CorrTransitionsPlusSum
  } # TransitionNo
  
  round(CorrTransitionsPlusFractions, digits = roundDigits)
  round(CorrTransitionsFractions, digits = roundDigits)
  
  names(CorrTransitions) <- colnames(ProbMatrix)
  names(CorrTransitionsPlus) <- colnames(ProbMatrix)
  names(ManResiduals) <- colnames(ProbMatrix)
  names(RelativeResiduals) <- colnames(ProbMatrix)
  names(CorrTransitionsFractions) <- colnames(ProbMatrix)
  names(CorrTransitionsPlusFractions) <- colnames(ProbMatrix)
  
  returnList <- list(
    CorrectedMonoisotopic = CorrTransitions, Corrected = CorrTransitionsPlus, CorrectedMonoisotopicFractions = CorrTransitionsFractions,
    CorrectedFractions = CorrTransitionsPlusFractions, Residuals = ManResiduals, RelativeResiduals = RelativeResiduals
  )
  
  return(returnList)
}


# DATA CORRECTION USING THE PROBABILITY MATRIX
#
# In this function, for each sample and molecule from the the input raw data the
# measured data will be corrected using the probability matrix generated before.

RawDataCorrectionOnCompleteDataset <- function(MeasurementInfo, MoleculeInfo, CombinedProbList, MoleculesTotal, SamplesTotal, NamesCorrectedOutput, roundDigits, logEnvironment, parallel, verbose) {
  if (verbose) {
    message(
      date(),
      " :: preparing ",
      MoleculesTotal,
      " molecules for correction [OK]"
    )
    
    message("\n", rep("#", options()$width))
    message(
      date(),
      " :: Starting data correction for <",
      SamplesTotal,
      " samples> ..."
    )
    message(rep("#", options()$width), "\n")
  }
  
  CorrectionResultList <- list()
  
  if (parallel){
    # Define a function to process each sample in parallel
    processSample <- function(SampleNo) {
      
      tmpCorrectionResult <- list()
      
      for (MoleculeNo in seq_len(MoleculesTotal)) {
        ProbMatrix <- CombinedProbList[[MoleculeNo]]
        
        MoleculeData <- MoleculeInfo[[MoleculeNo]]
        MoleculeName <- names(MoleculeInfo[MoleculeNo])
        
        TransitionLocation <- unlist(MoleculeData[["TransitionLocation"]])
        
        # For the molecule defined by MoleculeNo, get the vector of uncorrected values. The location of the value corresponding to a given transition is given by
        # TransitionLocation and based on finding measurement name tags in the input data (RawDataExtraction).
        
        UncorrectedData <- MeasurementInfo[TransitionLocation, SampleNo]
        
        # The function 'RawDataCorrection' performs the actual
        # data correction.
        CorrectionResult <- RawDataCorrection(
          UncorrectedData = UncorrectedData,
          MoleculeData = MoleculeData,
          MoleculeName = MoleculeName,
          ProbMatrix = ProbMatrix,
          MoleculeNo = MoleculeNo,
          SampleNo = SampleNo,
          SampleName = colnames(MeasurementInfo)[SampleNo],
          roundDigits = roundDigits,
          logEnvironment = logEnvironment, 
          verbose = verbose
        )
        
        tmpCorrectionResult[[MoleculeNo]] <- CorrectionResult
      } # MoleculeNo
      
      names(tmpCorrectionResult) <- names(MoleculeInfo)
      
      return(tmpCorrectionResult)
    }
    
    
    # Use mclapply to process samples in parallel
    CorrectionResultList <- parallel::mclapply(seq_len(SamplesTotal), processSample, mc.cores = (parallel::detectCores()-1))
  } else {
    for (SampleNo in seq_len(SamplesTotal)) {
      if (verbose) {
        message(
          date(),
          " :: processing Sample #",
          SampleNo,
          " [",
          colnames(MeasurementInfo)[SampleNo],
          "] "
        )
      }
      
      tmpCorrectionResult <- list()
      
      for (MoleculeNo in seq_len(MoleculesTotal)) {
        ProbMatrix <- CombinedProbList[[MoleculeNo]]
        
        MoleculeData <- MoleculeInfo[[MoleculeNo]]
        MoleculeName <- names(MoleculeInfo[MoleculeNo])
        
        TransitionLocation <- unlist(MoleculeData[["TransitionLocation"]])
        
        # For the molecule defined by MoleculeNo, get the vector of uncorrected values. The location of the value corresponding to a given transition is given by
        # TransitionLocation and based on finding measurement name tags in the input data (RawDataExtraction).
        
        UncorrectedData <- MeasurementInfo[TransitionLocation, SampleNo]
        
        # The function 'RawDataCorrection' performs the actual
        # data correction.
        CorrectionResult <- RawDataCorrection(
          UncorrectedData = UncorrectedData,
          MoleculeData = MoleculeData,
          MoleculeName = MoleculeName,
          ProbMatrix = ProbMatrix,
          MoleculeNo = MoleculeNo,
          SampleNo = SampleNo,
          SampleName = colnames(MeasurementInfo)[SampleNo],
          roundDigits = roundDigits,
          logEnvironment = logEnvironment, verbose = verbose
        )
        
        tmpCorrectionResult[[MoleculeNo]] <- CorrectionResult
      } # MoleculeNo
      
      names(tmpCorrectionResult) <- names(MoleculeInfo)
      
      CorrectionResultList[[as.character(SampleNo)]] <- tmpCorrectionResult
    } # SampleNo
  }
  
  names(CorrectionResultList) <- colnames(MeasurementInfo)
  
  x <- convert2array(input = CorrectionResultList,SamplesTotal = SamplesTotal)
  
  CorrectedDataOutput <- list()
  
  CorrectedDataOutput$RawData <- as.data.frame(MeasurementInfo)
  
  for (name in NamesCorrectedOutput$CorrectedDataNamesBase) {
    CorrectedDataOutput[[name]] <- as.data.frame(x[,,name])
  }
  
  return(list("CorrectedDataOutput" = CorrectedDataOutput, "CorrectionResultList" = CorrectionResultList))
}


# # no Parallel
# RawDataCorrectionOnCompleteDataset <- function(MeasurementInfo, MoleculeInfo, CombinedProbList, MoleculesTotal, SamplesTotal, NamesCorrectedOutput, roundDigits, logEnvironment, verbose) {
#   if (verbose) {
#     message(
#       date(),
#       " :: preparing ",
#       MoleculesTotal,
#       " molecules for correction [OK]"
#     )
#     
#     message("\n", rep("#", options()$width))
#     message(
#       date(),
#       " :: Starting data correction for <",
#       SamplesTotal,
#       " samples> ..."
#     )
#     message(rep("#", options()$width), "\n")
#   }
#   
#   CorrectionResultList <- list()
#   
#   for (SampleNo in seq_len(SamplesTotal)) {
#     if (verbose) {
#       message(
#         date(),
#         " :: processing Sample #",
#         SampleNo,
#         " [",
#         colnames(MeasurementInfo)[SampleNo],
#         "] "
#       )
#     }
#     
#     tmpCorrectionResult <- list()
#     
#     for (MoleculeNo in seq_len(MoleculesTotal)) {
#       ProbMatrix <- CombinedProbList[[MoleculeNo]]
#       
#       MoleculeData <- MoleculeInfo[[MoleculeNo]]
#       MoleculeName <- names(MoleculeInfo[MoleculeNo])
#       
#       TransitionLocation <- unlist(MoleculeData[["TransitionLocation"]])
#       
#       # For the molecule defined by MoleculeNo, get the vector of uncorrected values. The location of the value corresponding to a given transition is given by
#       # TransitionLocation and based on finding measurement name tags in the input data (RawDataExtraction).
#       
#       UncorrectedData <- MeasurementInfo[TransitionLocation, SampleNo]
#       
#       # The function 'RawDataCorrection' performs the actual
#       # data correction.
#       CorrectionResult <- RawDataCorrection(
#         UncorrectedData = UncorrectedData,
#         MoleculeData = MoleculeData,
#         MoleculeName = MoleculeName,
#         ProbMatrix = ProbMatrix,
#         MoleculeNo = MoleculeNo,
#         SampleNo = SampleNo,
#         SampleName = colnames(MeasurementInfo)[SampleNo],
#         roundDigits = roundDigits,
#         logEnvironment = logEnvironment, verbose = verbose
#       )
#       
#       tmpCorrectionResult[[MoleculeNo]] <- CorrectionResult
#     } # MoleculeNo
#     
#     names(tmpCorrectionResult) <- names(MoleculeInfo)
#     
#     CorrectionResultList[[SampleNo]] <- tmpCorrectionResult
#   } # SampleNo
#   
#   names(CorrectionResultList) <- colnames(MeasurementInfo)
# 
#   x <- convert2matrix(input = CorrectionResultList)
# 
#   CorrectedDataOutput <- list()
#   
#   CorrectedDataOutput$RawData <- as.data.frame(MeasurementInfo)
# 
#   for (name in NamesCorrectedOutput$CorrectedDataNamesBase) {
#     CorrectedDataOutput[[name]] <- as.data.frame(x[[name]])
#   }
#   
#   return(list("CorrectedDataOutput" = CorrectedDataOutput, "CorrectionResultList" = CorrectionResultList))
# }


# Calculates the mean isotopic enrichment for the corrected values.


CalcMeanEnrichment <- function(MoleculeInfo, MoleculesTotal, SamplesTotal, roundDigits, CorrectionResultList, UltraHighRes, verbose) {
  
  if(verbose){message("\n", date(), " :: calculating mean enrichment ...")}
  
  MeanEnrichmentResultList <- list()
  
  for (SampleNo in seq_len(SamplesTotal)) {
    
    MeanEnrichmentSample <- list()
    
    for (MoleculeNo in seq_len(MoleculesTotal)) {
      
      MoleculeData <- MoleculeInfo[[MoleculeNo]]
      MoleculeName <- names(MoleculeInfo[MoleculeNo])
      CorrectionData <- CorrectionResultList[[SampleNo]][[MoleculeNo]]
      Transitions <- MoleculeData[["Transitions"]]
      NumberTransitions <- nrow(Transitions)
      
      if (UltraHighRes == FALSE) {
        
        MeanEnrichmentTmp <- 0
        
        for (TransitionNo in seq_len(NumberTransitions)) {
          
          if (is.na(CorrectionData[["CorrectedFractions"]][TransitionNo]) == FALSE) {
            
            # It is ok to use the mass shift values and not neccessarily the label values here, because dividing by TotalMassShiftMax will in the end lead to correct
            # values
            
            TotalMassShift <- MoleculeData$Transitions$Precursor[[TransitionNo]]
            
            MeanEnrichmentTmp <- MeanEnrichmentTmp + as.numeric(CorrectionData[["CorrectedFractions"]][TransitionNo] * TotalMassShift)
            
          }
          
        }  #TransitionNo
        
        TotalMassShiftMax <- max(MoleculeData$TransitionsExpected$Precursor)
        
        MeanEnrichmentSample[[MoleculeName]] <- round(MeanEnrichmentTmp/TotalMassShiftMax, digits = roundDigits)
        
      } else if (UltraHighRes) {
        
        Tracers <- MoleculeData$Fragment_1$Tracer
        
        for (Tracer in names(Tracers)) {
          
          MoleculeNamePlusTracer <- paste0(MoleculeName, "_", Tracer)
          Label <- MoleculeData$Transitions[[Tracer]]
          MeanEnrichmentTmp <- 0
          
          for (TransitionNo in seq_len(NumberTransitions)) {
            
            if (is.na(CorrectionData[["CorrectedFractions"]][TransitionNo]) == FALSE) {
              
              MeanEnrichmentTmp <- MeanEnrichmentTmp + CorrectionData[["CorrectedFractions"]][TransitionNo] * Label[[TransitionNo]]
              
            }
            
          }
          
          MeanEnrichmentSample[[MoleculeNamePlusTracer]] <- round(MeanEnrichmentTmp/Tracers[[Tracer]], digits = roundDigits)
          
        }
        
      }
      
    }
    
    MeanEnrichmentResultList[[SampleNo]] <- MeanEnrichmentSample
    
  }
  
  names(MeanEnrichmentResultList) <- names(CorrectionResultList)
  
  MeanEnrichmentDataOutput <-
    as.data.frame(convertEnrichmentList2matrix(input = MeanEnrichmentResultList, verbose=verbose))
  
  if(verbose){message(date(), " :: calculating mean enrichment [OK]\n")}
  return(MeanEnrichmentDataOutput)
  
}  #MeanEnrichment


##### Final function ##########
NIAStripping <-
  function(MeasurementFile = NA,
           ElementFile = NA,
           MoleculeFile = NA,
           ObservationInfo = NULL,
           CorrectTracerImpurity = FALSE,
           CorrectTracerElementCore = TRUE,
           CalculateMeanEnrichment = TRUE,
           UltraHighRes = FALSE,
           DirOut = ".",
           FileOut = "result",
           FileOutFormat = "csv",
           ReturnResultsObject = TRUE,
           CorrectAlsoMonoisotopic = FALSE,
           CalculationThreshold = 10^-8,
           CalculationThreshold_UHR = 8,
           parallel = TRUE,
           verbose=FALSE,
           Testmode = FALSE) {
    message("\n", rep("#", options()$width))
    message(date(), " :: Starting IsoCorrection ...\n")
    # message(rep("#", options()$width))
    
    version <- "IsoCorrectoRBoost 0.1.0"
    
    timestamp <- Sys.time()
    timestampFMT <- strftime(timestamp, "%Y-%m-%d_%H%M%S")
    
    correctionOutput <- list()
    
    correctionOutput$success <- "FALSE"
    correctionOutput$results <- NA
    correctionOutput$log <- NA
    correctionOutput$error <- NA
    
    roundDigits <- 8
    
    # DEFINE LOG ENVIRONMENT
    
    log.env <- new.env()
    
    log.env$error <- character()
    log.env$warning <- character()
    log.env$general <- character()
    
    log.env$param <- list()
    
    log.env$param$Timestamp <- timestamp
    log.env$param$TimestampFMT <- timestampFMT
    
    log.env$param$version <- version
    log.env$param$date <- date()
    log.env$param$MeasurementFile <- MeasurementFile
    log.env$param$MoleculeFile <- MoleculeFile
    log.env$param$ElementFile <- ElementFile
    log.env$param$NamePrefix <- "IsoCorrectoR"
    log.env$param$WorkingDirectory <- getwd()
    log.env$param$LogfileName <- paste0(log.env$param$NamePrefix, ".log")
    log.env$param$FileOut <- FileOut
    log.env$param$FileOutFormat <- FileOutFormat
    
    log.env$param$CorrectTracerImpurity <- CorrectTracerImpurity
    log.env$param$CorrectTracerElementCore <- CorrectTracerElementCore
    log.env$param$UltraHighRes <- UltraHighRes
    log.env$param$CalculateMeanEnrichment <- CalculateMeanEnrichment
    
    if (!UltraHighRes) {
      log.env$param$threshold <- CalculationThreshold
    }
    else {
      log.env$param$threshold <- CalculationThreshold_UHR
    }
    
    tryCatch({
      
      # DEFINE OUTPUT DIRECTORY
      
      if (DirOut == "." || DirOut == "") {
        log.env$param$OutputDirectory <- log.env$param$WorkingDirectory
      }
      
      else {
        log.env$param$OutputDirectory <- DirOut
      }
      
      if (Testmode == FALSE) {
        log.env$param$OutputDirectory <-
          file.path(log.env$param$OutputDirectory, timestampFMT)
      }
      
      ifelse(
        !dir.exists(log.env$param$OutputDirectory),
        dir.create(log.env$param$OutputDirectory, recursive = TRUE),
        FALSE
      )
      
      # CHECK FUNCTION PARAMETERS
      
      checkIsoCorrectionParameters(
        MeasurementFile = MeasurementFile,
        ElementFile = ElementFile,
        MoleculeFile = MoleculeFile,
        CorrectTracerImpurity = CorrectTracerImpurity,
        CorrectTracerElementCore = CorrectTracerElementCore,
        CalculateMeanEnrichment = CalculateMeanEnrichment,
        UltraHighRes = UltraHighRes,
        FileOut = FileOut,
        FileOutFormat = FileOutFormat,
        ReturnResultsObject = ReturnResultsObject,
        CorrectAlsoMonoisotopic = CorrectAlsoMonoisotopic,
        CalculationThreshold = CalculationThreshold,
        CalculationThreshold_UHR = CalculationThreshold_UHR,
        logEnvironment = log.env,
        Testmode = Testmode
      )
      
      # DEFINE NAMES OF OUTPUT FILES
      
      NamesCorrectedOutput <- defineOutputFilenames(
        logEnvironment = log.env, FileOut = FileOut,
        FileOutFormat = FileOutFormat,
        CalculateMeanEnrichment = CalculateMeanEnrichment,
        CorrectAlsoMonoisotopic = CorrectAlsoMonoisotopic
      )
      
      # LOAD INPUT FILES, CHECK THEIR STRUCTURE, EXTRACT INFORMATION AND CHECK 
      # INPUT LOGIC
      
      InputFileInformation <- ExtractInputFileInformation(
        ElementFile = ElementFile, MoleculeFile = MoleculeFile,
        MeasurementFile = MeasurementFile, UltraHighRes = UltraHighRes,
        CorrectTracerImpurity = CorrectTracerImpurity,
        logEnvironment = log.env,
        verbose = verbose
      )
      
      MoleculeInfo <- InputFileInformation$MoleculeInfo
      ElementInfo <- InputFileInformation$ElementInfo
      MeasurementInfo <- InputFileInformation$MeasurementInfo
      
      MoleculesTotal <- length(MoleculeInfo)
      SamplesTotal <- ncol(MeasurementInfo)
      
      # CALCULATION OF PROBABILITY MATRIX
      
      # In this section, the probability matrix used for data correction 
      # is produced.
      
      if (verbose) {
        message(
          date(),
          " :: preparing ",
          MoleculesTotal,
          " molecules for correction ..."
        )
      }
      
      # if(verbose){message(date(), " :: UltraHighRes==", UltraHighRes)}
      if (verbose) {
        if (UltraHighRes == FALSE) {
          message(date(), " :: running in normal resolution mode.")
        } else {
          message(date(), " :: running in ultra high resolution mode,")
        }
      }
      
      
      # NORMAL RESOLUTION
      
      if (!UltraHighRes) {
        CombinedProbList <- ProbNormalRes(
          MoleculeInfo = MoleculeInfo,
          MoleculesTotal = MoleculesTotal,
          ElementInfo = ElementInfo,
          CorrectTracerElementCore = CorrectTracerElementCore,
          CorrectTracerImpurity = CorrectTracerImpurity,
          CalculationThreshold = CalculationThreshold,
          verbose = verbose
        )
        
        # HIGH RESOLUTION
      } else if (UltraHighRes) {
        CombinedProbList <- ProbUltraHighRes(
          MoleculeInfo = MoleculeInfo,
          MoleculesTotal = MoleculesTotal,
          ElementInfo = ElementInfo,
          CorrectTracerElementCore = CorrectTracerElementCore,
          CorrectTracerImpurity = CorrectTracerImpurity,
          CalculationThreshold_UHR = CalculationThreshold_UHR,
          verbose = verbose
        )
      } # (!UltraHighRes)
      
      # DATA CORRECTION USING THE PROBABILITY MATRIX
      
      CorrectedData <- RawDataCorrectionOnCompleteDataset(
        MeasurementInfo = MeasurementInfo,
        MoleculeInfo = MoleculeInfo,
        CombinedProbList = CombinedProbList,
        MoleculesTotal = MoleculesTotal,
        SamplesTotal = SamplesTotal,
        NamesCorrectedOutput = NamesCorrectedOutput,
        roundDigits = roundDigits,
        logEnvironment = log.env,
        parallel = parallel,
        verbose = verbose
      )

      CorrectedDataOutput <- CorrectedData$CorrectedDataOutput
      CorrectionResultList <- CorrectedData$CorrectionResultList

      # CALCULATION OF MEAN ENRICHMENT

      if (CalculateMeanEnrichment) {
        MeanEnrichmentDataOutput <-
          CalcMeanEnrichment(
            MoleculeInfo = MoleculeInfo,
            MoleculesTotal = MoleculesTotal,
            SamplesTotal = SamplesTotal,
            roundDigits = roundDigits,
            CorrectionResultList = CorrectionResultList,
            UltraHighRes = UltraHighRes,
            verbose = verbose
          )

        CorrectedDataOutput$MeanEnrichment <- MeanEnrichmentDataOutput
      }

      # TRANSPOSE BACK AND ADD OBSERVATION INFORMATION
      if (!is.null(ObservationInfo)){
        CorrectedDataOutput <- TransposeAddObservationInfo(CorrectedDataOutput, NamesCorrectedOutput, ObservationInfo)
      }

      # WRITE CORRECTED DATA TO FILE

      writeCorrectionResultsToFile(
        CorrectedDataOutput = CorrectedDataOutput,
        NamesCorrectedOutput = NamesCorrectedOutput,
        logEnvironment = log.env,
        verbose = verbose
      )

      # END OF CORRECTION

      if (ReturnResultsObject) {
        correctionOutput$results <- CorrectedDataOutput
      }

      if(length(log.env$warning) == 0) {

        correctionOutput$success <- "TRUE"

      } else {

        correctionOutput$success <- "WARNINGS"

      }

      # verbose independent message
      # message("\n\n", rep("#", options()$width))
      message(date(), " :: IsoCorrection successfully completed!")
      message(rep("#", options()$width), "\n")
    },

    # ERROR HANDLING
    error = function(e) {
      commonErrorString <-
        "\n\nTHE CORRECTION PROCEDURE WAS ABORTED BECAUSE AN ERROR HAS OCCURED.
      PLEASE CHECK YOUR LOG-FILE.\n\n"

      log.env$error <- c(log.env$error, e)

      # verbose independent message

      message(paste0(commonErrorString, e))

    }
    )

    tryCatch({
      writeLog(log.env)
    },
    error = function(e) {
      log.env$error <- c(log.env$error, e)

      # verbose independent message

      message(paste0("THE LOG-FILE COULD NOT BE WRITTEN. ERROR: ", e[[1]]))

    }
    )

    correctionOutput$log <- as.list(log.env)

    correctionOutput$error <- log.env$error[1][[1]]

    return(correctionOutput)
  } # IsoCorrection