# functions to process and load data easily


loadBRCADNAm <- function(signalMat=TRUE, signalCoord=TRUE, 
                         loadingMat=TRUE, pcScores=TRUE,
                         patientMetadata=TRUE,
                         .env=currentEnv) {
    # making sure parent.frame is evaluated inside function (not outside as 
    # when listed as default argument)
    currentEnv = parent.frame(n=1) # env where function was called
    .env = currentEnv
    
    if (signalMat || signalCoord) {
        simpleCache("combinedBRCAMethyl_noXY", assignToVariable = "brcaMList")
    }
    if (signalMat) {
        #restrict patients included in this analysis
        
        
        brcaMetadata = fread(paste0(Sys.getenv("CODE"), "COCOA_paper/metadata/brca_metadata.csv"))
        # only keep patients who have definitive status for ER and PGR
        brcaMetadata = brcaMetadata[brcaMetadata$ER_status %in% 
                                        c("Positive", "Negative"), ]
        brcaMetadata = brcaMetadata[brcaMetadata$PGR_status %in% 
                                        c("Positive", "Negative"), ]
        brcaMetadata = brcaMetadata[brcaMetadata$subject_ID %in% 
                                        colnames(brcaMList[["methylProp"]]), ]
        
        
        
        # brcaMetadata should have already screened out patients without ER/PGR status
        # resulting in 657 patients
        hasER_PGR_IDs = as.character(brcaMetadata[, subject_ID])
        filteredMData = brcaMList[["methylProp"]][, hasER_PGR_IDs] 
        
        assign("signalMat", filteredMData, envir=.env)
    }
    
    
    
    if (signalCoord) {
        assign("signalCoord", brcaMList$coordinates, envir=.env)
    }
    
    if (loadingMat || pcScores) {
        simpleCache("allMPCA_657", assignToVariable = "allMPCA")
    }
    if (loadingMat) {
        assign("loadingMat", allMPCA$rotation, envir=.env)
    }
    if (pcScores) {
        assign("pcScores", allMPCA$x, envir=.env)
    }
    
    #####
    if (patientMetadata) {
        brcaMetadata = fread(paste0(Sys.getenv("CODE"), 
                                    "COCOA_paper/metadata/brca_metadata.csv"))
        # only keep patients who have definitive status for ER and PGR
        brcaMetadata = brcaMetadata[brcaMetadata$ER_status %in% 
                                        c("Positive", "Negative"), ]
        brcaMetadata = brcaMetadata[brcaMetadata$PGR_status %in% 
                                        c("Positive", "Negative"), ]
        
        
        # indexed clinical data (brca_clinical_metadata.tsv) has more up to date follow up info
        brcaMetadata2 = read.table(file = paste0(Sys.getenv("CODE"), 
                                                 "COCOA_paper/metadata/brca_clinical_metadata.tsv"), 
                                   sep = "\t", header = TRUE)
        
        pMetadata = merge(brcaMetadata, brcaMetadata2[, c("bcr_patient_barcode", 
                                                          "vital_status", 
                                                          "days_to_death", 
                                                          "days_to_last_follow_up")],
                          by.x="subject_ID", 
                          by.y="bcr_patient_barcode", all.x=TRUE)
        
        row.names(pMetadata) <- pMetadata$subject_ID
        assign("patientMetadata", pMetadata, envir=.env)
    }
    
    message(paste0(paste(c("signalMat", "signalCoord", 
                           "loadingMat", "pcScores",
                           "patientMetadata")[c(signalMat, signalCoord,
                                                loadingMat, pcScores,
                                                patientMetadata)], 
                         collapse =" "), 
                   " loaded into the environment."))
    
}

loadBRCAatac <- function(signalMat=TRUE, signalCoord=TRUE, .env=currentEnv){
    
    # making sure parent.frame is evaluated inside function (not outside as 
    # when listed as default argument)
    currentEnv = parent.frame(n=1) # env where function was called
    .env = currentEnv
    
    ryan_brca_count_matrix <- fread(ffData("tcga/ATACseq/TCGA-ATAC_BRCA_peaks_counts.tsv"), 
                                    header=TRUE)
    
    if (signalMat) {
        
        counts   <- as.matrix(ryan_brca_count_matrix[,2:75])
        counts   <- counts[, order(colnames(counts))]
        # $sample is actually region ID (e.g. "BRCA_100")
        row.names(counts) <- ryan_brca_count_matrix$sample
        # tcount   <- t(counts)
        # colnames(tcount) <- ryan_brca_count_matrix$sample
        
        # # Add metadata for each sample that has it
        # metadata <- fread(ffProc("COCOA_paper/analysis/atac/scores/brca/tcga_brca_metadata.csv"))
        # metadata <- metadata[!duplicated(metadata$subject_ID),]
        # metadata <- metadata[order(subject_ID),]
        
        # merged    <- as.data.frame(tcount)
        # merged$id <- rownames(tcount)
        # merged    <- merge(merged, metadata, by.x="id", by.y="subject_ID")
        # row.names(merged) <- merged$id
        # merged = merged[, 2:215921]
        assign("signalMat", counts, envir = .env)
    }
    
    if (signalCoord) {
        peaks           <- ryan_brca_count_matrix[, .(Chromosome, Start, End)]
        colnames(peaks) <- c("chr","start","end")
        pGR             <- makeGRangesFromDataFrame(peaks)
        assign("signalCoord", pGR, envir = .env)
    }
    
    
    message(paste0(paste(c("signalMat", "signalCoord")[c(signalMat, signalCoord)], 
                         collapse =" "), 
                   " loaded into the environment."))
}

# gets coordinates for CLL methyl
# returns a list with methylProp that has methylation and "methylCoord"
# that has corresponding genomic coordinates
# filters out X and Y chromosomes
prepareCLLMethyl = function(removeXY=TRUE) {
    
    # get microarray data
    eh = ExperimentHub()
    meth = eh[[names(query(eh, "CLLmethylation"))]] # EH1071
    # rows are cpgs, columns are samples
    methData = assay(meth)
    dataProbeNames = row.names(methData)
    
    
    # match probe names with genomic coordinates
    # ls('package:FDb.InfiniumMethylation.hg19')
    m450kAnno = get450k()
    length(m450kAnno)
    
    # get coordinates in same order as CLL data
    methCoord = m450kAnno[dataProbeNames]
    all(names(methCoord) == dataProbeNames)
    methCoordDT = COCOA:::grToDt(methCoord)
    # keep start coordinate as CpG site
    methCoordDT = methCoordDT[, .(chr, start)]
    
    if (removeXY) {
        xyInd = methCoordDT$chr %in% c("chrX", "chrY")
        methCoordDT = methCoordDT[!xyInd, ]
        methData = methData[!xyInd, ]
    }
    
    methCoord = COCOA:::dtToGr(methCoordDT)
    if (nrow(methCoordDT) != nrow(methData)) {
        stop("error matching probes to coordinates")
    }
    
    return(list(methylProp = methData, methylCoord = methCoord))
    
}


loadMOFAData <- function(methylMat=TRUE, signalCoord=TRUE, latentFactorMat=TRUE, 
                         factorWeights=FALSE, cllMultiOmics=FALSE, .env=currentEnv) {
    
    # making sure parent.frame is evaluated inside function (not outside as 
    # when listed as default argument)
    currentEnv = parent.frame(n=1) # env where function was called
    .env = currentEnv
    
    library(ExperimentHub)
    library("SummarizedExperiment")
    # library(MOFAtools)
    library(MOFA)
    library("MOFAdata")
    library("FDb.InfiniumMethylation.hg19")
    library(MultiAssayExperiment)
    
    if (signalCoord) {
        cllMethyl = prepareCLLMethyl()
        methCoord = cllMethyl$methylCoord
        assign(x = "signalCoord", value = methCoord, envir = .env)
        
    }
    
    # if both are loaded, put samples in same order
    if (methylMat || latentFactorMat) {
        
        # load but don't yet assign
        if (methylMat) {
            cllMethyl = prepareCLLMethyl()
            methData = cllMethyl$methylProp
        }
        
        # load but don't yet assign
        if (latentFactorMat) {
            # Loading an existing trained model
            filepath <- system.file("extdata", "CLL_model.hdf5",
                                    package = "MOFAdata")
            
            MOFAobject <- loadModel(filepath)
            
            LFs <- getFactors(
                MOFAobject,
                as.data.frame = FALSE
            )
        }
        if (methylMat && latentFactorMat) {
            sharedNames = row.names(LFs)[row.names(LFs) %in% colnames(methData)]
            
            LFs = LFs[sharedNames, ]
            methData = methData[, sharedNames]
            
            assign("latentFactorMat", LFs, envir = .env)
            assign("methylMat", methData, envir = .env)
        } else if (methylMat) {
            assign("methylMat", methData, envir = .env)
        } else if (latentFactorMat) {
            assign("latentFactorMat", LFs, envir = .env)
        }
    }
    
    if (factorWeights) {
        # Loading an existing trained model
        filepath <- system.file("extdata", "CLL_model.hdf5",
                                package = "MOFAdata")
        
        MOFAobject <- loadModel(filepath)
        
        MOFAweights <- getWeights(
            MOFAobject, 
            views = "Methylation", 
            factors = "all", 
            as.data.frame = TRUE
        )
        assign("factorWeights", MOFAweights, envir = .env)
    }
    
    if (cllMultiOmics) {
        
        data("CLL_data", package="MOFAdata") 
        
        # only reorder if both are being loaded
        if (cllMultiOmics && latentFactorMat) {
            # Loading an existing trained model
            filepath <- system.file("extdata", "CLL_model.hdf5",
                                    package = "MOFAdata")
            
            MOFAobject <- loadModel(filepath)
            
            LFs <- getFactors(
                MOFAobject,
                as.data.frame = FALSE
            )
            
            # put multiOmicsData in same order as latentFactors
            for (i in seq_along(CLL_data)) {
                
                sharedNames = row.names(LFs)[row.names(LFs) %in% colnames(CLL_data[[i]])]
                CLL_data[[i]] <- CLL_data[[i]][, sharedNames]
            }
        }
        
        assign("cllMultiOmics", CLL_data, envir = .env)
    }
    
    message(paste0(paste(c("methylMat", "signalCoord", 
                           "latentFactorMat", "factorWeights", "cllMultiOmics")[c(methylMat, signalCoord,
                                                                                  latentFactorMat, factorWeights, cllMultiOmics)], 
                         collapse =" "), 
                   " loaded into the environment."))
    
}

# loads GRList, rsName, rsDescription

loadGRList <- function(genomeV = "hg38", .env=currentEnv) {
    # making sure parent.frame is evaluated inside function (not outside as 
    # when listed as default argument)
    currentEnv = parent.frame(n=1) # env where function was called
    .env = currentEnv
    
    if (genomeV == "hg38") {
        source(paste0(Sys.getenv("CODE"), "COCOA_paper/src/load_process_regions_brca.R"))
        names(GRList) = rsName
        assign(x = "GRList", GRList, envir = .env)
        assign(x = "rsName", rsName, envir = .env)
        assign(x = "rsDescription", rsDescription, envir = .env)
    } else if (genomeV == "hg19") {
        source(paste0(Sys.getenv("CODE"), "COCOA_paper/src/load_process_regionDB_hg19.R"))
        names(GRList) = rsName
        assign(x = "GRList", GRList, envir = .env)
        assign(x = "rsName", rsName, envir = .env)
        assign(x = "rsDescription", rsDescription, envir = .env)
    } else {
        stop("Only hg38 and hg19 are available in this function currently.")
    }
    
}

# function for processing DNA methylation data from curatedTCGAdata
# also could work for other methylation arrays
# filter out sex chromosomes
# for methylMat and covMat (if provided), rows are CpGs and columns are samples 
# signalCoord is DF
filtMethylMat = function(signalCoord, methylMat, covMat = NULL, removeXY=TRUE, removeCpGWithNA=TRUE, genomeV = "hg19") {
    
    
    library(FDb.InfiniumMethylation.hg19)
    methylList = list()
    
    if (genomeV != "hg19") {
        stop("only hg19 is supported.")
    }
    
    if (is.null(signalCoord)) {
        dataProbeNames = row.names(methylMat)
        # match probe names with genomic coordinates
        # ls('package:FDb.InfiniumMethylation.hg19')
        m450kAnno = get450k()
        length(m450kAnno)
        # get coordinates in same order as CLL data
        methCoord = m450kAnno[dataProbeNames]
        all(names(methCoord) == dataProbeNames)
        methCoordDT = COCOA:::grToDt(methCoord)
        # keep start coordinate as CpG site
        methCoordDT = methCoordDT[, .(chr, start)]
        signalCoord = as.data.frame(methCoordDT)
    }
    
    
    if (removeXY) {
        xyInd = signalCoord$chr %in% c("chrX", "chrY")
        signalCoord = signalCoord[!xyInd, ]
        methylMat = methylMat[!xyInd, ]
    }
    
    if (removeCpGWithNA) {
        cpgHasNA = apply(methylMat, 1, function(x) any(is.na(x)))
        methylMat = methylMat[!cpgHasNA, ]
        signalCoord = signalCoord[!cpgHasNA, ]
        
    }
    
    if (nrow(signalCoord) != nrow(methylMat)) {
        stop("error matching probes to coordinates")
    }
    
    methylList[["coordinates"]] = signalCoord
    methylList[["methylProp"]] = methylMat
    
    
    return(methylList)
}

# load list that has coordinates and methylation matrix for a single cancer
# 450k microarrays only
# assigns methylList, pMeta to environment
loadTCGAMethylation <- function(cancerID, methylList=TRUE, pMeta=TRUE,
                                removeXY=TRUE, removeCpGWithNA=TRUE, 
                                genomeV="hg19", .env=currentEnv) {
    # making sure parent.frame is evaluated inside function (not outside as 
    # when listed as default argument)
    currentEnv = parent.frame(n=1) # env where function was called
    .env = currentEnv
    
    library(curatedTCGAData)
    library(TCGAutils)

    mData = curatedTCGAData(diseaseCode = cancerID, assays = c("*methyl450*"), dry.run = FALSE)

    if (methylList) {
        testM = assays(mData)[[1]]
        testM = as.matrix(testM)
        ##### get CpG coordinates
        # match probe names to coordinates
        # screen out CpGs that have any NAs and the XY chromosomes
        mList = filtMethylMat(signalCoord=NULL, methylMat=testM, covMat = NULL, 
                              removeXY=removeXY, removeCpGWithNA=removeCpGWithNA, 
                              genomeV = genomeV)
        
        assign("methylList", mList, envir = .env)
    }
    
    if (pMeta) {
        # get patient metadata, stored in colData(curatedTCGAData())
        metaDataCols = getClinicalNames(cancerID)
        allMeta = colData(mData) 
        tcgaMetadata = allMeta[, metaDataCols]
        assign("pMeta", tcgaMetadata, envir = .env)
    }
    
}
