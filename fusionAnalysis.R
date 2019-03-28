rm(list = ls())
workDir = "/Users/sindiris/R Scribble/RNASeq.Fusion.data"
setwd(paste0(workDir,"./AnnotatedFusionCalls"))

## REad well annotated DataFrame
metaData <- read.csv(paste0(workDir,"/MetadataMapper.v3.txt"), header = T, sep = "\t", stringsAsFactors = FALSE); 
# metaData <- metaData %<>% dplyr::rename(SampleName = Biowulf_SampleName)
head(metaData)

cols <- c("case_id", "patient_id", "left_gene", "right_gene", "left_chr", "left_position", "right_chr", "right_position",
          "sample_id", "tool", "spanreadcount", "type", "var_level","ACMG_v2" ,"Brain_Cancer","CGCensus_Hereditary","Colon_Cancer",
          "Fusion", "Fusion_Genes","Fusion_Genes_Illumina","Hotspot_Genes","Inherited_Diseases","Kidney_Cancer","Loss_Of_Function_Genes",
          "Melanoma_Germline","NIAID_Panel","Neuroblastoma_Germline","PAX3-FOXO1_Targets_ChipSeq","Tumor_Suppressor_Genes","clinomics_gene_list")

## Empty data matrix ####
emptyDF <- data.frame(matrix(ncol = 32)); colnames(emptyDF) <- cols

## function to concatenate files to a single matrix 
makResultMatrix <- function(x) {
  #print(x)
  data <- tryCatch({
    
    read.csv(paste0("./fusions.annotated/",x),sep="\t", header = TRUE, stringsAsFactors = FALSE )
    
  }, error = function(e){
    
    emptyDF
  })
  #print(ncol(data))
  # fileName <- gsub(".annotated_fusions.txt", "",x ); samplenameList <- unlist(strsplit(fileName,"_")); 
  # samplename <- paste0(samplenameList[3:length(samplenameList)], collapse="_")
  # data <- cbind(Sample.Biowulf.ID =  samplename, data )
  # colnames(data) <- cols
  return(data)
}

### List all the files ####
fusionFiles <- list.files("./fusions.annotated/"); length(fusionFiles)

### concatenate files and change "sample_id" to "Sample.ID"
finalResultMatrix <- rbindlist(lapply(fusionFiles, makResultMatrix), fill = TRUE)[,-c(32,33)]
dim(finalResultMatrix); View(head(finalResultMatrix))
finalResultMatrix <- finalResultMatrix %>% dplyr::rename(Sample.ID=sample_id)

### Append metadata to the finalResultMatrix ###
finalResultMatrixJoin <- dplyr::left_join(finalResultMatrix, metaData[,
                                                 c("Patient.ID", "Sample.ID", "Sample.ID.Alias", "DIAGNOSIS.Alias", "LIBRARY_TYPE", "Case.ID")], 
                                                 by="Sample.ID")
head(finalResultMatrixJoin)

## Sanity Check1 : To evaluate the join // Check if join produced any NA
df <- finalResultMatrixJoin %>% data.frame()
df$Any_NA <- apply(df[,grep(paste( c("Patient.ID", "Sample.ID", "Sample.ID.Alias", "DIAGNOSIS.Alias", "LIBRARY_TYPE"),collapse = "|"), names(df) )],1, function(x) anyNA(x))
df_NA_TRUE <- df[which(df$Any_NA == "TRUE"),] %>% dplyr::distinct(Sample.ID,DIAGNOSIS.Alias,LIBRARY_TYPE);dim(df_NA_TRUE); View(df_NA_TRUE)

## Get only complecases 
finalResultMatrixJoin <- finalResultMatrixJoin %>% filter(!Sample.ID %in% c(df_NA_TRUE$Sample.ID))
dim(finalResultMatrixJoin) ; length(unique(finalResultMatrixJoin$Sample.ID))

## Make a key 
finalResultMatrixJoin$key <- paste(finalResultMatrixJoin$left_chr, finalResultMatrixJoin$left_gene, finalResultMatrixJoin$right_chr, 
                                   finalResultMatrixJoin$right_gene, sep="_")
## Count samples with Defuse
finalfusionResultMatrix.Defuse <- finalResultMatrixJoin %>% filter( grepl('Defuse', tool)) %>%
                                          dplyr::distinct(Patient.ID, Sample.ID, Case.ID)
dim(finalfusionResultMatrix.Defuse); View(finalfusionResultMatrix.Defuse)
write.table(finalResultMatrixJoin, "../FinalfusionResultMatrix.txt", sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
## Keep only max spanning read count.

##Read Annotation Data
Annotation  <- readRDS("C:/Users/sindiris/R Scribble/Annotation RDS/annotation_ENSEMBL_gene.RDS") ; Annotation$Start <- as.numeric(as.character(Annotation$Start)) ; Annotation$End <- as.numeric(as.character(Annotation$End))

##Filter Bogus Genes
geneFusionRight <-  as.character( unique(finalResultMatrixJoin$right_gene) )
rightGenePair   <- geneFusionRight[which(!geneFusionRight %in% Annotation$GeneName)]
geneFusionLeft  <- as.character(unique(finalResultMatrixJoin$left_gene) )
leftGenePair    <- geneFusionLeft[which(!geneFusionLeft %in% Annotation$GeneName)]

## Remove disgusting samples . "NCIEWS5000muscle_T_D23N9ACXX" , "NCI0228normal_T_C3JV2ACXX" . I dont understand why tumor patient normal is used in
## RNASeq cohort experiment
finalResultMatrixJoin.NoPatientNormal <-  finalResultMatrixJoin %>% filter()
dim(finalResultMatrixJoin.NoPatientNormal)

##Filter Fusion File
dim(finalResultMatrixJoin)
fusionFile.NS  <- finalResultMatrixJoin.NoPatientNormal %>% filter( DIAGNOSIS.Alias %in% c("NS") ) ; dim(fusionFile.NS )
finalResultMatrixJoin.NoNS <- finalResultMatrixJoin.NoPatientNormal %>% filter(! key %in% fusionFile.NS$key ); dim(finalResultMatrixJoin.NoNS)
## Remove bogus genes
fusionFile  <- finalResultMatrixJoin.NoNS %>% filter( !(left_gene %in% leftGenePair) ) %>%  filter( !(right_gene %in% rightGenePair) ) ; dim(fusionFile)
## Keep the selected tier
fusionFile <- fusionFile %>% filter(var_level %in% c("Tier 1.1", "Tier 1.2", "Tier 1.3", "Tier 2.1", "Tier 2.2")) ; dim(fusionFile)

### Extra Filtering and Condensing
fusionFileSplits <- fusionFile %>% tidyr::separate(spanreadcount, c("SPTool1", "SPTool2", "SPTool3"), " ")
fusionFileSplits[is.na(fusionFileSplits)] <- 0
head(fusionFileSplits)

fusionFileFilt.v1 <- fusionFileSplits %>% filter( !(tool %in% c("STAR-fusion", "FusionCatcher")) & SPTool1 >= 5 | SPTool2 >= 5 | SPTool3 >= 5)
fusionFileFilt.v2 <- fusionFileFilt.v1 %>% filter( !(var_level %in% c("Tier 2.1")) & SPTool1 >= 10 | SPTool2 >= 10 | SPTool3 >= 10) ; dim(fusionFileFilt.v2)
write.table(fusionFileFilt.v2 , "../FinalFilteredfusionResultMatrix.txt", sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)