rm(list = ls())
workDir = "/Users/sindiris/R Scribble/RNASeq.Fusion.data"
setwd(paste0(workDir,"./AnnotatedFusionCalls"))

## REad well annotated DataFrame
metaData <- read.csv(paste0(workDir,"/MetadataMapper.v3.txt"), header = T, sep = "\t", stringsAsFactors = FALSE); 
# metaData <- metaData %<>% dplyr::rename(SampleName = Biowulf_SampleName)
head(metaData)

cols <- c("Sample.ID", "case_id", "patient_id", "left_gene", "right_gene", "left_chr", "left_position", "right_chr", "right_position",
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
  fileName <- gsub(".annotated_fusions.txt", "",x ); samplenameList <- unlist(strsplit(fileName,"_")); 
  samplename <- paste0(samplenameList[3:length(samplenameList)], collapse="_")
  data <- cbind(Sample.Biowulf.ID =  samplename, data )
  colnames(data) <- cols
  return(data)
}

### List all the files ####
fusionFiles <- list.files("./fusions.annotated/"); length(fusionFiles)

### concatenate files
finalResultMatrix <- rbindlist(lapply(fusionFiles, makResultMatrix), fill = TRUE)[,-c(32,33)]
dim(finalResultMatrix); View(head(finalResultMatrix))

### Append metadata to the finalResultMatrix ###
finalResultMatrixJoin <- dplyr::left_join(finalResultMatrix, metaData[,
                                                 c("Patient.ID", "Sample.ID", "Sample.ID.Alias", "DIAGNOSIS.Alias", "LIBRARY_TYPE", "Case.ID")], 
                                                 by="Sample.ID")
head(finalResultMatrixJoin)

## Sanity Check1 : To evaluate the join // Check if join produced any NA
df <- finalResultMatrixJoin %>% data.frame()
df$Any_NA <- apply(df[,grep(paste( c("Patient.ID", "Sample.ID", "Sample.ID.Alias", "DIAGNOSIS.Alias", "LIBRARY_TYPE"),collapse = "|"), names(df) )],1, function(x) anyNA(x))
df_NA_TRUE <- df[which(df$Any_NA == "TRUE"),] %>% dplyr::distinct(Sample.ID,DIAGNOSIS.Alias,LIBRARY_TYPE);dim(df_NA_TRUE); View(df_NA_TRUE)

## Count samples with Defuse
finalResultMatrixJoinDefuse <- finalResultMatrixJoin %>% filter( grepl('Defuse', tool)) %>%
                                          dplyr::distinct(Patient.ID, Sample.ID, Case.ID)
dim(finalResultMatrixJoinDefuse); View(finalResultMatrixJoinDefuse)
write.table(finalResultMatrixJoinDefuse, "../finalResultMatrixJoinDefuse.txt", sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
## Keep only max spanning read count.