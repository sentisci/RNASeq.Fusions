rm(list = ls())
workDir = "/Users/sindiris/R Scribble/RNASeq.Fusion.data"
setwd(paste0(workDir,"./AnnotatedFusionCalls"))

## REad well annotated DataFrame
metaData <- read.csv(paste0(workDir,"/MetadataMapper.v2.txt"), header = T, sep = "\t", stringsAsFactors = FALSE); 
# metaData <- metaData %<>% dplyr::rename(SampleName = Biowulf_SampleName)
head(metaData)

## Empty data matrix ####
emptyDF <- data.frame( 
  "case_id"                       = "",
  "patient_id"                    = "",
  "left_gene"                     = "",
  "right_gene"                    = "",  
  "left_chr"                      = "", 
  "left_position"                 = "",
  "right_chr"                     = "",
  "right_position"                = "", 
  "sample_id"                     = "",
  "tool"                          = "",   
  "spanreadcount"                 = "",
  "type"                          = "",
  "var_level"                     = "",
  "ACMG_v1"                       = "", 
  "ACMG_v2"                       = "",
  "Brain_Cancer"                  = "",  
  "CGCensus_Hereditary"           = "",
  "Colon_Cancer"                  = "",
  "Fusion"                        = "",
  "Fusion_Genes"                  = "",   
  "Fusion_Genes_Illumina"         = "",
  "Hotspot_Genes"                 = "", 
  "Inherited_Diseases"            = "",
  "Kidney_Cancer"                 = "",
  "Loss_Of_Function_Genes"        = "",   
  "Melanoma_Germline"             = "",
  "NIAID_Panel"                   = "",
  "Neuroblastoma_Germline"        = "", 
  "PAX3-FOXO1_Targets_ChipSeq"    = "", 
  "Tumor_Suppressor_Genes"        = "", 
  "clinomics_gene_list"           = "" 
)

## function to concatenate files to a single matrix ####
makResultMatrix <- function(x) {
  #print(x)
  data <- tryCatch({
    
    read.table(paste0(x,"/result.bedpe"),sep="\t", header = FALSE, stringsAsFactors = FALSE )
    
  }, error = function(e){
    
    emptyDF
  })
  #print(ncol(data))
  data <- cbind(SampleName = unlist(strsplit(x, "/"))[3], data )
  colnames(data) <- cols
  return(data)
}

### List all the files
bedpeDirs <- list.dirs("./")[-1]
bedpeDirsFusion <- bedpeDirs[grepl("fusion_antigen_out",bedpeDirs)]; head(bedpeDirsFusion)

