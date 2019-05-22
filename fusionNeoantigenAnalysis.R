rm(list = ls())
workDir = "/Users/sindiris/R Scribble/RNASeq.Fusion.data"
setwd(paste0(workDir,"/NeoantigensFromFusions/"))

## Read well annotated DataFrame
metaData <- read.csv(paste0(workDir,"/MetadataMapper.v3.txt"), header = T, sep = "\t", stringsAsFactors = FALSE); 
# metaData <- metaData %<>% dplyr::rename(SampleName = Biowulf_SampleName)
head(metaData)

### Standard Colnames
cols <- c("SampleName", "chr.5p", "start.5p", "end.5p", "chr.3p", "start.3p", "end.3p", "fusion", "TierOfFusion", "leftStrand", "RightStand",
          "MaxSpanningReads", "Epitope.Sequence", "Epitope.Affinity", "HLAallele", "HLAcategory", "HLA.score", "HLA.e-value",
          "HLA.confidence")

### Required functions
makResultBedpeMatrix <- function(x) {
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

## Empty data frame incase file is empty ####
emptyDF <- data.frame( 
  "chr.5p"     = "",
  "start.5p"   = "",
  "end.5p"     ="",
  "chr.3p"     ="",
  "start.3p"   ="",
  "end.3p"     ="",
  "fusion"     = "NF",
  "TierOfFusion" = "",
  "leftStrand" ="",
  "RightStand" = "",
  "MaxSpanningReads" ="",
  "Epitope.Sequence" = "",
  "Epitope.Affinity(nanoMolar)" = "",
  "HLAallele" ="",
  "HLAcategory" = "",
  "HLA.score" ="",
  "HLA.e-value" ="",
  "HLA.confidence" =""
)



### read all the files and build a matrix ####
finalResultBedpeMAtrix <- rbindlist(lapply(bedpeDirsFusion, makResultBedpeMatrix), fill = TRUE);
dim(finalResultBedpeMAtrix);head(finalResultBedpeMAtrix)
finalResultBedpeMAtrix <- finalResultBedpeMAtrix %>%  dplyr::rename(Sample.Data.ID = SampleName )
finalResultBedpeMAtrix$Epitope.Affinity <- as.numeric(as.character(finalResultBedpeMAtrix$Epitope.Affinity))
dim(finalResultBedpeMAtrix);head(finalResultBedpeMAtrix)

## Adding metadata
finalResultBedpeMAtrixJoin <- dplyr::left_join(finalResultBedpeMAtrix, metaData[,
                                              c("Sample.Data.ID", "Sample.ID.Alias", "DIAGNOSIS.Alias", "LIBRARY_TYPE")], by="Sample.Data.ID")
head(finalResultBedpeMAtrixJoin)

## Sanity Check1 : To evaluate the join // Check if join produced any NA
df <- finalResultBedpeMAtrixJoin %>% data.frame()
df$Any_NA <- apply(df[,grep(paste(c("Sample.Data.ID", "Sample.ID.Alias", "DIAGNOSIS.Alias", "LIBRARY_TYPE"),collapse = "|"), names(df) )],1, function(x) anyNA(x))
df_NA_TRUE <- df[which(df$Any_NA == "TRUE"),] %>% dplyr::distinct(Sample.Data.ID,DIAGNOSIS.Alias,LIBRARY_TYPE);dim(df_NA_TRUE); View(df_NA_TRUE)

## Sanity Check2 : match the samples & patient in the vcf files with the metadata

## Remove CellLines
finalResultBedpeMAtrixJoin <- finalResultBedpeMAtrixJoin %>% filter(!LIBRARY_TYPE %in% c("CellLine"))
metaDataNoCellLine         <- metaData %>% filter(!LIBRARY_TYPE %in% c("CellLine"))

## Check for the sample match / mismatch
samplesFromBedpe <- unique(finalResultBedpeMAtrixJoin$Sample.Data.ID) ; samplesFromMetadata <- unique(metaDataNoCellLine$Sample.Data.ID)
patientsEqual = length(samplesFromBedpe) == length(samplesFromMetadata)
if(!patientsEqual) { 
    
    stop(paste0("Following samples data couldn't be retirieved from files",
                paste(samplesFromMetadata[which( ! samplesFromMetadata %in% samplesFromBedpe)],collapse = " ; ") )) 
}

### Write the data
write.table(finalResultBedpeMAtrixJoin, paste0(workDir,"allNeoFusedEpitopes.txt"), sep="\t", row.names = FALSE, quote = FALSE)

## Grouping
## Normal Fusion genes
NormalEpitopes <- finalResultBedpeMAtrixJoin %>% dplyr::filter( grepl("Normal", LIBRARY_TYPE) ) %>% 
                                                 dplyr::filter( !Sample.Data.ID %in% c("Sample_NCIEWS5000muscle_T_D23N9ACXX", "Sample_NCI0228normal_T_C3JV2ACXX"))
dim(NormalEpitopes); head(NormalEpitopes)

## For.Each.Sample.And.Fusion.maxEffinityEpitope.Filtered.Normal.Fusions
TumorEpitopes <- finalResultBedpeMAtrixJoin %>% dplyr::filter( grepl("Tumor", LIBRARY_TYPE)) %>% 
                                                dplyr::filter( !fusion %in% c(NormalEpitopes$fusion)) %>% 
                                                dplyr::filter( !fusion %in% c("NF")); dim(TumorEpitopes); head(TumorEpitopes)
write.table(TumorEpitopes, paste0(workDir,"/NoMaxEffinityEpitope.Filtered.Normal.Fusions.txt"), sep="\t", row.names = FALSE, quote = FALSE)

TumorEpitopes$MaxSpanningReads <- as.numeric(as.character(TumorEpitopes$MaxSpanningReads))

## Max affinity epitopes
TumorEpitopes_MaxAffinity <- TumorEpitopes %>% 
                            dplyr::group_by(Sample.Data.ID, fusion, Epitope.Sequence) %>%
                            dplyr::mutate(maxAffinity = max(Epitope.Affinity)) %>%
                            filter(Epitope.Affinity == maxAffinity) #%>%
                            #filter( !(Epitope.Sequence %in% unique(NormalEpitopes$Epitope.Sequence)) )
write.table(TumorEpitopes_MaxAffinity, paste0(workDir,"/MaxEffinityEpitope.Filtered.Normal.Fusions.txt"), sep="\t", row.names = FALSE, quote = FALSE)

## For.Each.Sample.And.Fusion.maxEffinityEpitope.Filtered.Normal.Fusions + Stats for each epitopes at Diagnosis and Sample Level
finalGroupBy2 <- TumorEpitopes %>% dplyr::group_by(fusion )           %>%
  mutate(SamplesCountFusion   = length( unique(as.character(Sample.Data.ID))) ) %>%
  mutate(DiagnosisListFusion      = paste( unique(as.character(DIAGNOSIS.Alias)), collapse="," ) ) %>%
  mutate(DiagnosisCountFusion = length( unique(as.character(DIAGNOSIS.Alias))) ) %>%
  ungroup() %>% dplyr::group_by(fusion, Epitope.Sequence)           %>%
  mutate(SamplesCountEpitope   = length( unique(as.character(Sample.Data.ID))) ) %>%
  mutate(DiagnosisListEpitope      = paste( unique(as.character(DIAGNOSIS.Alias)), collapse="," ) ) %>%
  mutate(DiagnosisCountEpitope = length( unique(as.character(DIAGNOSIS.Alias))) ) %>% 
  mutate(MaxSpanningReads = as.numeric(as.character(MaxSpanningReads)))
dim(finalGroupBy2)
head(finalGroupBy2)
View(finalGroupBy2)
#write.table(finalGroupBy2, "../For.Each.Sample.And.Fusion.maxEffinityEpitope.Filtered.Normal.Fusions.Epitopes.STATS.txt", sep="\t", row.names = FALSE, quote = FALSE)

# ### Filter the data frame to keep only valid read-throughs ####
# finalGroupBy2ReadThrough <- finalGroupBy2 %>% dplyr::filter( chr.5p == chr.3p )
# CheckReadThroughGenes <- function(x){
#   
#   geneList <- unlist(strsplit(x[8], ">>")[[1]])
#   rightGeneIndex <- which(EnsemblAnnotationSorted$GeneName == geneList[2])
#   #print(paste(x[8], rightGeneIndex, length(geneList) ))
#   
#   data <- tryCatch({
#     
#     if(EnsemblAnnotationSorted[rightGeneIndex , "GeneName"] == geneList[2] ) { return(FALSE)}
#     else { return(TRUE)}
#     
#   }, error = function(e){
#     
#     return(FALSE)
#     
#   })
#   
# }
# 
# finalGroupBy2ReadThrough <-  finalGroupBy2ReadThrough %>% data.frame()
# binaryVector <- apply(finalGroupBy2ReadThrough, 1, CheckReadThroughGenes)
# finalGroupBy2ReadThrough <- cbind(finalGroupBy2ReadThrough, binaryVector) 
# finalGroupBy2ReadThroughFocal <- finalGroupBy2ReadThrough %>% filter(binaryVector)


 


### Notes: Following are the two ways to manage read-through fusions. Use these files (before grouping) to add Neoantigen filters.
## Merge this file with the fusions file from "FusionAnalysis.R" to implement neoantigen filters.

################################## METHOD 1 #######################
### Add more filters to the above step USING ONLY TRANSGENIC FUSIONS ####
finalGroupBy2.Filt <- finalGroupBy2 %>% dplyr::filter( !(chr.5p == chr.3p) )
dim(finalGroupBy2.Filt)
#head(finalGroupBy2)
#write.table(finalGroupBy2.Filt, "../For.Each.Sample.And.Fusion.maxEffinityEpitope.Filtered.Normal.Fusions.Epitopes.STATS.Filt.NoReadThrough.txt", sep="\t", row.names = FALSE, quote = FALSE)
write.table(finalGroupBy2.Filt, paste0(workDir,"/NoMaxEffinityEpitope.Filtered.Normal.NoReadThrough.txt"), sep="\t", row.names = FALSE, quote = FALSE)

## Count NeoEpitopes per sample
EpitopesGroupBySample <- finalGroupBy2.Filt %>% dplyr::group_by(Sample.Data.ID ) %>% mutate(FusionNeoAntigenCount = n()) %>% 
  dplyr::select(Sample.Data.ID, FusionNeoAntigenCount) %>% 
  dplyr::distinct()

EpitopesGroupBySample.ToPrint <- dplyr::full_join(metaData, EpitopesGroupBySample, by="Sample.Data.ID") %>% 
  mutate(FusionNeoAntigenCountMod = ifelse(is.na(FusionNeoAntigenCount), 0, FusionNeoAntigenCount)) %>%
  dplyr::filter( !is.na(FusionNeoAntigenCount))
dim(EpitopesGroupBySample)
View(EpitopesGroupBySample.ToPrint)
#write.table(EpitopesGroupBySample.ToPrint, "../EpitopesGroupBySampleFilt.NoReadThrough.txt", sep="\t", row.names = FALSE, quote = FALSE)
write.table(EpitopesGroupBySample.ToPrint, paste0(workDir,"/EpitopesGroupBySampleFilt.NoReadThrough.noMaxAffinityHLA.txt"), sep="\t", row.names = FALSE, quote = FALSE)

################################## METHOD 2 #######################

### Add more filters to the above step USING readthroughs with Spanning reads >=3
finalGroupBy2.Filt <- finalGroupBy2 %>% dplyr::filter( !(chr.5p == chr.3p | MaxSpanningReads < 3) )
dim(finalGroupBy2.Filt)
#head(finalGroupBy2)
#write.table(finalGroupBy2.Filt, "../For.Each.Sample.And.Fusion.maxEffinityEpitope.Filtered.Normal.Fusions.Epitopes.STATS.Filt.SP3.txt", sep="\t", row.names = FALSE, quote = FALSE)
write.table(finalGroupBy2.Filt, paste0(workDir,"/NoMaxEffinityEpitope.Filtered.Normal.NoReadThrough.SR.GTE.3.txt"), sep="\t", row.names = FALSE, quote = FALSE)


## Count NeoEpitopes per sample 
EpitopesGroupBySample <- finalGroupBy2.Filt %>% dplyr::group_by(Sample.Data.ID ) %>% mutate(FusionNeoAntigenCount = n()) %>% 
  dplyr::select(Sample.Data.ID, FusionNeoAntigenCount) %>% 
  dplyr::distinct()

EpitopesGroupBySample.ToPrint <- dplyr::full_join(metaData, EpitopesGroupBySample, by="Sample.Data.ID") %>% 
  mutate(FusionNeoAntigenCountMod = ifelse(is.na(FusionNeoAntigenCount), 0, FusionNeoAntigenCount))%>%
  dplyr::filter( !is.na(FusionNeoAntigenCount))
dim(EpitopesGroupBySample)
View(EpitopesGroupBySample.ToPrint)
write.table(EpitopesGroupBySample.ToPrint, paste0(workDir,"/EpitopesGroupBySampleFilt.NoReadThrough.NomaxEffinityEpitope.SP3.txt"), sep="\t", row.names = FALSE, quote = FALSE)


