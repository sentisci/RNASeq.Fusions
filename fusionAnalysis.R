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


### List all the files ; Remove unconfirmed genes; Remove fusions found in Normal samples ####
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

## Get only completecases 
finalResultMatrixJoin <- finalResultMatrixJoin %>% filter(!Sample.ID %in% c(df_NA_TRUE$Sample.ID))
dim(finalResultMatrixJoin) ; length(unique(finalResultMatrixJoin$Sample.ID))

## Sanity Check
finalResultMatrixJoin %>% filter(grepl("CREM|INO80D", right_gene)) %>% dim()

## Make a key 
finalResultMatrixJoin$key <- paste(finalResultMatrixJoin$left_chr, finalResultMatrixJoin$left_gene, finalResultMatrixJoin$right_chr, 
                                   finalResultMatrixJoin$right_gene, sep="_")
finalResultMatrixJoin$key2 <- paste(finalResultMatrixJoin$left_chr, finalResultMatrixJoin$left_gene, finalResultMatrixJoin$left_position, finalResultMatrixJoin$right_chr, 
                                   finalResultMatrixJoin$right_gene, finalResultMatrixJoin$right_position, finalResultMatrixJoin$Sample.ID, sep="_")
finalResultMatrixJoin.NoDup <- finalResultMatrixJoin[!duplicated(finalResultMatrixJoin$key2), ]
finalResultMatrixJoin.Dup <- finalResultMatrixJoin[duplicated(finalResultMatrixJoin$key2), ]

## Sanity Check
finalResultMatrixJoin.NoDup %>% filter(grepl("CREM|INO80D", right_gene)) %>% dim()

## Count samples with Defuse
finalfusionResultMatrix.Defuse <- finalResultMatrixJoin.NoDup %>% filter( grepl('Defuse', tool)) %>%
                                          dplyr::distinct(Patient.ID, Sample.ID, Case.ID)
dim(finalfusionResultMatrix.Defuse)
## Keep only max spanning read count.

##Read Annotation Data
Annotation  <- readRDS("C:/Users/sindiris/R Scribble/Annotation RDS/annotation_ENSEMBL_gene.RDS") ; Annotation$Start <- as.numeric(as.character(Annotation$Start)) ; Annotation$End <- as.numeric(as.character(Annotation$End))

##Filter Bogus Genes
geneFusionRight <-  as.character( unique(finalResultMatrixJoin.NoDup$right_gene) )
rightGenePair   <- geneFusionRight[which(!geneFusionRight %in% Annotation$GeneName)]
geneFusionLeft  <- as.character(unique(finalResultMatrixJoin.NoDup$left_gene) )
leftGenePair    <- geneFusionLeft[which(!geneFusionLeft %in% Annotation$GeneName)]

## Remove disgusting samples . "NCIEWS5000muscle_T_D23N9ACXX" , "NCI0228normal_T_C3JV2ACXX" . I dont understand why tumor patient normal is used in
## RNASeq cohort experiment
finalResultMatrixJoin.NoDup.NoPatientNormal <-  finalResultMatrixJoin.NoDup %>% filter()
dim(finalResultMatrixJoin.NoDup.NoPatientNormal)
## Sanity Check
finalResultMatrixJoin.NoDup.NoPatientNormal %>% filter(grepl("CREM|INO80D", right_gene)) %>% dim()

##Filter Fusion File
fusionFile.NS  <- finalResultMatrixJoin.NoDup.NoPatientNormal %>% filter( DIAGNOSIS.Alias %in% c("NS") ) ; dim(fusionFile.NS )
finalResultMatrixJoin.NoDup.NoNS <- finalResultMatrixJoin.NoDup.NoPatientNormal %>% filter(! key %in% fusionFile.NS$key ); dim(finalResultMatrixJoin.NoDup.NoNS)
## Remove bogus genes
fusionFile  <- finalResultMatrixJoin.NoDup.NoNS %>% filter( !(left_gene %in% leftGenePair) ) %>%  filter( !(right_gene %in% rightGenePair) ) ; dim(fusionFile)
fusionFile <- data.table(fusionFile)
## Sanity Check
fusionFile %>% filter(grepl("CREM|INO80D", right_gene)) %>% dim()
dim(fusionFile)


### Extra decondensing and filtering For Fusion ONLY ####

### Step 1 Keep the selected tier
fusionFile <- fusionFile[ grepl("Tier 1.1|Tier 1.2|Tier 1.3|Tier 2.1", var_level) ]; 
## Sanity Check
fusionFile %>% filter(grepl("CREM|INO80D", right_gene)) %>% dim()
dim(fusionFile)

### Step 2 decondensing
newCols <- c("SPTool1", "SPTool2", "SPTool3")
fusionFileSplits.v1 <- fusionFile %>% tidyr::separate( spanreadcount, newCols , " " )
fusionFileSplits.v1[is.na(fusionFileSplits.v1)] <- 0 
fusionFileSplits.v1[,(newCols):= lapply(.SD, as.numeric), .SDcols = newCols]
## Sanity Check
fusionFileSplits.v1 %>% filter(grepl("CREM|INO80D", right_gene)) %>% dim()
dim(fusionFileSplits.v1); View(fusionFileSplits.v1)

### Step 3 Remove fusions called by Star-fusions or FusionCatcher
fusionFileFilt.v2 <- fusionFileSplits.v1[ !grepl("^FusionCatcher$|^STAR-fusion$", tool) ] ;
fusionFileFilt.v2 %>% filter(grepl("CREM|INO80D", right_gene)) %>% dim()
dim(fusionFileFilt.v2); View(fusionFileFilt.v2)

### Step 4 Keep if any two or rmore callers regardless of spanning reads
fusionFileFilt.v3 <- fusionFileFilt.v2[ grepl("^tophatFusion$", tool) &  SPTool1 >= 5 |
                                             SPTool1 > 0  & SPTool2 > 0 |
                                             SPTool1 > 0  & SPTool3 > 0 |
                                             SPTool2 > 0  & SPTool3 > 0 
                                        ];
fusionFileFilt.v3 %>% filter(grepl("CREM|INO80D", right_gene)) %>% dim()
dim(fusionFileFilt.v3); View(fusionFileFilt.v3)

### Step 5 Keep Tier 2.1 fusions with spanning reads >= 10
fusionFileFilt.v3.ToRemove <- fusionFileFilt.v3[ var_level %in% c("Tier 2.1") &
                                !(SPTool1 >= 10 & SPTool2 >= 10 & SPTool3 >= 10) & 
                                !grepl("^tophatFusion$", tool) &  SPTool1 >= 5   &
                                !grepl("^FusionCatcher STAR-fusion STAR-fusion tophatFusion$", tool) &  SPTool3 < 5 &
                                !grepl("^FusionCatcher tophatFusion$", tool) &  SPTool2 < 5 &
                                !grepl("^STAR-fusion tophatFusion$", tool) &  SPTool2 < 5
                              ]
dim(fusionFileFilt.v3.ToRemove); View(fusionFileFilt.v3.ToRemove)

### Step 6 Final filter to generate the file
fusionFileFilt.v4 <- fusionFileFilt.v3[ !(fusionFileFilt.v3$key %in% fusionFileFilt.v3.ToRemove$key) ]
fusionFileFilt.v4 %>% filter(grepl("CREM|INO80D", right_gene)) %>% dim()
dim(fusionFileFilt.v4); View(fusionFileFilt.v4)

write.table(finalResultMatrixJoin.NoDup, "../FinalFusionResultMatrix.txt", sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE )
write.table(fusionFileFilt.v4 , "../FinalFilteredfusionResultMatrix.txt", sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)



### Extra decondensing and filtering for Fusion Neoantigens ####

### Step 1 decondensing
newCols <- c("SPTool1", "SPTool2", "SPTool3")
fusionFileSplits.v1 <- fusionFile %>% tidyr::separate( spanreadcount, newCols , " " )
fusionFileSplits.v1[is.na(fusionFileSplits.v1)] <- 0 
fusionFileSplits.v1[,(newCols):= lapply(.SD, as.numeric), .SDcols = newCols]
## Sanity Check
fusionFileSplits.v1 %>% filter(grepl("^CREM$|^INO80D$|^PAX3$|^PAX7$|^EWSR1$|^FLI1$|^SSX18$", left_gene)) %>% dim()
dim(fusionFileSplits.v1); View(fusionFileSplits.v1)

### Step 2 Remove fusions called by Star-fusions or FusionCatcher
fusionFileFilt.v2 <- fusionFileSplits.v1[ !grepl("^FusionCatcher$|^STAR-fusion$", tool) ] ;
fusionFileFilt.v2 %>% filter(grepl("^CREM$|^INO80D$|^PAX3$|^PAX7$|^EWSR1$|^FLI1$|^SSX18$", right_gene)) %>% dim()
dim(fusionFileFilt.v2); View(fusionFileFilt.v2)

### Step 3 Keep if any two or rmore callers regardless of spanning reads
fusionFileFilt.v3 <- fusionFileFilt.v2[ grepl("^tophatFusion$", tool) &  SPTool1 >= 5 |
                                          SPTool1 > 0  & SPTool2 > 0 |
                                          SPTool1 > 0  & SPTool3 > 0 |
                                          SPTool2 > 0  & SPTool3 > 0  ];
dim(fusionFileFilt.v3); View(fusionFileFilt.v3)

#write.table(fusionFileFilt.v3, "../FinalFusionResultMatrixForNeoantigens.txt", sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE )

## merge this with the final fusion bedpe file.Remove any fusions in bedpe file that aren't present in the above file.
## file "finalGroupBy2.Filt" from "fusionNeoantigenAnalysis.R"

neoepitopesFromFusions <- finalGroupBy2.Filt %>% data.table() ; dim(neoepitopesFromFusions) ; View(neoepitopesFromFusions)
validFusionsForNeoantigen <- fusionFileFilt.v3 %>% dplyr::filter(type=="In-frame")  ; dim(validFusionsForNeoantigen)  ; View(validFusionsForNeoantigen)

## Make Primary ID
colsToChange <- c("start.5p", "end.5p", "start.3p", "end.3p")
neoepitopesFromFusions[,(colsToChange):= lapply(.SD, as.character), .SDcols = colsToChange]
neoepitopesFromFusions[,(colsToChange):= lapply(.SD, as.numeric), .SDcols = colsToChange]

# ## Make new ID
neoepitopesFromFusionsID <- neoepitopesFromFusions %>% dplyr::mutate(left_position = ifelse(neoepitopesFromFusions$start.5p > neoepitopesFromFusions$end.5p,
                                                           neoepitopesFromFusions$start.5p,
                                                           neoepitopesFromFusions$end.5p),
                                         right_position = ifelse(neoepitopesFromFusions$start.3p > neoepitopesFromFusions$end.3p,
                                                                 neoepitopesFromFusions$start.3p,
                                                                 neoepitopesFromFusions$end.3p) )
# ## Make PrimaryID for neoepitopes
neoepitopesFromFusionsID %<>% dplyr::mutate( primaryID = paste0("chr", neoepitopesFromFusionsID$chr.5p,
                                                                "_", neoepitopesFromFusionsID$left_position,
                                                                "_", neoepitopesFromFusionsID$fusion, "_",
                                                                "chr", neoepitopesFromFusionsID$chr.3p, "_",
                                                                neoepitopesFromFusionsID$right_position) ) %>% data.table()
# View(neoepitopesFromFusionsID); dim(neoepitopesFromFusionsID)
## There is an issue in INTEGRATENEO output. Some times the coordinates (Start and End) positions are reduced by one. Just to get those fusions Adding one to
## the INTEGRATENEO output.
neoepitopesFromFusionsID.AddOne <- neoepitopesFromFusionsID %>% dplyr::mutate(left_position = gsub('.{1}$', '', left_position),
                                                                              right_position = gsub('.{1}$', '', right_position))
neoepitopesFromFusionsID.AddOne %<>% dplyr::mutate( primaryID = paste0("chr", neoepitopesFromFusionsID.AddOne$chr.5p, "_",
                                                                              neoepitopesFromFusionsID.AddOne$left_position, "_",
                                                                              neoepitopesFromFusionsID.AddOne$fusion, "_",
                                                                       "chr", neoepitopesFromFusionsID.AddOne$chr.3p, "_",
                                                                              neoepitopesFromFusionsID.AddOne$right_position) ) %>% data.table()
# View(neoepitopesFromFusionsID.AddOne); dim(neoepitopesFromFusionsID.AddOne)

## Make PrimaryID for fusions
validFusionsForNeoantigenID <- validFusionsForNeoantigen %>% dplyr::mutate(primaryID.AddOne =  paste0(validFusionsForNeoantigen$left_chr, "_",
                                                                                               gsub('.{1}$', '', validFusionsForNeoantigen$left_position), "_",
                                                                                              paste0(validFusionsForNeoantigen$left_gene, ">>", validFusionsForNeoantigen$right_gene), "_",
                                                                                               validFusionsForNeoantigen$right_chr, "_",
                                                                                              gsub('.{1}$', '', validFusionsForNeoantigen$right_position)) ) %>% data.table()
View(validFusionsForNeoantigenID) 

################### Sanity Check #################################
neoepitopesFromFusionsID[grep("ASPSCR1", neoepitopesFromFusionsID$primaryID)]
neoepitopesFromFusionsID[grep("SS18", neoepitopesFromFusionsID$primaryID)]
neoepitopesFromFusionsID[grep("EWSR1", neoepitopesFromFusionsID$primaryID)]
neoepitopesFromFusionsID[grep("PAX7", neoepitopesFromFusionsID$primaryID)]
neoepitopesFromFusionsID[grep("PAX3", neoepitopesFromFusionsID$primaryID)]
neoepitopesFromFusionsID[grep("INO80D", neoepitopesFromFusionsID$primaryID)]
neoepitopesFromFusionsID[grep("CREM", neoepitopesFromFusionsID$primaryID)]

validFusionsForNeoantigenID[grep("ASPSCR1", validFusionsForNeoantigenID$primaryID)]
validFusionsForNeoantigenID[grep("SS18", validFusionsForNeoantigenID$primaryID),]
validFusionsForNeoantigenID[grep("EWSR1", validFusionsForNeoantigenID$primaryID)]
validFusionsForNeoantigenID[grep("PAX7", validFusionsForNeoantigenID$primaryID)]
validFusionsForNeoantigenID[grep("PAX3", validFusionsForNeoantigenID$primaryID),]
validFusionsForNeoantigenID[grep("INO80D", validFusionsForNeoantigenID$primaryID),]
validFusionsForNeoantigenID[grep("CREM", validFusionsForNeoantigenID$primaryID),]
##################################################################

## Filter only valid fusion-neoepitopes
neoepitopesFromFusionsID.AddOne.Filtered <- neoepitopesFromFusionsID.AddOne[primaryID %in% validFusionsForNeoantigenID$primaryID.AddOne]; dim(neoepitopesFromFusionsID.AddOne.Filtered)
neoepitopesFromFusionsID.Filtered.Final <- neoepitopesFromFusionsID.AddOne.Filtered; dim(neoepitopesFromFusionsID.Filtered.Final)
##

#################### Sanity Check ################################
neoepitopesFromFusionsID.Filtered.Final[grep("ASPSCR1", neoepitopesFromFusionsID.Filtered.Final$primaryID)]
neoepitopesFromFusionsID.Filtered.Final[grep("SS18", neoepitopesFromFusionsID.Filtered.Final$primaryID)]
neoepitopesFromFusionsID.Filtered.Final[grep("EWSR1", neoepitopesFromFusionsID.Filtered.Final$primaryID)]
neoepitopesFromFusionsID.Filtered.Final[grep("PAX7", neoepitopesFromFusionsID.Filtered.Final$primaryID)]
neoepitopesFromFusionsID.Filtered.Final[grep("PAX3", neoepitopesFromFusionsID.Filtered.Final$primaryID)]
neoepitopesFromFusionsID.Filtered.Final[grep("INO80D", neoepitopesFromFusionsID.Filtered.Final$primaryID)]
neoepitopesFromFusionsID.Filtered.Final[grep("CREM", neoepitopesFromFusionsID.Filtered.Final$primaryID)]
##################################################################

write.table(neoepitopesFromFusionsID.Filtered, paste0(workDir,"/NoMaxEffinityEpitope.Filtered.Normal.ReadThrough.GTE.3.ValidFusions.txt"), sep="\t", row.names = FALSE, quote = FALSE)



## Count NeoEpitopes per sample
EpitopesGroupBySample <- neoepitopesFromFusionsID.Filtered %>% dplyr::group_by(Sample.Data.ID ) %>% mutate(FusionNeoAntigenCount = n()) %>% 
  dplyr::select(Sample.Data.ID, FusionNeoAntigenCount) %>% 
  dplyr::distinct()

EpitopesGroupBySample.ToPrint <- dplyr::full_join(metaData, EpitopesGroupBySample, by="Sample.Data.ID") %>% 
  mutate(FusionNeoAntigenCountMod = ifelse(is.na(FusionNeoAntigenCount), 0, FusionNeoAntigenCount)) %>%
  dplyr::filter( !is.na(FusionNeoAntigenCount))
dim(EpitopesGroupBySample)
View(EpitopesGroupBySample.ToPrint)
write.table(EpitopesGroupBySample.ToPrint, paste0(workDir,"/EpitopesGroupBySampleFilt.NoReadThrough.noMaxAffinityHLA.ValidFusions.txt"), sep="\t", row.names = FALSE, quote = FALSE)

## Make plots for Differential gene expression analysis




