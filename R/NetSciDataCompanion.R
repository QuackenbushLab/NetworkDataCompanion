NetSciDataCompanion=setRefClass("NetSciDataCompanion",

         fields = list(TCGA_purities= "data.frame",
                       clinical_patient_data = "data.frame",
                       project_name = "character",
                       gene_mapping = "data.frame",
                       sample_type_mapping = "data.frame"),
         methods = list(


           ## Extract experiment specific information and metadata from ranged summarized experiment object
           ## Returns a named list with rds_sample_info corresponding to meta information about the samples (columns)
           ##                       and rds_gene_info corresponding to meta information about genes (rows)
           ## 20220913 man page done
           extractSampleAndGeneInfo = function(expression_rds_obj){
             return(list(rds_sample_info=as.data.frame(colData(test_exp_rds)), rds_gene_info=as.data.frame(rowRanges(test_exp_rds))))
           },


           ## Maps two sets of barcodes
           ## There are 4 different pieces of information returned in a named list that are all useful depending on the context they are used in
           ## is_inter1 is an indicator (boolean) vector of the same length as bc1 that indicates which elements of bc1 are present in bc2
           ## idcs1 indicate where to find each barcode of bc1 in bc2, NA if missing. That is, bc1[i] == bc2[idcs1[i]] (if idcs1[i] != NA)
           ## The same information is provided for bc2
           ## For example, if you want to map experiment 1 on experiment to, keeping only the information for samples that are present in both,
           ## and reordering the first experiment to match the samples of the second, you can do
           ## exp1[,is_inter1]                --- this will remove samples that are not in exp2
           ## exp2[,idcs1]                    --- this will remove samples that are not in exp1 and reorder to match exp1
           ## 20220920 man page done
           mapBarcodeToBarcode = function(bc1, bc2){
             if(class(bc1) != "character" | class(bc2) != "character"){
               stop("Error: barcodes need to be vectors of strings")
             }
             m1 <- match(bc1, bc2)
             m1 <- m1[!is.na(m1)]
             m2 <- match(bc2, bc1)
             m2 <- m2[!is.na(m2)]
             return(list(is_inter1=(bc1 %in% bc2), idcs1=m1, is_inter2=(bc2 %in% bc1), idcs2=m2))
           },

           ## A convenience wrapper function for mapBarcodeToBarcode that applies the function directly to two data frames
           ## returns a list of the two argument data frames, intersected, and the second frame ordered to match the first
           ## NOTE: Ordering is done based on columns, which are expected to be named by TCGA barcodes
           ## 20220920 man page done
           filterBarcodesIntersection = function(exp1, exp2){
             if("data.frame" %in% class(exp1) & "matrix" %in% class(exp1) ){
               stop("Error: argument 1 needs to be data.frame or matrix")
             }
             if("data.frame" %in% class(exp2) & "matrix" %in% class(exp2) ){
               stop("Error: argument 2 needs to be data.frame or matrix")
             }
             map <- mapBarcodeToBarcode(extractSampleOnly(colnames(exp1)), extractSampleOnly(colnames(exp2)))
             return(list("mappedExp1"=exp1[,map$is_inter1], "mappedExp2"=exp2[,map$idcs1]))
           },


           ## Computes the log TPM normalization based on an expression RDS object
           ## Returns a named list with the count data.frame (useful for duplicate filtering based on sequencing depth, see filterDuplicatesSeqDepth)
           ##                               TPM data.frame (useful for TPM based filtering, see filterGenesByTPM)
           ##                and the actual logTPM which corresponds to log(TPM + 1)
           ## 20220920 man page done
           logTPMNormalization = function(expression_rds_obj){
             if(class(expression_rds_obj) != "RangedSummarizedExperiment"){
               stop("Error: expression matrices need to be an RSE object")
             }
             assays(expression_rds_obj)$counts <- recount3::transform_counts(expression_rds_obj)
             assays(expression_rds_obj)$TPM <- recount::getTPM(expression_rds_obj)
             return(list(counts=assays(expression_rds_obj)$counts,
                         TPM=assays(expression_rds_obj)$TPM,
                         logTPM=log(assays(expression_rds_obj)$TPM + 1)))
           },

           # Computes the log transformed CPM normalization based on an expression RDS object
           # This is following the recipe provided by edgeR package to get TMM valyes
           ## Returns a named list with the count data.frame (useful for duplicate filtering based on sequencing depth, see filterDuplicatesSeqDepth)
           ##                               CPM data.frame (useful for CPM based filtering, see filterGenesByCPM)
           ##                and the actual logCPM which corresponds to log(CPM + 1)
           logCPMNormalization = function(expression_rds_obj){
             if(class(expression_rds_obj) != "RangedSummarizedExperiment"){
               stop("Error: expression matrices need to be an RSE object")
             }
             # get counts through recount
             assays(expression_rds_obj)$counts <- recount3::transform_counts(expression_rds_obj)
             # apply edgeR function to get differential gene lists
             dge <- DGEList(assays(expression_rds_obj)$counts)
             # get the normalizing factors from that list
             dge <- calcNormFactors(dge)
             # get the log tranformed cpms (aka: get TMM)
             assays(expression_rds_obj)$CPM <- cpm(dge, log=F)
             return(list(counts=assays(expression_rds_obj)$counts,
                         CPM=assays(expression_rds_obj)$CPM,
                         logCPM=log(assays(expression_rds_obj)$CPM+1)))
           },

           #### more methods go here

           # maybe have this presaved in class
           extractSampleOnly = function(TCGA_barcodes){
             return(sapply(TCGA_barcodes, substr, 1, 12))
           },

           extractVialOnly = function(TCGA_barcodes){
              return(sapply(TCGA_barcodes, substr, 1, 16))
           },

           extractSampleType = function(TCGA_barcodes){
              return(sapply(TCGA_barcodes, substr, 14, 15))
           },

           findDuplicates = function(TCGA_barcodes){
              return(duplicated(extractVialOnly(TCGA_barcodes)))
           },

           mapUUIDtoTCGA = function(UUID, useLegacy = F){
              if(class(UUID) != "character"){
                stop("Error: Expected UUID argument to be vector of strings")
              }
              info = files(legacy = useLegacy) %>%
               GenomicDataCommons::filter( ~ file_id %in% UUID) %>%
               GenomicDataCommons::select('cases.samples.submitter_id') %>%
               results_all()
              # The mess of code below is to extract TCGA barcodes
              # id_list will contain a list (one item for each file_id)
              # of TCGA barcodes of the form 'TCGA-XX-YYYY-ZZZ'
              id_list = lapply(info$cases,function(a) {
               a[[1]][[1]][[1]]})
              # so we can later expand to a data.frame of the right size
              barcodes_per_file = sapply(id_list,length)
              # sort to match input UUID order
              file_id <- rep(ids(info),barcodes_per_file)
              reord <- match(UUID, file_id)
              # And build the data.frame
              return(data.frame(file_id = file_id[reord],
                               submitter_id = unlist(id_list)[reord]))
           },


           # rangeUp and rangeDown should both be non-negative numbers
           # indicating the distance to look upstream and downstream
           # from a probe for a TSS.
           # the function will convert this to negative number.
           # note this implementation ONLY looks for a TSS and not
           # for gene body or other regions.
           # long form: returns one row per gene, so if a probe  maps
           # to the TSS of two different genes, each of those gets a row

           mapProbesToGenes = function(probelist,
                                       rangeUp = 200,
                                       rangeDown = 0,
                                       localManifestPath=NA,
                                       longForm = F){

             if(is.na(localManifestPath))
             {
               print("[NetSciDataCompanion::mapProbesToGenes] Sourcing probe annotation from https://zwdzwd.github.io/InfiniumAnnotation ...")

               # source hg38 with gencode 36 from https://zwdzwd.github.io/InfiniumAnnotation
               download.file('https://zhouserver.research.chop.edu/InfiniumAnnotation/20210615/HM450/HM450.hg38.manifest.gencode.v36.tsv.gz',
                             destfile = "./HM450.hg38.manifest.gencode.v36.tsv.gz")

               # unzip
               system2(command="gunzip",args=c("./HM450.hg38.manifest.gencode.v36.tsv.gz"))

               # load into memory
               manifest = data.frame(fread("./HM450.hg38.manifest.gencode.v36.tsv",sep="\t",header=T))

               # remove from local storage
               system2(command="rm",args="./HM450.hg38.manifest.gencode.v36.tsv")
             }

             if(!is.na(localManifestPath))
             {
               print(paste("[NetSciDataCompanion::mapProbesToGenes] Loading probe file from:",localManifestPath))
               manifest = read.table(localManifestPath,sep="\t",header=T)
             }

             # get indices matching probes
             smallManifest = manifest %>% dplyr::filter(probeID %in% probelist)
             rm(manifest)
             gc()

             # define empty map
             mymap = matrix(rep(NA,4*nrow(smallManifest)),ncol=4)
             mymap[,1] = probelist
             mymap = as.data.frame(mymap)

             # iterate through map with for loop
             # please feel free to vectorize this etc
             for(i in 1:nrow(smallManifest))
             {
               if(i %% 10000 == 0) print(paste("[NetSciDataCompanion::mapProbesToGenes] Processing probe number:",i))

               x = smallManifest[i,]
               genes = str_split(x$geneNames,";",simplify=T)
               tssDist = as.numeric(str_split(x$distToTSS,";",simplify=T))
               ensemblIDs = str_split(x$transcriptIDs,";",simplify=T)
               inRegion = which(tssDist > -1*rangeUp & tssDist < rangeDown)
               if(length(inRegion) > 0)
               {
                 genesInRegion = genes[inRegion]
                 ensemblInRegion = ensemblIDs[inRegion]
                 tssDist = tssDist[inRegion]
                 mymap[i,] =c (x$probeID,
                               paste(genesInRegion,collapse=";"),
                               paste(ensemblInRegion,collapse=";"),
                               paste(tssDist,collapse=";"))
               }
             }

             colnames(mymap) = c("probeID","geneName","ensemblID","distToTSS")

             if(!longForm)
              return(mymap)
             else
             {
               # find any probe that has more than one gene mapped to it
               if(length(mymap[grep(";",mymap[,2])]) == 0)
                 return(mymap)

               doubleGenes = data.frame(mymap[grep(";",mymap[,2]),])
               for(i in 1:nrow(doubleGenes))
               {
                 if(i %% 1000 == 0) print(i)
                 # get every split gene
                 theseSplitGenes = str_split(doubleGenes[i,2],";",simplify=T)[1,]
                 theseSplitEns = str_split(doubleGenes[i,3],";",simplify=T)[1,]
                 theseSplitTSS = str_split(doubleGenes[i,4],";",simplify=T)[1,]
                 splitInfo = data.frame("probeID"=doubleGenes[i,1],
                                        "geneName"=theseSplitGenes,
                                        "ensemblID"=theseSplitEns,
                                        "distToTSS"=theseSplitTSS)
                 mymap = rbind.data.frame(mymap,splitInfo)
               }

               # now remove all the original entries that had > 1 gene
               mymap_long = mymap[-grep(";",mymap$geneName),]
               # this map will have multiple rows for a single geneName
               # in some cases (e.g. splice isoform, different ensemblIDs)
               # handling this is a downstream decision, as the map
               # from this point is "long form" in the sense that every
               # row has only one gene
               return(mymap_long)
             }
           },

           # Function to map to probes to a gene-level measurement
           # probe_gene_map is in the format output from the mapProbesToGenes function
           # not all tfGenes need to be in probe_gene_map, but if none are, then this is meaningless
           probeToMeanTFMethylation = function(methylation_betas, probe_gene_map, tfGenes){

             # merge probe_gene_map with beta values
             # use a left join to keep only probes that mapped to genes of interest
             mappedBetas = left_join(probe_gene_map, methylation_betas, by="probeID")

             betaMeans = matrix(NA,nrow = length(tfGenes), ncol = ncol(methylation_betas)-1)
             betaSDs = matrix(NA,nrow = length(tfGenes), ncol = ncol(methylation_betas)-1)

             for(i in 1:length(tfGenes))
             {
               if(i %% 300 == 0) print(i)
               thisGene = tfGenes[i]
               theseProbes = probe_gene_map %>% dplyr::filter(gene == thisGene) %>% dplyr::select(probeID)
               theseBetas = mappedBetas %>% dplyr::filter(probeID %in% theseProbes$probeID) %>% dplyr::select(-c(probeID,gene,ensemblID,distToTSS))
               # any probe that was missed just doesn't contribute to the average
               betaMeans[i,] = apply(theseBetas,2,mean,na.rm=T)
             }

             betaMeansDF = data.frame(t(betaMeans))
             colnames(betaMeansDF) = tfGenes
             row.names(betaMeansDF) = names(methylation_betas)[-1]
             return(betaMeansDF)
           },
           # Input to convertBetaToM is a vector of methylation betas
           # User should use this function with `apply` to convert a matrix
           # 20220920 man page done
           convertBetaToM = function(methylation_betas){
              M = log2(methylation_betas/(1-methylation_betas))
              return(M)
           },

           ## Filter out all duplicates based on sequencing depth
           ## Returns indices about which samples to KEEP
           ## 20220920 man page done
           filterDuplicatesSeqDepth = function(expression_count_matrix){
             sample_barcodes <- extractSampleOnly(colnames(expression_count_matrix))
             seq_depth <- colSums(expression_count_matrix)
             duplicate_throwout <- rep(NA, ncol(expression_count_matrix))
             for (idx in 1:ncol(expression_count_matrix))
             {
               if (is.na(duplicate_throwout[idx]))
               {
                 ## find all vials and replicates of current sample
                 rep_idcs <- sample_barcodes[idx] == sample_barcodes
                 ## get those with highest seq depth
                 max_sample_idx <- which.max(seq_depth * rep_idcs)
                 ## throw out all but maximum one
                 duplicate_throwout[rep_idcs] <- T
                 duplicate_throwout[max_sample_idx] <- F
               }
             }
            return(which(!duplicate_throwout))
           },

           ## Filter out all duplicates based on sequencing depth, take random one if no info on seq depth for all vials
           ## Returns indices in given tcga barcodes to KEEP
           ## 20220920 man page done
           filterDuplicatesSeqDepthOther = function(expression_count_matrix, tcga_barcodes){
             sample_vials_ge <- extractVialOnly(colnames(expression_count_matrix))
             seq_depth <- colSums(expression_count_matrix)
             duplicate_throwout <- rep(NA, length(tcga_barcodes))
             for (idx in 1:length(tcga_barcodes))
             {
               if (is.na(duplicate_throwout[idx]))
               {
                 ## find all vials and replicates of current barcode
                 rep_idcs <- which(extractSampleOnly(tcga_barcodes[idx]) == extractSampleOnly(tcga_barcodes))
                 rep_vials <- extractVialOnly(tcga_barcodes[rep_idcs])
                 ## match with vials in expression matrix
                 mIdx <- match(rep_vials, sample_vials_ge)
                 ## get matched vial with highest seqdepth, just take first one if no match at all
                 max_sample_idx <- 1
                 curr_max <- 0
                 for (j in 1:length(mIdx))
                 {
                   if (!is.na(mIdx[j]))
                   {
                     if (seq_depth[mIdx[j]] > curr_max)
                     {
                       curr_max <- seq_depth[mIdx[j]]
                       max_sample_idx <- rep_idcs[j]
                     }
                   }
                 }
                 ## throw out all but maximum one
                 duplicate_throwout[rep_idcs] <- T
                 duplicate_throwout[max_sample_idx] <- F
               }
             }
             return(which(!duplicate_throwout))
           },

           ## Filter samples indicated by *TCGA_barcodes* based on the method *method* and threshold *threshold*
           ## Returns a list of indices indicating which samples should be kept
           ## 20220920 Man page done
           filterPurity = function(TCGA_barcodes, method="ESTIMATE", threshold=.6){
             if(class(TCGA_barcodes) != "character"){
               stop("Error: Expected TCGA_barcodes argument to be vector of strings")
             }
             if (!(method %in% c("ESTIMATE",
                                 "ABSOLUTE", "LUMP", "IHC", "CPE")))
             {
               stop("Error: Expected method name should be ESTIMATE, ABSOLUTE, LUMP, IHC, CPE")
             }
             sample_names <- extractVialOnly(TCGA_barcodes)
             purity_names <- extractVialOnly(rownames(TCGA_purities))
             name_matching <- match(purity_names, sample_names)
             cut <- TCGA_purities[,method] > threshold
             cut[is.na(cut)] <- F
             cut_idcs <- name_matching[cut]
             cut_idcs <- cut_idcs[!is.na(cut_idcs)]
             return(cut_idcs)
           },

           # return tissue type given an input barcode
           getTissueType = function(TCGA_barcode)
           {
             this_sample = substr(str_split(TCGA_barcode,"-",simplify=T)[1,4],1,2)
             print(this_sample)
             print(as.numeric(this_sample))
             print(sample_type_mapping$numcode)
             print(as.numeric(sample_type_mapping$numcode))
             return(sample_type_mapping[which(as.numeric(sample_type_mapping$numcode) == as.numeric(this_sample)),])
           },

           ## Filtering samples in an rds with a particular sample type (e.g., "Primary Tumor", "Solid Tissue Normal", "Primary Blood Derived Cancer - Peripheral Blood")
           ## 20220920 Man page done
           filterTumorType = function(TCGA_barcodes, type_of_tumor, rds_info){
             if(class(TCGA_barcodes) != "character"){
               stop("Error: TCGA_barcodes argument needs to be a character vector")
             }
             if(class(type_of_tumor) != "character"){
               stop("Error: type_of_tumor argument needs to be a string")
             }
             if(class(rds_info) != "data.frame"){
               stop("Error: expression matrices need to be an RSE object")
             }
             sample_names <- sapply(TCGA_barcodes, substr, 1, 16)
             ## use column names of original object
             type_names <- sapply(rds_info$tcga.tcga_barcode, substr, 1, 16)
             sample_type <- rds_info$tcga.cgc_sample_sample_type == type_of_tumor
             sample_type[is.na(sample_type)] <- F

             return(which(sample_type[match(sample_names, type_names)]))
           },

           ## Filtering all tumor samples (e.g. barcode sample types {01,..09})
           filterTumorSamples = function(TCGA_barcodes){
             if(class(TCGA_barcodes) != "character"){
               stop("Error: TCGA_barcodes argument needs to be a character vector")
             }
             sample_types <- extractSampleType(TCGA_barcodes)
             ## tumor samples have codes between 01 and 09
             tumors = c('01','02','03','04','05','06','07','08','09')
             sample_tumor <- sample_types %in% tumors
             sample_tumor[is.na(sample_tumor)] <- F

             return(which(sample_tumor))
           },

           ## Filtering all normal samples (e.g. barcode sample types {10,..19})
           filterNormalSamples = function(TCGA_barcodes){
             if(class(TCGA_barcodes) != "character"){
               stop("Error: TCGA_barcodes argument needs to be a character vector")
             }
             sample_types <- extractSampleType(TCGA_barcodes)
             ## tumor samples have codes between 01 and 09
             normals = c('10','11','12','13','14','15','16','17','18','19')
             sample_normal <- sample_types %in% normals
             sample_normal[is.na(sample_normal)] <- F

             return(which(sample_normal))
           },

           ## Filtering all control samples (e.g. barcode sample types {20,..29})
           filterControlSamples = function(TCGA_barcodes){
             if(class(TCGA_barcodes) != "character"){
               stop("Error: TCGA_barcodes argument needs to be a character vector")
             }
             sample_types <- extractSampleType(TCGA_barcodes)
             ## tumor samples have codes between 01 and 09
             controls = c('20','21','22','23','24','25','26','27','28','29')
             sample_control <- sample_types %in% controls
             sample_control[is.na(sample_control)] <- F

             return(which(sample_control))
           },

           ## Filter out protein coding genes based on rds info
           filterGenesProteins = function(rds_gene_info){
             if(class(rds_gene_info) != "data.frame"){
               stop("Error: gene info argument should be a data.frame. Best \
                    use extractSampleAndGeneInfo function to retrieve this information \
                    from an rds expression object.")
             }
             return(which(rds_gene_info$gene_type == "protein_coding"))
           },

           ## Filter all genes which have at least *norm_threshold* scores (of normalized gene expression) in at least *sample_fraction* of samples
           ## expression_matrix should be TPM or CPM values (NOT log scaled)
           ## sample_fraction should be in [0,1]
           filterGenesByNormExpression = function(expression_matrix, norm_threshold, sample_fraction){
             if(class(expression_matrix) != "data.frame"){
               stop("Error: expression_matrix argument should be a data.frame")
             }
             if(class(norm_threshold) != "numeric" | norm_threshold <= 0){
               stop("Error: norm_threshold argument should be a numeric > 0")
             }
             if(class(sample_fraction) != "numeric" | sample_fraction <= 0 | sample_fraction >= 1){
               stop("Error: sample_fraction argument should be a numeric in [0,1]")
             }
             minSamples = sample_fraction*ncol(expression_matrix)
             keep = rowSums(expression_matrix >= norm_threshold) >= minSamples
             return(which(keep))
           },

           ## 20220921 man page done
           filterChromosome = function(rds_gene_info, chroms){
             if(class(rds_gene_info) != "data.frame"){
               stop("Error: gene info argument should be a data.frame. Best \
                    use extractSampleAndGeneInfo function to retrieve this information \
                    from an rds expression object.")
             }
             if (class(chroms) != "character"){
               stop("Error: expected chroms argument to be a vector of strings.")
             }
             return(which(rds_gene_info$seqnames %in% chroms))
           },

           ##gets gene information from gencode given a list of genes names or ids
           ##it is supposed to automatically infer wheter id or name
           ##it is supposed to automatically infer whether . exists in id
           ## 20220921 man page done
           getGeneInfo = function(gene_names_or_ids){
             is_id <-  grepl("ENSG", gene_names_or_ids, fixed=TRUE)
             if(any(is_id == TRUE)){
               version <- grepl(".", gene_names_or_ids, fixed=TRUE)
               if(any(version == TRUE)){
                 to_return <- subset(gene_mapping, gene_names_or_ids %in% gene_mapping$gene_id)
               }
               else{
                 to_return <- subset(gene_mapping,  gene_names_or_ids %in% gene_mapping$gene_id_no_ver)
               }
             }
             else{
               to_return <- subset(gene_mapping, gene_names_or_ids %in% gene_mapping$gene_name)
             }
             return(to_return)
           },

           ## the version corresponds to whether we want the . and number after from gene ids
           ## 20220921 man page done
           geneNameToENSG = function(gene_names, version = FALSE){
             to_return <- getGeneInfo(gene_names)
             if(version == TRUE){
               to_return <- to_return$gene_id
             }
             else{
               to_return <- to_return$gene_id_no_ver
             }
             return(to_return)
           },

           geneENSGToName = function(gene_ids){
             to_return <- getGeneInfo(gene_ids)
             return(to_return$gene_name)
           },

           # getGeneAliases(gene_names)
           # # TODO: return all alias names
           # #

           getGeneIdcs = function(gene_names, rds_gene_info){
             if(class(rds_gene_info) != "data.frame"){
               stop("Error: gene info argument should be a data.frame. Best \
                    use extractSampleAndGeneInfo function to retrieve this information \
                    from an rds expression object.")
             }
             if (class(gene_names) != "character"){
               stop("Error: expected gene_names argument to be a vector of strings.")
             }
             return(match(gene_names, rds_gene_info$gene_name))
           },

           getStage = function(TCGA_barcodes){
             sample_names <- extractSampleOnly(TCGA_barcodes)
             stage_names <- clinical_patient_data$bcr_patient_barcode
             stages <- clinical_patient_data$ajcc_pathologic_tumor_stage[match(sample_names, stage_names)]
             return(stages)
           },

           getSex = function(TCGA_barcodes){
             sample_names <- extractSampleOnly(TCGA_barcodes)
             sex_names <- clinical_patient_data$bcr_patient_barcode
             sex <- clinical_patient_data$gender[match(sample_names, sex_names)]
             return(sex)
           }
         )
)

### constructors for NetSciDataCompanion class
### like preparing and creating your object before you can use the methods above
### the export decorator is for roxygen to know which methods to export

#' @export "CreateNetSciDataCompanionObject"
CreateNetSciDataCompanionObject <- function(clinical_patient_file=NULL, project_name="default_project"){

  ## Load purities for purity package
  obj <- CreateTCGAPurityFilteringObject()

  #this is an easy hack for not breaking, but something smarter would be great
  #TODO skip purityfiltering completely and do it here instead
  with_purity = c("ACC","BLCA","BRCA","CESC","COAD","GBM",
                  "HNSC","KIRC","KIRP","KICH","LGG","LIHC",
                  "LUAD","LUSC","OV","PRAD","READ","SKCM",
                  "THCA","UCEC","UCS")
  purities <- data.frame()

  if(project_name %in%  with_purity){
    purities <- obj$get_tissue_purities(cancer_type = project_name)
  }

  ## Load patient's clinical data
  if(!is.null(clinical_patient_file))
  {
    patient_data <- read.table(clinical_patient_file, header=T, sep=",")
    ## maybe we want to keep the alternative column names later? For now this is discarded
    alt_colnames <- patient_data[1:2,]
    patient_data <- patient_data[-c(1,2),]
  }

  else
  {
    patient_data = data.frame()
  }

  fpath <- system.file("extdata", "gen_v26_mapping.csv", package="NetSciDataCompanion")
  gene_mapping <- read.csv(file = fpath, sep=",", header=TRUE, row.names = 1)
  gene_mapping$gene_id_no_ver <- gsub("\\..*","",gene_mapping[,"gene_id"])

  fpath_sample <- system.file("extdata", "TCGA_sample_type.csv", package="NetSciDataCompanion")
  sample_type_mapping <- read.csv(file = fpath_sample, header=T, sep=",")


  s <- NetSciDataCompanion$new(TCGA_purities = purities,
                               clinical_patient_data = patient_data,
                               project_name = project_name,
                               gene_mapping = gene_mapping,
                               sample_type_mapping = sample_type_mapping)
}





