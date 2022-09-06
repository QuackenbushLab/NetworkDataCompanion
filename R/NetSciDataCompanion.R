NetSciDataCompanion=setRefClass("NetSciDataCompanion",
         fields = list(TCGA_purities= "data.frame",
                       clinical_patient_data = "data.frame",
                       project_name = "character"),
         methods = list(


           ## Extract experiment specific information and metadata from ranged summarized experiment object
           ## Returns a named list with rds_sample_info corresponding to meta information about the samples (columns)
           ##                       and rds_gene_info corresponding to meta information about genes (rows)
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
           ## exp2[,idcs1[!is.na(idcs1)]]     --- this will remove samples that are not in exp1 and reorder to match exp1
           mapBarcodeToBarcode = function(bc1, bc2){
             if(class(bc1) != "character" | class(bc2) != "character"){
               stop("Error: barcodes need to be vectors of strings")
             }
             m1 <- match(bc1, bc2)
             m1 <- m1[!is.na(m1)]
             m2 <- match(bc2, bc1)
             m2 <- m2[!is.na(m2)]
             return(list(is_inter1=(bc1 %in% bc2), idcs1=m1, isinter2=(bc2 %in% bc1), idcs2=m2))
           },


           ## Computes the log TPM normalization based on an expression RDS object
           ## Returns a named list with the count data.frame (useful for duplicate filtering based on sequencing depth, see filterDuplicatesSeqDepth)
           ##                               TPM data.frame (useful for TPM based filtering, see filterGenesByTPM)
           ##                and the actual logTPM which corresponds to log(TPM + 1)
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

           #### more methods go here

           # maybe have this presaved in class
           extractSampleOnly = function(TCGA_barcodes){
             return(sapply(TCGA_barcodes, substr, 1, 12))
           },

           extractVialOnly = function(TCGA_barcodes){
              return(sapply(TCGA_barcodes, substr, 1, 16))
           },

           findDuplicates = function(TCGA_barcodes){
              return(duplicated(extractVialOnly(TCGA_barcodes)))
           },

           mapUUIDtoTCGA = function(UUID){
              if(class(UUID) != "character"){
                stop("Error: Expected UUID argument to be vector of strings")
              }
              info = files(legacy = T) %>%
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
              # And build the data.frame
              return(data.frame(file_id = rep(ids(info),barcodes_per_file),
                               submitter_id = unlist(id_list)))
           },

           mapProbesToGenes = function(probelist, rangeUp, rangeDown){
              return("matrix of probes to gene name")
           },

           # Input to convertBetaToM is a vector of methylation betas
           # User should use this function with `apply` to convert a matrix
           convertBetaToM = function(methylation_betas){
              M = log2(methylation_betas/(1-methylation_betas))
              return(M)
           },

           ## Filter out all duplicates based on sequencing depth
           ## Returns indices about which samples to KEEP
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

           ## Filtering samples with a particular sample type (e.g., "Primary Tumor", "Solid Tissue Normal", "Primary Blood Derived Cancer - Peripheral Blood")
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

           ## Filter out protein coding genes based on rds info
           filterGenesProteins = function(rds_gene_info){
             if(class(rds_gene_info) != "data.frame"){
               stop("Error: gene info argument should be a data.frame. Best \
                    use extractSampleAndGeneInfo function to retrieve this information \
                    from an rds expression object.")
             }
             return(which(rds_gene_info$gene_type == "protein_coding"))
           },

           ## Filter all genes which have at least *tpm_threshold* TPM scores in at least *sample_fraction* of samples
           ## expression_tpm_matrix should be TPM values (NOT log scaled)
           ## sample_fraction should be in [0,1]
           filterGenesByTPM = function(expression_tpm_matrix, tpm_threshold, sample_fraction){
             if(class(expression_tpm_matrix) != "data.frame"){
               stop("Error: expression_tpm_matrix argument should be a data.frame")
             }
             if(class(tpm_threshold) != "numeric" | tpm_threshold <= 0){
               stop("Error: tpm_threshold argument should be a numeric > 0")
             }
             if(class(sample_fraction) != "numeric" | sample_fraction <= 0 | sample_fraction >= 1){
               stop("Error: sample_fraction argument should be a numeric in [0,1]")
             }
             minSamples = sample_fraction*ncol(expression_tpm_matrix)
             keep = rowSums(expression_tpm_matrix >= tpm_threshold) >= minSamples
             return(which(keep))
           },

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

           ######the following 3 could be implemented in Gene2Gene2Gene and used directly from there

           # geneNameToENSG(gene_names)
           # # use Panos genetogenetogene mapping
           #
           # geneENSGToName(gene_names)
           # # use Panos genetogenetogene mapping
           #
           # getGeneAliases(gene_names)
           # # return all alias names
           # # use Panos genetogenetogene mapping

           getGeneIdcs = function(gene_names, rds_gene_info){
             if(class(rds_gene_info) != "data.frame"){
               stop("Error: gene info argument should be a data.frame. Best \
                    use extractSampleAndGeneInfo function to retrieve this information \
                    from an rds expression object.")
             }
             if (class(gene_names) != "character"){
               stop("Error: expected gene_names argument to be a vector of strings.")
             }
             return(which(gene_names) %in% rds_gene_info)
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
CreateNetSciDataCompanionObject <- function(clinical_patient_file, project_name){

  ## Load purities for purity package
  obj <- CreateTCGAPurityFilteringObject()
  purities <- obj$get_tissue_purities(cancer_type = project_name)

  ## Load patient's clinical data
  patient_data <- read.table(clinical_patient_file, header=T, sep=",")
  ## maybe we want to keep the alternative column names later? For now this is discarded
  alt_colnames <- patient_data[1:2,]
  patient_data <- patient_data[-c(1,2),]

  s <- NetSciDataCompanion$new(TCGA_purities = purities,
                               clinical_patient_data = patient_data,
                               project_name = project_name)
}





