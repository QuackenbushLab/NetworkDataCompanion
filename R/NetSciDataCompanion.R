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
             return(list(rds_sample_info=as.data.frame(colData(expression_rds_obj)),
                         rds_gene_info=as.data.frame(rowRanges(expression_rds_obj))))
           },


           ## Maps two sets of barcodes
           ## There are 4 different pieces of information returned in a named list that are all useful depending on the context they are used in
           ## is_inter1 is an indicator (boolean) vector of the same length as bc1 that indicates which elements of bc1 are present in bc2
           ## idcs1 indicate where to find each barcode of bc1 in bc2, NA if missing. That is, bc1[i] == bc2[idcs1[i]] (if idcs1[i] != NA)
           ## The same information is provided for bc2
           ## For example, if you want to map experiment 1 on experiment two, keeping only the information for samples that are present in both,
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
             if(!("data.frame" %in% class(exp1) | "matrix" %in% class(exp1)) ){
               stop("Error: argument 1 needs to be data.frame or matrix")
             }
             if(!("data.frame" %in% class(exp2) | "matrix" %in% class(exp2)) ){
               stop("Error: argument 2 needs to be data.frame or matrix")
             }
             map <- mapBarcodeToBarcode(extractSampleOnly(colnames(exp1)), extractSampleOnly(colnames(exp2)))
             return(list("mappedExp1"=exp1[,map$is_inter1], "mappedExp2"=exp2[,map$idcs1]))
           },


           ## Computes the log TPM normalization based on an expression RDS object
           ## Returns a named list with the count data.frame (useful for duplicate filtering based on sequencing depth, see filterDuplicatesSeqDepth)
           ##                               TPM data.frame (useful for TPM based filtering, see filterGenesByNormExpression)
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

           # Computes the log transformed CPM normalization based on raw expression (count) data
           # This is following the recipe provided by edgeR package to get TMM valyes
           ## Returns a named list with the count data.frame (useful for duplicate filtering based on sequencing depth, see filterDuplicatesSeqDepth)
           ##                               CPM data.frame (useful for CPM based filtering, see filterGenesByCPM)
           ##                and the actual logCPM which corresponds to log(CPM + 1)
           logCPMNormalization = function(exp_count_mat){
             if(sum(class(exp_count_mat) %in% c("data.frame", "matrix")) == 0){
               stop("Error: expression matrices need to be an RSE object")
             }
             # apply edgeR function to get differential gene lists
             dge <- DGEList(exp_count_mat)
             # get the normalizing factors from that list
             dge <- calcNormFactors(dge)
             # get the log tranformed cpms (aka: get TMM)
             cpm_vals <- cpm(dge, log=F)
             return(list(counts=exp_count_mat,
                         CPM=cpm_vals,
                         logCPM=log(cpm_vals+1)))
           },

           #### more methods go here

           # maybe have this presaved in class
           extractSampleOnly = function(TCGA_barcodes){
             return(sapply(TCGA_barcodes, substr, 1, 12))
           },

           extractSampleAndType = function(TCGA_barcodes){
             return(sapply(TCGA_barcodes, substr, 1, 15))
           },

           extractVialOnly = function(TCGA_barcodes){
              return(sapply(TCGA_barcodes, substr, 1, 16))
           },


           extractSampleType = function(TCGA_barcodes){
              return(sapply(TCGA_barcodes, substr, 14, 15))
           },

           findDuplicates = function(TCGA_barcodes){
             dupPos = duplicated(extractVialOnly(TCGA_barcodes))
             return(TCGA_barcodes[dupPos])
           },

           mapUUIDtoTCGA = function(UUID){
              if(class(UUID) != "character"){
                stop("Error: Expected UUID argument to be vector of strings")
              }
              info = files() %>%
               GenomicDataCommons::filter( ~ file_id %in% UUID) %>%
               GenomicDataCommons::select('cases.samples.submitter_id') %>%
               results_all()
              if(length(info)==0)
              {
                stop("Error: No UUIDs were matched by TCGA ids. Perhaps you have input legacy UUIDs?")

              }
              # The mess of code below is to extract TCGA barcodes
              # id_list will contain a list (one item for each file_id)
              # of TCGA barcodes of the form 'TCGA-XX-YYYY-ZZZ'
              id_list = lapply(info$cases,function(a) {
               a[[1]][[1]][[1]]})
              # so we can later expand to a data.frame of the right size
              barcodes_per_file = sapply(id_list,length)
              # sort to match input UUID order
              file_id <- rep(GenomicDataCommons::ids(info),barcodes_per_file)
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
           # default if probelist is NULL is to map all the probes in the manifest
           mapProbesToGenes = function(probelist = NULL,
                                       rangeUp = 200,
                                       rangeDown = 0,
                                       localManifestPath=NA,
                                       mapToNearest = F){

             if(is.na(localManifestPath))
             {
               print("[NetSciDataCompanion::mapProbesToGenes] Sourcing 450k probe annotation from https://zwdzwd.github.io/InfiniumAnnotation ...")
               print("[NetSciDataCompanion::mapProbesToGenes] https://zhouserver.research.chop.edu/InfiniumAnnotation/20210615/HM450/HM450.hg38.manifest.gencode.v36.tsv.gz")
               print("[NetSciDataCompanion::mapProbesToGenes] HG38, Gencode v36")

               # source hg38 with gencode 36 from https://zwdzwd.github.io/InfiniumAnnotation
               download.file('https://zhouserver.research.chop.edu/InfiniumAnnotation/20210615/HM450/HM450.hg38.manifest.gencode.v36.tsv.gz',
                             destfile = "./HM450.hg38.manifest.gencode.v36.tsv.gz")

               # unzip
               # system2(command="gunzip",args=c("./HM450.hg38.manifest.gencode.v36.tsv.gz"))
               R.utils::gunzip("./HM450.hg38.manifest.gencode.v36.tsv.gz", overwrite=T)

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

             if(is.null(probelist))
               probelist = manifest$probeID
             # get indices matching probes
             smallManifest = manifest %>% dplyr::filter(probeID %in% probelist)
             rm(manifest)
             gc()

             # define empty map
             mymap = matrix(rep(NA,4*nrow(smallManifest)),ncol=4)
             colnames(mymap) = c("probeID","geneNames","ensemblID","distToTSS")

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


             if(!mapToNearest)
              return(mymap)
             if(mapToNearest)
             {
               processRow = function(x) # x is one row of mymap
               {
                 x = as.data.frame(x)
                 genes = str_split(x$geneName,";",simplify=T)
                 ensemblID = str_split(x$ensemblID,";",simplify=T)
                 tssDist = as.numeric(str_split(x$distToTSS,";",simplify=T))
                 if(length(unique(genes)) > 1)
                 {
                   index = which.min(abs(tssDist))
                   thisRow = c(x$probeID,genes[index],ensemblID[index],tssDist[index])
                   return(thisRow)
                 }

                 return(x[1,])
               }

               myNearMap = mymap
               for(i in 1:nrow(mymap))
               {
                 if(i %% 10000 == 0) print(paste("probe number:",i))
                 myNearMap[i,]=processRow(mymap[i,])
               }
               return(myNearMap)
             }
           },

           # Function to map to probes to a gene-level measurement
           # probe_gene_map is in the format output from the mapProbesToGenes function
           # not all genesOfInterest need to be in probe_gene_map, but if none are, then this is meaningless
           probeToMeanPromoterMethylation = function(methylation_betas, probe_gene_map, genesOfInterest){

             # merge probe_gene_map with beta values
             # use a left join to keep only probes that mapped to genes of interest

             mappedBetas = probe_gene_map %>%
               dplyr::filter(geneNames %in% genesOfInterest) %>%
               dplyr::select(geneNames,probeID) %>%
               left_join(methylation_betas, by="probeID")

             ## split the probe map to a long form where there are multiple genes mapped to the same probe
             mappedBetasLong = mappedBetas %>%
               separate_rows(geneNames, sep = ";") %>%
               drop_na(geneNames) %>%
               data.frame(check.names = F)

             ## map probe-level methylation to the mean for each gene
             betaMeans = mappedBetasLong %>%
               group_by(geneNames) %>%
               summarise_at(colnames(mappedBetasLong)[3:(ncol(mappedBetasLong))], mean) %>%
               data.frame(check.names = F, row.names=1) %>% t()
             return(betaMeans)
           },

           # Input to convertBetaToM is a vector of methylation betas
           # User should use this function with `apply` to convert a matrix
           # 20220920 man page done
           convertBetaToM = function(methylation_betas){
              M = log2(methylation_betas/(1-methylation_betas))
              return(M)
           },

           ## Run EPISCORE to estimate cell counts
           estimateCellCountsEpiSCORE = function(methylation_betas, tissue, array = "450k"){
             tissue_options = c("Bladder",
                                "Brain",
                                "Breast",
                                "Colon",
                                "Heart",
                                "Kidney",
                                "Liver",
                                "Lung",
                                "OE", # olfactory epithelial
                                "Pancreas_6ct", # defined over 6 cell types
                                "Pancreas_9ct", # defined over 9 cell types
                                "Prostate",
                                "Skin")
             if(!tissue %in% tissue_options)
             {
               print(paste("[NetSciDataCompanion::estimateCellCountsEpiSCORE()] EpiSCORE is not implemented for tissue:", tissue))
             }

             # map methylation to gene level
             # methylation betas should be samples in rows
             # and CGs in columns
             # luad for testing
             # methylation_betas = fread("~/Desktop/tcga_luad_methylations.txt",data.table=F)
             # row.names(methylation_betas) = methylation_betas$probeID
             # array = "450k"

	     row.names(methylation_betas) = methylation_betas$probeID

             geneLevelMeth = methylation_betas %>%
               dplyr::select(-probeID) %>%
               as.matrix() %>%
               constAvBetaTSS(type = array)

             cellEst = NULL

             if(tissue == "Bladder"){
               #data(BladderRef)
               cellEst = wRPC(data = geneLevelMeth, ref.m = mrefBladder.m)
             }
             if(tissue == "Brain"){
               #data(BrainRef)
               cellEst = wRPC(data = geneLevelMeth, ref.m = mrefBrain.m)
             }
             if(tissue == "Breast"){
               #data(BreastRef)
               cellEst = wRPC(data = geneLevelMeth, ref.m = mrefBreast.m)
             }
             if(tissue == "Colon")
             {
               #data(ColonRef)
               cellEst = wRPC(data = geneLevelMeth, ref.m = Colon_Mref.m)
             }
             if(tissue == "Heart"){
               #data(HeartRef)
               cellEst = wRPC(data = geneLevelMeth, ref.m = mrefHeart.m)
             }
             if(tissue == "Kidney"){
               #data(KidneyRef)
               cellEst = wRPC(data = geneLevelMeth, ref.m = Kidney_Mref.m)
             }
             if(tissue == "Liver"){
               #data(LiverRef)
               cellEst = wRPC(data = geneLevelMeth, ref.m = mrefLiver.m)
             }
             if(tissue == "Lung")
             {
               #data(LungRef)
               cellEst = wRPC(data = geneLevelMeth, ref.m = mrefLung.m)
             }
             if(tissue == "OE"){
               #data(OEref)
               cellEst = wRPC(data = geneLevelMeth, ref.m = mrefOE.m)
             }
             if(tissue == "Pancreas_6ct"){
               cellEst = wRPC(data = geneLevelMeth, ref.m = mrefPancreas.m)
             }
             if(tissue == "Pancreas_9ct"){
               cellEst = wRPC(data = geneLevelMeth, ref.m = mrefPancreas9ct.m)
             }
             if(tissue == "Prostate"){
               cellEst = wRPC(data = geneLevelMeth, ref.m = mrefProstate.m)
             }
             if(tissue == "Skin"){
               cellEst = wRPC(data = geneLevelMeth, ref.m = mrefSkin.m)
             }


	          # from cellEst, take estF
	          outData = data.frame(cellEst$estF)
	          outData$UUID = row.names(cellEst$estF)

	          # get TCGA barcodes
	          outTCGA = mapUUIDtoTCGA(outData$UUID)

	          # merge and relabel
	          outLabeled = outTCGA %>% inner_join(outData,by=c("file_id"="UUID")) %>%
	            dplyr::rename("TCGAbarcode"="submitter_id","UUID"="file_id")

            return(outLabeled)
           },

           ## Extract AHRR methylation at probe site cg05575921 as a proxy for smoking status
           extractAHRRMethylation = function(methylation_betas)
           {
             ahrr = methylation_betas %>%
                    dplyr::filter(probeID == "cg05575921") %>%
                    dplyr::select(-probeID) %>%
                    t()
             colnames(ahrr)[1]="ahrr_cg05575921_beta"
             return(as.data.frame(ahrr))
           },

           ## Filter out all duplicates based on sequencing depth
           ## Returns indices about which samples to KEEP
           ## 20220920 man page done
           filterDuplicatesSeqDepth = function(expression_count_matrix){
             sample_barcodes <- extractSampleAndType(colnames(expression_count_matrix))
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
             this_sample = as.numeric(substr(str_split(TCGA_barcode,"-",simplify=T)[1,4],1,2))
             if(!this_sample%in% sample_type_mapping$numcode)
             {
               print(paste("[NetSciDataCompanion::getTissueType()] Error: unknown sample type:",this_sample))
               return(NA)
             }

             this_map = data.frame("TCGA_barcode"=TCGA_barcode)
             row.names(this_map) = TCGA_barcode
             this_match = sample_type_mapping[which(as.numeric(sample_type_mapping$numcode) == as.numeric(this_sample)),]

             return(cbind.data.frame(this_map,this_match))
           },
           ## access sample type map
	         getSampleTypeMap = function(){
	           return(sample_type_mapping)
	         },
           ## Filtering samples in an rds with particular sample types (e.g., "Primary Tumor", "Solid Tissue Normal", "Primary Blood Derived Cancer - Peripheral Blood")
           filterSampleType = function(TCGA_barcodes, types_of_samples){
             if(class(TCGA_barcodes) != "character"){
               stop("Error: TCGA_barcodes argument needs to be a character vector")
             }
             if(class(types_of_samples) != "character"){
               stop("Error: types_of_sample argument needs to be a character vector. Use NetSciDataCompanion::getSampleTypeMap() to see available types.")
             }

             observed_sample_types = extractSampleType(TCGA_barcodes)
             nonExistTypes <- which(!(types_of_samples %in% observed_sample_types))
             if (length(nonExistTypes) > 0) {

               if (length(nonExistTypes) == length(types_of_samples)){
                  stop("Error: No specified types in types_of_sample argument exist in sample info.\n Use NetSciDataCompanion::getSampleTypeMap() to see available types.")
               }
               print(paste0("Warning: sample types ", types_of_samples[nonExistTypes], " are not present in sample info."))
             }

             return(data.frame("sample_type"=observed_sample_types[which(observed_sample_types %in% types_of_samples)],
                               "index"=which(observed_sample_types %in% types_of_samples)))

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
             ## normal samples have codes between 10 and 19
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
             ## control samples have codes between 20 and 29
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
             if(sum(class(expression_matrix) %in% c("data.frame","matrix")) == 0) {
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

    	     ##gets gene information from gencode given a list of genes names or ids from the gene mapping variable
    	     ##automatically infers whether ensmbl ID or name or entrez
    	     ##automatically infers whether version (i.e., the dot) exists in ensemble ID
    	     getGeneInfo = function(gene_names_or_ids_or_entrezs){
    	       is_id <-  any(grepl("ENSG", gene_names_or_ids_or_entrezs, fixed=TRUE))
    	       is_entrez <- any(grepl("^\\d+$", gene_names_or_ids_or_entrezs))

    	       if(is_id){
    	         version <- grepl(".", gene_names_or_ids_or_entrezs, fixed=TRUE)
    	         if(any(version == TRUE)){
    	           to_return <- subset(gene_mapping, gene_mapping$gene_id %in% gene_names_or_ids_or_entrezs )
    	         }
    	         else{
    	           to_return <- subset(gene_mapping, gene_mapping$gene_id_no_ver %in% gene_names_or_ids_or_entrezs)
    	         }
    	       }
    	       else if (is_entrez){
    	         to_return <- subset(gene_mapping, gene_mapping$gene_entrez %in% gene_names_or_ids_or_entrezs)
    	       }
    	       else{
    	         to_return <- subset(gene_mapping, gene_mapping$gene_name %in% gene_names_or_ids_or_entrezs)
    	       }

    	       if(nrow(to_return)!=length(gene_names_or_ids_or_entrezs)){
    	         print('There was at least one one-to-many mapping (most probably from multiple ensembl IDs for the input)')
    	       }
    	       return(to_return)
    	     },

    	     ## the version corresponds to whether we want the . and number after from gene ids
    	     geneEntrezToENSG = function(gene_entrezs, version = FALSE){
    	       if(!("gene_entrez" %in% colnames(gene_mapping)))
    	       {
    	         stop('Column gene_entrez not found in gene mapping.')
    	       }
    	       to_return <- getGeneInfo(gene_names)
    	       if(version == TRUE){
    	         to_return <- to_return[c('gene_entrez','gene_id')]
    	       }
    	       else{
    	         to_return <- to_return[c('gene_entrez','gene_id_no_ver')]
    	       }
    	       return(to_return)
    	     },

    	     geneENSGToName = function(gene_ids){
    	       to_return <- getGeneInfo(gene_ids)
    	       if(anyNA(to_return$gene_name)){
    	         print('Not all ensembl IDs were mapped to names')
    	       }
    	       return(to_return[c('gene_id_no_ver','gene_name')])
    	     },

    	     geneENSGToEntrez = function(gene_ids){
    	       if(!("gene_entrez" %in% colnames(gene_mapping)))
    	       {
    	         stop('Column gene_entrez not found in gene mapping.')
    	       }
    	       to_return <- getGeneInfo(gene_ids)
    	       if(anyNA(to_return$gene_entrez)){
    	         print('Not all ensembl IDs were mapped to entrez')
    	       }
    	       return(to_return[c('gene_id_no_ver','gene_entrez')])
    	     },

    	     geneNameToEntrez = function(gene_names){
    	       if(!("gene_entrez" %in% colnames(gene_mapping)))
    	       {
    	         stop('Column gene_entrez not found in gene mapping')
    	       }
    	       to_return <- getGeneInfo(gene_names)
    	       if(anyNA(to_return$gene_entrez)){
    	         print('Not all names were mapped to entrez')
    	       }
    	       to_return <- to_return[c('gene_name','gene_entrez')]
    	       to_return <- to_return[!duplicated(to_return),]
    	       return(to_return)
    	     },

    	     geneEntrezToName = function(gene_entrezs){
    	       if(!("gene_entrez" %in% colnames(gene_mapping)))
    	       {
    	         stop('Column gene_entrez not found in gene mapping')
    	       }
    	       to_return <- getGeneInfo(gene_entrezs)
    	       if(anyNA(to_return$gene_name)){
    	         print('Not all entrez were mapped to names')
    	       }
    	       to_return <- to_return[c('gene_entrez','gene_name')]
    	       to_return <- to_return[!duplicated(to_return),]
    	       return(to_return)
    	     },

    	     ## the version corresponds to whether we want the . and number after from gene ids
    	     geneNameToENSG = function(gene_names, version = FALSE){
    	       to_return <- getGeneInfo(gene_names)
    	       if(version == TRUE){
    	         to_return <- to_return[c('gene_name','gene_id')]
    	       }
    	       else{
    	         to_return <- to_return[c('gene_name','gene_id_no_ver')]
    	       }
    	       return(to_return)
    	     },

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

  fpath <- system.file("extdata", "gen_v26_mapping_w_entrez.csv", package="NetSciDataCompanion")
  gene_mapping <- read.csv(file = fpath, sep=",", header=TRUE, row.names = 1)
  gene_mapping$gene_id_no_ver <- gsub("\\..*","",gene_mapping[,"gene_id"])

  fpath_sample <- system.file("extdata", "TCGA_sample_type.csv", package="NetSciDataCompanion")
  sample_type_mapping <- read.csv(file = fpath_sample, header=T, sep=",",
                                  colClasses = "character") # read codes as characters so that 01, 02, etc. are read properly

  s <- NetSciDataCompanion$new(TCGA_purities = purities,
                               clinical_patient_data = patient_data,
                               project_name = project_name,
                               gene_mapping = gene_mapping,
                               sample_type_mapping = sample_type_mapping)
}





