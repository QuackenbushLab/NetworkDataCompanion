library(NetworkDataCompanion)
# library(TCGAPurityFiltering)
# library(recount)
# library(recount3)
# library(GenomicDataCommons)
# library(magrittr)
###detach("package:NetworkDataCompanion")

datapath = "testdata/luad/"

project_name <- "LUAD"
patient_data <- paste0(datapath,"clinical_patient_luad.csv")
exp_data <- paste0(datapath,"tcga_luad.rds")
mut_data <- paste0(datapath,"tcga_luad_mutations.txt")
mut_pivot_data <- paste0(datapath, "tcga_luad_mutations_pivot.csv")
meth_data <- paste0(datapath, "tcga_luad_methylations.txt")


# Create Companion Object
obj <- CreateNetworkDataCompanionObject(clinical_patient_file = patient_data,
                                       project_name = project_name)

# Read RDS expression data
test_exp_rds <- readRDS(exp_data)

# Extract included meta info from RDS expression object
rds_info <- obj$extractSampleAndGeneInfo(test_exp_rds)
rds_sample_info <- rds_info$rds_sample_info
rds_gene_info <- rds_info$rds_gene_info

# Some gene mapping
gene_info_from_names <- obj$getGeneInfo(rds_gene_info$gene_name)
genes_info_from_ids <- obj$getGeneInfo(rds_gene_info$gene_id)

genes_from_names_to_ids <- data.frame(obj$geneNameToENSG(rds_gene_info$gene_name))
genes_from_ids_to_names <- data.frame(obj$geneENSGToName(rds_gene_info$gene_id))

gene_mapping <- obj$gene_mapping
library(dplyr)
to_return <- subset(gene_mapping,   rds_gene_info$gene_name %in%  gene_mapping$gene_name)

# Normalize data
test_exp_all <- obj$logTPMNormalization(test_exp_rds)
test_exp_count <- test_exp_all$counts
test_exp_tpm <- test_exp_all$TPM
test_exp_logtpm <- test_exp_all$logTPM


# Map column names to TCGA barcodes
newcolnames <- obj$mapUUIDtoTCGA(colnames(test_exp_logtpm))
colnames(test_exp_count) <- newcolnames[,2]
colnames(test_exp_tpm) <- newcolnames[,2]
colnames(test_exp_logtpm) <- newcolnames[,2]

# Get indices of samples passing purity filtering
idcs_purity <- obj$filterPurity(colnames(test_exp_count))
# Get indices of samples that are primary tumor
idcs_tumor <- obj$filterTumorType(colnames(test_exp_count), "Primary Tumor", rds_sample_info)
idcs_normal <- obj$filterTumorType(colnames(test_exp_count), "Solid Tissue Normal", rds_sample_info)
# Get indices of nonduplicates
idcs_nonduplicate_tumor <- obj$filterDuplicatesSeqDepth(expression_count_matrix = test_exp_count[,idcs_tumor])
idcs_nonduplicate_normal <- obj$filterDuplicatesSeqDepth(expression_count_matrix = test_exp_count[,idcs_normal])

# Get indices of genes that have a minimum TPM in the data
test_exp_tpm_df <- data.frame(test_exp_tpm)
idcs_genes_mintpm <- obj$filterGenesByTPM(test_exp_tpm_df, 1, 0.5)




# read mutations and pivot
mutations <- read.table(mut_data, header=T, sep=" ")
mutations_pivot <- read.table(mut_pivot_data, header=T, sep=",")
# NOTE: colnames (barcodes) are with . instead of -, replace
colnames(mutations_pivot) <- gsub('.','-',colnames(mutations_pivot), fixed=T)
# filter duplicates based on sequencing depth of expression in corresponding vials
idcs_nonduplicate_muta <- obj$filterDuplicatesSeqDepthOther(expression_count_matrix = test_exp_count[,idcs_tumor],
                                                       obj$extractVialOnly(colnames(mutations_pivot[,-c(1,2)]))) + 2



# read methylation
#methylations <- read.table(meth_data, header=T, sep=" ")


# ASSEMBLE FILTERED DATA

# obtain preprocessed expression data for filtered samples and filtered genes that are also present in the other omics or for which we have normal
tumor_exp <- test_exp_logtpm[idcs_genes_mintpm, idcs_tumor][, intersect(idcs_nonduplicate_tumor, idcs_purity[idcs_tumor])]
normal_exp <- test_exp_logtpm[idcs_genes_mintpm, idcs_normal][, idcs_nonduplicate_normal]
mutations_filtered <- mutations_pivot[, c(1,2,idcs_nonduplicate_muta)]

## match to normal
common_samples <- obj$mapBarcodeToBarcode(obj$extractSampleOnly(colnames(tumor_exp)),
                                          obj$extractSampleOnly(colnames(normal_exp)))
tumor_exp_matchednormal <- tumor_exp[,common_samples$is_inter1]
normal_exp_matchednormal <- normal_exp[,common_samples$idcs1]
## alternatively you can use a convenience wrapper
matched_exps <- obj$filterBarcodesIntersection(tumor_exp, normal_exp)
tumor_exp_matchednormal <- matched_exps[[1]]
normal_exp_matchednormal <- matched_exps[[2]]

# Get for example the sex of the matched samples
obj$getSex(colnames(tumor_exp_matchednormal))



## match to mutations
common_samples <- obj$mapBarcodeToBarcode(obj$extractSampleOnly(colnames(tumor_exp)),
                                          obj$extractSampleOnly(colnames(mutations_filtered[,-c(1,2)])))

tumor_exp_matchedmutations <- tumor_exp[,common_samples$is_inter1]
mutations_exp_matchedmutations <- cbind(mutations_filtered[,c(1,2)], mutations_filtered[,-c(1,2)][,common_samples$idcs1])

## test probe map generation
## this is a test only in the sense that no
## unrecoverable errors happen
## it is not testing correctness of the probe mapping

obj <- CreateNetworkDataCompanionObject()
myProbeList = c("cg14008030","cg12045430","cg03130891")
shortMap = obj$mapProbesToGenes(myProbeList,rangeUp = 1500, rangeDown = 0,localManifestPath = NA)
longMap = obj$mapProbesToGenes(myProbeList,rangeUp = 1500, rangeDown = 0,
                               localManifestPath = NA , longForm = T)

