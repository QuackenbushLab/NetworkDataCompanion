library(NetSciDataCompanion)
library(TCGAPurityFiltering)
library(recount)
library(recount3)
library(GenomicDataCommons)
library(magrittr)
###detach("package:NetSciDataCompanion")


setwd("~/Projects/TCGA_data/preprocessing/")


project_name <- "LUAD"
patient_data <- "~/Projects/TCGA_data/preprocessing/luad/clinical_patient_luad.csv"
exp_data <- "~/Projects/TCGA_data/preprocessing/luad/tcga_luad.rds"
mut_data <- "~/Projects/TCGA_data/preprocessing/luad/tcga_luad_mutations.txt"
mut_pivot_data <- "~/Projects/TCGA_data/preprocessing/luad/tcga_luad_mutations_pivot.csv"
meth_data <- "~/Projects/TCGA_data/preprocessing/luad/tcga_luad_methylations.txt"


# Create Companion Object
## TODO: this requires TCGAPurityFiltering to be loaded in this file
obj <- CreateNetSciDataCompanionObject(clinical_patient_file = patient_data,
                                       project_name = project_name)

# Read RDS expression data
test_exp_rds <- readRDS(exp_data)

# Extract included meta info from RDS expression object
## TODO: this requires recount to be loaded by this file
rds_info <- obj$extractSampleAndGeneInfo(test_exp_rds)
rds_sample_info <- rds_info$rds_sample_info
rds_gene_info <- rds_info$rds_gene_info

# Normalize data
test_exp_all <- obj$logTPMNormalization(test_exp_rds)
test_exp_count <- test_exp_all$counts
test_exp_tpm <- test_exp_all$TPM
test_exp_logtpm <- test_exp_all$logTPM


# Map column names to TCGA barcodes
## TODO: this is not working. Can't update package.
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
idcs_genes_mintpm <- obj$filterGenesByTPM(expression_tpm_matrix = test_exp_tpm)




# read mutations and pivot
mutations <- read.table(mut_data, header=T, sep=" ")
mutations_pivot <- read.table(mut_pivot_data, header=T, sep=",")
# NOTE: colnames are with . instead of -, replace
colnames(mutations_pivot) <- gsub('.','-',colnames(mutations_pivot), fixed=T)
# TODO: filter duplicates.



# read methylation
#methylations <- read.table(meth_data, header=T, sep=" ")


# ASSEMBLE FILTERED DATA

# obtain preprocessed expression data for filtered samples and filtered genes that are also present in the other omics
tumor_exp <- test_exp_logtpm[idcs_genes_mintpm, idcs_tumor][, intersect(idcs_nonduplicate_tumor, idcs_purity[idcs_tumor])]
normal_exp <- test_exp_logtpm[idcs_genes_mintpm, idcs_normal][, intersect(idcs_nonduplicate_normal, idcs_purity[idcs_normal])]
common_samples <- obj$mapBarcodeToBarcode(obj$extractVialOnly(colnames(tumor_exp)),
                                          obj$extractVialOnly(colnames(mutations_pivot)[-c(1,2)]))
tumor_exp[,common_samples$is_inter1][,common_samples$idcs1[common_samples$is_inter1]]
mutations_exp <- cbind(mutations_pivot[,c(1,2)], mutations_pivot[,-c(1,2)][,common_samples$is_inter2])

print(head(tumpr_exp))
print(head(mutations_exp))
