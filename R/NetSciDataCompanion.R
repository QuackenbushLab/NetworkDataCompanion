NetSciDataCompanion=setRefClass("NetSciDataCompanion",
         fields = list(TCGA_purities= "data.frame"),
         methods = list(
           mapBarcodeToBarcode = function(bc1, bc2){
             ### SOME example ASSERTS
             if(class(bc1) != "data.frame" | class(bc2) != "data.frame"){
               stop("Error: barcodes need to be data.frames")
             }


             return("return me baby")
           },
           
           logTPMNormalization = function(expression_matrix){
             ### SOME example ASSERTS
             if(class(expression_matrix) != "data.frame"){
               stop("Error: expression matrices need to be data.frames")
             }
             return("log(TPM(expression_matrix)+1)")
           },
           
           #### more methods go here
           
           # maybe have this presaved in class
           extractSampleOnly = function(TCGA_barcodes){
              return("barcode up to sample")
           },
           
           extractVialOnly = function(TCGA_barcodes){
              return("barcode up to vial")
           },
           
           findDuplicates = function(TCGA_barcodes){
              return("positions with same vial")
           },
           
           mapUUIDtoTCGA = function(UUID){
              return("internal lookup of UUID")
           },
           
           mapProbesToGenes = function(probelist, rangeUp, rangeDown){
              return("matrix of probes to gene name")
           },
           
           convertBetaToM = function(methylation_betas){
              return("M")
           },
           
           filterDuplicatesSeqDepth = function(expression_matrix, TCGA_barcodes){
           ## internally: look for same vials
              return("indices of max seqdepth columns")
           },
           
           filterPurity = function(TCGA_barcodes, method, threshold){
              return("filtered indices")
           },
           
           filterTumorType = function(type_of_tumor, clinical_data){
              return("filtered indices")
           },
           
           filterGenesProteins = function(gene_expression, gene_info){
              return("filteredGenesProteins")
           },
           
           filterGenesByTPM = function(gene_expression_tpm, tpm_threshold, sample_fraction){
             return("filteredGenesByTPM")
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
           
           getGeneIdcs = function(gene_names, data){
            return("returns indices of gene names in gene list")
           }
         )
)

### constructors for NetSciDataCompanion class
### like preparing and creating your object before you can use the methods above
### the export decorator is for roxygen to know which methods to export

#' @export "CreateNetSciDataCompanionObject"
CreateNetSciDataCompanionObject <- function(argument_one = "An argument"){
  
  ##if external data are to be loaded, they should be placed in the inst/extdata folder
  ##for example, these could be the gene2gene2gene table computed already
  fpath <- system.file("extdata", "purities.csv", package="NetSciDataCompanion")
  
  external_data <- read.csv(file = fpath, sep="\t", header=TRUE, row.names = 1)
  
  
  ### do a bunch of stuff and arrive at the variables to be the fields of the object
  TCGA_purities <- data.frame()
  print(paste("I just created an object with argument: ",argument_one))
  s <- NetSciDataCompanion$new(TCGA_purities = TCGA_purities)
}





