# Test probeToMeanTFMethylation
context("[NetSciDataCompanion] Testing probeToMeanTFMethylation function ... ")

# make toy beta values
my_betas = data.frame("probeID"=c("cg00000001",
                                  "cg00000002",
                                  "cg00000003",
                                  "cg00000004",
                                  "cg00000005"))

my_betas$minerva = c(0.1,0.2,0.3,0.4,0.5)
my_betas$albus = c(0.9,0.8,0.7,0.6,0.5)

test_that("Mean calculation is correct with default mapToNearest = F, including probes mapping to multiple genes and/or isoforms.",{

  my_friend = NetSciDataCompanion::CreateNetSciDataCompanionObject()
  my_map = my_friend$mapProbesToGenes(probelist = c("cg00000001","cg00000002","cg00000003","cg00000004","cg00000005"),
                                      rangeUp = 201,
                                      rangeDown = 0,
                                      localManifestPath = "testdata/manifest_test.tsv")
  names(my_map)[2] = "geneNames"
  mean_meth = my_friend$probeToMeanPromoterMethylation(methylation_betas = my_betas,
                                           genesOfInterest = c("HARRY","SEVERUS"),
                                           probe_gene_map = my_map)
  minerva_mean_harry = mean(my_betas$minerva[1:4])
  minerva_mean_severus = 1/3*(my_betas$minerva[4] + 2*my_betas$minerva[5])
  albus_mean_harry = mean(my_betas$albus[1:4])
  albus_mean_severus = 1/3*(my_betas$albus[4] + 2*my_betas$albus[5])

  manual_mean = matrix(c(minerva_mean_harry, minerva_mean_severus,
                        albus_mean_harry, albus_mean_severus),
                      ncol = 2, byrow=T)
  row.names(manual_mean) = c("minerva","albus")
  colnames(manual_mean)=c("HARRY","SEVERUS")

  expect_equal(mean_meth, manual_mean)

})

# Add test for mapToNearest option

# Add test for NA option once developed
# test_that("if there are NA values for some probes mapping to genes, \n
#          but not NA values for others, those probes with NAs are excluded \n
#          from the calculation of the mean if the excludeNA option is set in \n
#          the function call.",{})
