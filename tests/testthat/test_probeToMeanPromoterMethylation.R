# Test probeToMeanTFMethylation
context("[NetSciDataCompanion] Testing probeToMeanTFMethylation function ... ")

test_that("all probes that map to a gene are extracted from the map",{
  # make a toy manifest with three probes that map to the same gene
  # and one that does not
  manifest_line_1 = data.frame("CpG_chrm"="chr42",
                               "CpG_beg"=3141,
                               "CpG_end"=3142,
                               "probe_strand"="+",
                               "probeID"="cg00000001",
                               "genesUniq"="HARRY",
                               "geneNames"="HARRY",
                               "transcriptTypes"=NA,
                               "transcriptIDs"=NA,
                               "distToTSS"=-150,
                               "CGI"=NA,
                               "CGIposition"=NA)
  manifest_line_2 = data.frame("CpG_chrm"="chr42",
                               "CpG_beg"=3143,
                               "CpG_end"=3144,
                               "probe_strand"="+",
                               "probeID"="cg00000002",
                               "genesUniq"="HARRY",
                               "geneNames"="HARRY",
                               "transcriptTypes"=NA,
                               "transcriptIDs"=NA,
                               "distToTSS"=-125,
                               "CGI"=NA,
                               "CGIposition"=NA)
  manifest_line_3 = data.frame("CpG_chrm"="chr42",
                               "CpG_beg"=3145,
                               "CpG_end"=3146,
                               "probe_strand"="+",
                               "probeID"="cg00000003",
                               "genesUniq"="HARRY",
                               "geneNames"="HARRY",
                               "transcriptTypes"=NA,
                               "transcriptIDs"=NA,
                               "distToTSS"=-120,
                               "CGI"=NA,
                               "CGIposition"=NA)
  manifest_line_4 = data.frame("CpG_chrm"="chr42",
                               "CpG_beg"=3245,
                               "CpG_end"=3246,
                               "probe_strand"="+",
                               "probeID"="cg00000004",
                               "genesUniq"="SEVERUS",
                               "geneNames"="SEVERUS",
                               "transcriptTypes"=NA,
                               "transcriptIDs"=NA,
                               "distToTSS"=100,
                               "CGI"=NA,
                               "CGIposition"=NA)
  manifest = rbind.data.frame(manifest_line_1, manifest_line_2, manifest_line_3,manifest_line_4)

  # make toy beta values
  my_betas = data.frame("probeID"=c("cg00000001",
                                    "cg00000002",
                                    "cg00000003",
                                    "cg00000004"))

  set.seed(42)
  my_betas$minerva = runif(4)
  my_betas$albus = runif(4)

  write.table(manifest,file="testdata/manifest_test.tsv",sep="\t",row.names=F,quote=F)
  my_friend = NetSciDataCompanion::CreateNetSciDataCompanionObject()
  my_map = my_friend$mapProbesToGenes(probelist = c("cg00000001","cg00000002","cg00000003","cg00000004"),
                                      rangeUp = 200,
                                      rangeDown = 200,
                                      localManifestPath = "testdata/manifest_test.tsv")
  names(my_map)[2] = "geneNames" #TODO Add switch case for column names depending on manifest type
  my_friend$probeToMeanPromoterMethylation(methylation_betas = my_betas,
                                     genesOfInterest = c("HARRY","SEVERUS"),
                                     probe_gene_map = my_map)
})

test_that("mean calculation is correct",{

  # make toy beta values
  my_betas = data.frame("probeID"=c("cg00000001",
                                    "cg00000002",
                                    "cg00000003",
                                    "cg00000004"))

  set.seed(42)
  my_betas$minerva = runif(4)
  my_betas$albus = runif(4)

  my_friend = NetSciDataCompanion::CreateNetSciDataCompanionObject()
  my_map = my_friend$mapProbesToGenes(probelist = c("cg00000001","cg00000002","cg00000003","cg00000004"),
                                      rangeUp = 200,
                                      rangeDown = 200,
                                      localManifestPath = "testdata/manifest_test.tsv")
  names(my_map)[2] = "geneNames"
  mean_meth = my_friend$probeToMeanPromoterMethylation(methylation_betas = my_betas,
                                           genesOfInterest = c("HARRY","SEVERUS"),
                                           probe_gene_map = my_map)
  minerva_mean_harry = mean(my_betas$minerva[1:3])
  minerva_mean_severus = my_betas$minerva[4]
  albus_mean_harry = mean(my_betas$albus[1:3])
  albus_mean_severus = my_betas$albus[4]

  manual_mean = matrix(c(minerva_mean_harry, minerva_mean_severus,
                        albus_mean_harry, albus_mean_severus),
                      ncol = 2, byrow=T)
  row.names(manual_mean) = c("minerva","albus")
  colnames(manual_mean)=c("HARRY","SEVERUS")

  expect_equal(mean_meth, manual_mean)
})

