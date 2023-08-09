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
                               "probeID"="cg00000005",
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
  my_map_200_0 = my_friend$mapProbesToGenes(probelist = NULL, # NULL probelist maps everything
                                      rangeUp = 200,
                                      rangeDown = 0,
                                      localManifestPath = "testdata/manifest_test.tsv")

  my_map_200_101 = my_friend$mapProbesToGenes(probelist = NULL, # NULL probelist maps everything
                                      rangeUp = 200,
                                      rangeDown = 101,
                                      localManifestPath = "testdata/manifest_test.tsv")

  # test that all probes mapped to the right genes

  names(my_map_200_0)[2] = "geneNames"
  meanMeth = my_friend$probeToMeanPromoterMethylation(methylation_betas = my_betas,
                                     genesOfInterest = c("HARRY","SEVERUS"),
                                     probe_gene_map = my_map_200_0)

})

