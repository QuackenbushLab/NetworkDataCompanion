context("[NetSciDataCompanion] Testing mapProbesToGenes function ... ")

test_that("probes are mapped correctly to TSS200",{
  # make a toy manifest
  manifest_line_1 = data.frame("CpG_chrm"="chr42",
                               "CpG_beg"=3141,
                               "CpG_end"=3142,
                               "probe_strand"="+",
                               "probeID"="cg00000001",
                               "genesUniq"="HARRY",
                               "geneNames"="HARRY",
                               "transcriptTypes"=NA,
                               "transcriptIDs"=NA,
                               "distToTSS"=1,
                               "CGI"=NA,
                               "CGIposition"=NA)
  manifest_line_2 = data.frame("CpG_chrm"="chr42",
                               "CpG_beg"=5926,
                               "CpG_end"=5927,
                               "probe_strand"="+",
                               "probeID"="cg00000002",
                               "genesUniq"="RON",
                               "geneNames"="RON",
                               "transcriptTypes"=NA,
                               "transcriptIDs"=NA,
                               "distToTSS"=-150,
                               "CGI"=NA,
                               "CGIposition"=NA)

  # @Jonas: is this right that distToTSS is already calculated
  # based on stranding?
  manifest = rbind.data.frame(manifest_line_1,manifest_line_2)
  write.table(manifest,file="testdata/manifest_test.tsv",sep="\t",row.names=F,quote=F)
  my_friend = NetSciDataCompanion::CreateNetSciDataCompanionObject()
  my_map_calc = my_friend$mapProbesToGenes(probelist = c("cg00000001","cg00000002"),
                             localManifestPath = "testdata/manifest_test.tsv")

  # expected output
  my_map_hand = data.frame("probeID"=c("cg00000001","cg00000002"),
                           "geneName"=c(NA,"RON"),
                           "ensemblID"=NA,
                           "distToTSS"=c(NA,-150))
  expect_equal(my_map_hand$probeID, my_map_calc$probeID)
  expect_equal(my_map_hand$geneName, my_map_calc$geneName)
  expect_equal(my_map_hand$ensemblID, as.logical(my_map_calc$ensemblID))
  expect_equal(my_map_hand$distToTSS, as.double(my_map_calc$distToTSS))

})

# TODO: Implement more flexible slicing
# (i.e. slicing that allows user to find probes in intervals that don't contain the TSS)

test_that("probes are mapped correctly to a custom region, [TSS - 10: TSS + 10]",{
  # make a toy manifest
  manifest_line_1 = data.frame("CpG_chrm"="chr42",
                               "CpG_beg"=3141,
                               "CpG_end"=3142,
                               "probe_strand"="+",
                               "probeID"="cg00000001",
                               "genesUniq"="HARRY",
                               "geneNames"="HARRY",
                               "transcriptTypes"=NA,
                               "transcriptIDs"=NA,
                               "distToTSS"=1,
                               "CGI"=NA,
                               "CGIposition"=NA)
  manifest_line_2 = data.frame("CpG_chrm"="chr42",
                               "CpG_beg"=5926,
                               "CpG_end"=5927,
                               "probe_strand"="+",
                               "probeID"="cg00000002",
                               "genesUniq"="RON",
                               "geneNames"="RON",
                               "transcriptTypes"=NA,
                               "transcriptIDs"=NA,
                               "distToTSS"=-150,
                               "CGI"=NA,
                               "CGIposition"=NA)
  manifest_line_3 = data.frame("CpG_chrm"="chr42",
                               "CpG_beg"=535,
                               "CpG_end"=536,
                               "probe_strand"="+",
                               "probeID"="cg00000003",
                               "genesUniq"="HERMIONE",
                               "geneNames"="HERMIONE",
                               "transcriptTypes"=NA,
                               "transcriptIDs"=NA,
                               "distToTSS"=-9,
                               "CGI"=NA,
                               "CGIposition"=NA)

  # @Jonas: is this right that distToTSS is already calculated
  # based on stranding?
  manifest = rbind.data.frame(manifest_line_1,manifest_line_2,manifest_line_3)
  write.table(manifest,file="testdata/manifest_test.tsv",sep="\t",row.names=F,quote=F)
  my_friend = NetSciDataCompanion::CreateNetSciDataCompanionObject()
  my_map_calc = my_friend$mapProbesToGenes(probelist = c("cg00000001","cg00000002","cg00000003"),
                                           rangeUp=10,
                                           rangeDown=10,
                                           localManifestPath = "testdata/manifest_test.tsv")

  # expected output
  my_map_hand = data.frame("probeID"=c("cg00000001","cg00000002","cg00000003"),
                           "geneName"=c("HARRY",NA,"HERMIONE"),
                           "ensemblID"=NA,
                           "distToTSS"=c(1,NA,-9))
  expect_equal(my_map_hand$probeID, my_map_calc$probeID)
  expect_equal(my_map_hand$geneName, my_map_calc$geneName)
  expect_equal(my_map_hand$ensemblID, as.logical(my_map_calc$ensemblID))
  expect_equal(my_map_hand$distToTSS, as.double(my_map_calc$distToTSS))
})

# test download, gunzip, and load of manifest
test_that("Illumina manifest is correctly downloaded and gunzipped...",
{
  my_friend = NetSciDataCompanion::CreateNetSciDataCompanionObject()
  nsdc_map = my_friend$mapProbesToGenes(probelist="cg05575921", rangeUp = 0, rangeDown = 51000)
  # true location hard-coded from the manifest located at
  # https://zhouserver.research.chop.edu/InfiniumAnnotation/20210615/HM450/HM450.hg38.manifest.gencode.v36.tsv.gz
  #  probeID geneNames                            ensemblID
  #1 cg05575921 AHRR;AHRR ENST00000316418.10;ENST00000510400.5
  #distToTSS
  #1 50676;50693
  true_map = data.frame("probeID"="cg05575921",
                        "geneNames"="AHRR;AHRR",
                        "ensemblID"="ENST00000316418.10;ENST00000510400.5",
                        "distToTSS"="50676;50693")
  expect_equal(nsdc_map,true_map)
})
