context("[NetSciDataCompanion] Testing mapProbesToGenes function ... ")

# begin by making a toy manifest, outside of the test_that functions
# the manifest includes the following special cases:
# 1. A probe that maps to two different genes within 200BP downstream (cg00000004)
# 2. A probe that maps to two different isoforms of the same gene within 200BP downstream
# and a different gene within 200BP upstream (cg00000005)
manifest_line_1 = data.frame("CpG_chrm"="chr42",
                             "CpG_beg"=3141,
                             "CpG_end"=3142,
                             "probe_strand"="+",
                             "probeID"="cg00000001",
                             "genesUniq"="HARRY",
                             "geneNames"="HARRY",
                             "transcriptTypes"=NA,
                             "transcriptIDs"=NA,
                             "distToTSS"=-124,
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
                             "distToTSS"=-122,
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
                             "genesUniq"="HARRY;SEVERUS",
                             "geneNames"="HARRY;SEVERUS",
                             "transcriptTypes"=NA,
                             "transcriptIDs"=NA,
                             "distToTSS"="-20;-200",
                             "CGI"=NA,
                             "CGIposition"=NA)
manifest_line_5 = data.frame("CpG_chrm"="chr42",
                             "CpG_beg"=3440,
                             "CpG_end"=3441,
                             "probe_strand"="+",
                             "probeID"="cg00000005",
                             "genesUniq"="HARRY;SEVERUS;SEVERUS",
                             "geneNames"="HARRY;SEVERUS;SEVERUS",
                             "transcriptTypes"=NA,
                             "transcriptIDs"=NA,
                             "distToTSS"="195;-5;-25",
                             "CGI"=NA,
                             "CGIposition"=NA)
manifest = rbind.data.frame(manifest_line_1, manifest_line_2, manifest_line_3,manifest_line_4,manifest_line_5)
write.table(manifest,file="testdata/manifest_test.tsv",sep="\t",row.names=F,quote=F)

test_that("probes are mapped correctly to TSS200 if no other bound is specified",{
  my_friend = NetSciDataCompanion::CreateNetSciDataCompanionObject()
  my_map_calc = my_friend$mapProbesToGenes(probelist = c("cg00000001",
                                                         "cg00000002",
                                                         "cg00000003",
                                                         "cg00000004",
                                                         "cg00000005"),
                             localManifestPath = "testdata/manifest_test.tsv")

  # expected output
  my_map_hand = data.frame("probeID"=c("cg00000001",
                                       "cg00000002",
                                       "cg00000003",
                                       "cg00000004",
                                       "cg00000005"),
                           "geneName"=c("HARRY",
                                        "HARRY",
                                        "HARRY",
                                        "HARRY", # will not get SEVERUS until TSS201
                                        "SEVERUS;SEVERUS"),
                           "ensemblID"=NA,
                           "distToTSS"=c("-124",
                                         "-122",
                                         "-120",
                                         "-20",
                                         "-5;-25")) # TSS for HARRY is 3265; TSS for SEVERUS.1 is 3445; TSS for SEVERUS.2 is 3465

  expect_equal(my_map_hand$probeID, my_map_calc$probeID)
  expect_equal(my_map_hand$geneName, my_map_calc$geneName)
  expect_equal(my_map_hand$ensemblID, as.logical(my_map_calc$ensemblID))
  expect_equal(my_map_hand$distToTSS, my_map_calc$distToTSS)

})

test_that("probes are mapped correctly to a custom region, [TSS - 201: TSS + 200]",{
  my_friend = NetSciDataCompanion::CreateNetSciDataCompanionObject()
  my_map_calc = my_friend$mapProbesToGenes(probelist = c("cg00000001",
                                                         "cg00000002",
                                                         "cg00000003",
                                                         "cg00000004",
                                                         "cg00000005"),
                                           rangeUp = 201,
                                           rangeDown = 200,
                                           localManifestPath = "testdata/manifest_test.tsv")
  # expected output
  my_map_hand = data.frame("probeID"=c("cg00000001",
                                       "cg00000002",
                                       "cg00000003",
                                       "cg00000004",
                                       "cg00000005"),
                           "geneName"=c("HARRY",
                                        "HARRY",
                                        "HARRY",
                                        "HARRY;SEVERUS", # will not get SEVERUS until TSS201
                                        "HARRY;SEVERUS;SEVERUS"), # gets HARRY upstream
                           "ensemblID"=NA,
                           "distToTSS"=c("-124",
                                         "-122",
                                         "-120",
                                         "-20;-200",
                                         "195;-5;-25")) # TSS for HARRY is 3265; TSS for SEVERUS.1 is 3445; TSS for SEVERUS.2 is 3465
  expect_equal(my_map_hand$probeID, my_map_calc$probeID)
  expect_equal(my_map_hand$geneName, my_map_calc$geneName)
  expect_equal(my_map_hand$ensemblID, as.logical(my_map_calc$ensemblID))
  expect_equal(my_map_hand$distToTSS, my_map_calc$distToTSS)
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
