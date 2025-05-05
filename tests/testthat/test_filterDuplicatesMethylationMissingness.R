context("[NetworkDataCompanion] Testing filterDuplicatesMethylationMissingness function ... ")

test_that("Testing filterDuplicatesMethylationMissingness",{
  # samples loaded from quickstart, Harvard dataverse archive
  # data_download/methylation/tcga_coad_cms1_methylations.txt
  id_map = read.csv("testdata/duplicate_meth_uuids.csv", row.names=1)
  # fake methylation data with various missingness
  beta_vals = matrix(nrow=3,ncol=6)
  beta_vals[,3] = rep(1,3) # this should be picked
  beta_vals[,4] = c(NA,rep(1,2))
  beta_vals[,5] = rep(1,3) 
  beta_df = data.frame(beta_vals)
  beta_df[,1] = c("cg00000001",
                  "cg00000002",
                  "cg00000003")
  names(beta_df) =c("probeID",id_map$file_id)
  
  my_friend = NetworkDataCompanion::CreateNetworkDataCompanionObject()
  keep_uuids = my_friend$filterDuplicatesMethylationMissingness(beta_df)
  expect_equal(keep_uuids, c(id_map$file_id[2],id_map$file_id[4]))
})
