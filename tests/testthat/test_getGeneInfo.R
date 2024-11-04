context("[NetworkDataCompanion] Testing getGeneInfo functions ... ")

test_that("getGeneInfo functions correctly gets gene info from genecode V26",{

  my_friend = NetworkDataCompanion::CreateNetworkDataCompanionObject()

  ## test genes: MIF has two ids, TP53 has one id, and WNT3 has three ids.
  gene_name = c("MIF","TP53","WNT3")
  gene_entrez = c(4282,7157,7473)
  gene_id = c("ENSG00000240972.1","ENSG00000276701.2",
              "ENSG00000141510.16",
              "ENSG00000277626.1", "ENSG00000108379.9", "ENSG00000277641.2")
  gene_id_no_ver = c("ENSG00000240972","ENSG00000276701",
                     "ENSG00000141510",
                     "ENSG00000277626", "ENSG00000108379", "ENSG00000277641")

  ## input gene name
  gene_info = my_friend$getGeneInfo(gene_name)
  expect_equal(gene_info$gene_name,rep(gene_name, times = c(2, 1, 3)))
  expect_equal(gene_info$gene_entrez,rep(gene_entrez, times = c(2, 1, 3)))
  expect_equal(gene_info$gene_id,gene_id)
  expect_equal(gene_info$gene_id_no_ver,gene_id_no_ver)

  ## input ensemble id
  rm(gene_info)
  gene_info = my_friend$getGeneInfo(gene_id)
  expect_equal(gene_info$gene_name,rep(gene_name, times = c(2, 1, 3)))
  expect_equal(gene_info$gene_entrez,rep(gene_entrez, times = c(2, 1, 3)))
  expect_equal(gene_info$gene_id,gene_id)
  expect_equal(gene_info$gene_id_no_ver,gene_id_no_ver)

  ## input ensemble id without version number
  rm(gene_info)
  gene_info = my_friend$getGeneInfo(gene_id_no_ver)
  expect_equal(gene_info$gene_name,rep(gene_name, times = c(2, 1, 3)))
  expect_equal(gene_info$gene_entrez,rep(gene_entrez, times = c(2, 1, 3)))
  expect_equal(gene_info$gene_id,gene_id)
  expect_equal(gene_info$gene_id_no_ver,gene_id_no_ver)

  ## input entrez id
  rm(gene_info)
  gene_info = my_friend$getGeneInfo(gene_entrez)
  expect_equal(gene_info$gene_name,rep(gene_name, times = c(2, 1, 3)))
  expect_equal(gene_info$gene_entrez,rep(gene_entrez, times = c(2, 1, 3)))
  expect_equal(gene_info$gene_id,gene_id)
  expect_equal(gene_info$gene_id_no_ver,gene_id_no_ver)
})
