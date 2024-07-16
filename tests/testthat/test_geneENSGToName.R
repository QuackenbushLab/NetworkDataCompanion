context("[NetZooDataCompanion] Testing geneENSGToName functions ... ")

test_that("geneENSGToName functions correctly converts ENSG to Name",{

  my_friend = NetZooDataCompanion::CreateNetZooDataCompanionObject()

  ## test genes: MIF has two ids, TP53 has one id, and WNT3 has three ids.
  gene_name = c("MIF","TP53","WNT3")
  gene_entrez = c(4282,7157,7473)
  gene_id = c("ENSG00000240972.1","ENSG00000276701.2",
              "ENSG00000141510.16",
              "ENSG00000277626.1", "ENSG00000108379.9", "ENSG00000277641.2")
  gene_id_no_ver = c("ENSG00000240972","ENSG00000276701",
                     "ENSG00000141510",
                     "ENSG00000277626", "ENSG00000108379", "ENSG00000277641")

  ## with version number
  out = my_friend$geneENSGToName(gene_id)
  expect_equal(out$gene_name,rep(gene_name, times = c(2, 1, 3)))

  ## without version number
  out = my_friend$geneENSGToName(gene_id_no_ver)
  expect_equal(out$gene_name,rep(gene_name, times = c(2, 1, 3)))

})
