context("[NetworkDataCompanion] Testing geneNameToENSG functions ... ")

test_that("geneNameToENSG functions correctly converts Name to ENSG",{

  my_friend = NetworkDataCompanion::CreateNetworkDataCompanionObject()

  ## test genes: MIF has two ids, TP53 has one id, and WNT3 has three ids.
  gene_name = c("TP53","WNT3","MIF")
  gene_entrez = c(7157,7473,4282)
  gene_id = c("ENSG00000141510.16",
              "ENSG00000277626.1", "ENSG00000108379.9", "ENSG00000277641.2",
              "ENSG00000240972.1","ENSG00000276701.2")
  gene_id_no_ver = c("ENSG00000141510",
                     "ENSG00000277626", "ENSG00000108379", "ENSG00000277641",
                     "ENSG00000240972","ENSG00000276701")

  ## with version number
  out = my_friend$geneNameToENSG(gene_name,include_no_version = FALSE)
  expect_equal(out$gene_id,gene_id)

  ## include without version number
  out = my_friend$geneNameToENSG(gene_name,include_no_version = TRUE)
  expect_equal(out$gene_id_no_ver,gene_id_no_ver)

})
