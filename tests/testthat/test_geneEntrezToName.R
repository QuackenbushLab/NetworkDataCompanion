context("[NetworkDataCompanion] Testing geneEntrezToName functions ... ")

test_that("geneEntrezToName functions correctly converts Entrez to Name",{

  my_friend = NetworkDataCompanion::CreateNetworkDataCompanionObject()

  ## test genes:
  ## MIF has two ids, TP53 has one id, and WNT3 has three ids.
  gene_name = c("MIF","TP53","WNT3")
  gene_entrez = c(4282,7157,7473)

  out = my_friend$geneEntrezToName(gene_entrez)
  expect_equal(out$gene_name,gene_name)

})
