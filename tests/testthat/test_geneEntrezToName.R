context("[NetworkDataCompanion] Testing geneEntrezToName functions ... ")

test_that("geneEntrezToName functions correctly converts Entrez to Name",{

  my_friend = NetworkDataCompanion::CreateNetworkDataCompanionObject()

  ## test genes:
  gene_name = c("TP53","WNT3","MIF")
  gene_entrez = c(7157,7473,4282)

  out = my_friend$geneEntrezToName(gene_entrez)
  expect_equal(out$gene_name,gene_name)

})


