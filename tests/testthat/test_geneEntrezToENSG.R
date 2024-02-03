context("[NetSciDataCompanion] Testing geneEntrezToENSG functions ... ")

test_that("geneEntrezToENSG functions correctly converts Entrez to ENSG",{

  my_friend = NetSciDataCompanion::CreateNetSciDataCompanionObject()

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
  out = my_friend$geneEntrezToENSG(gene_entrez,version = TRUE)
  expect_equal(out$gene_id,gene_id)

  ## without version number
  out = my_friend$geneEntrezToENSG(gene_entrez,version = FALSE)
  expect_equal(out$gene_id,gene_id_no_ver)

})
