context("[NetworkDataCompanion] Testing geneENSGToName functions ... ")

test_that("geneENSGToName functions correctly converts ENSG to Name",{

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
  out = my_friend$geneENSGToName(gene_id)
  expect_equal(out$gene_name,rep(gene_name, times = c(1, 3, 2)))

  ## without version number
  out = my_friend$geneENSGToName(gene_id_no_ver)
  expect_equal(out$gene_name,rep(gene_name, times = c(1, 3, 2)))
  
  ## assert stop if mixed input
  gene_ids_mixed = c(gene_id[1:5],gene_id_no_ver[6])
  expect_error(my_friend$geneENSGToName(gene_ids_mixed),
               regexp="\\[NetworkDataCompanion::getGeneInfo\\] Currently, ensembl IDs must either all have versions or all have no versions. \n Please adjust your input accordingly.")
})
