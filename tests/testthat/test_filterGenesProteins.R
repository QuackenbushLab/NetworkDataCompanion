context("[NetSciDataCompanion] Testing filterGenesProteins function ... ")

test_that("Testing filterGenesProteins",{
  
  my_friend = NetSciDataCompanion::CreateNetSciDataCompanionObject()
  human_projects <- available_projects()
  proj_info <- subset(
    human_projects,
    project == "SRP009615" & project_type == "data_sources"
  )
  expr <- create_rse(proj_info)
  
  # Test integration with rds
  gene_info <- my_friend$extractSampleAndGeneInfo(expr)$rds_gene_info
  res <- my_friend$filterGenesProteins(gene_info)
  expect_true(length(res) > 0)
  
  # Test functionality
  expect_equal(unique(gene_info[res,]$gene_type), 'protein_coding')
})