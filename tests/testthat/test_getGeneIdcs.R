context("[NetworkDataCompanion] Testing getGeneIdcs function ... ")

test_that("Testing getGeneIdcs",{

  my_friend = NetworkDataCompanion::CreateNetworkDataCompanionObject()
  human_projects <- available_projects()
  proj_info <- subset(
    human_projects,
    project == "SRP009615" & project_type == "data_sources"
  )
  expr <- create_rse(proj_info)

  # Test integration with rds
  gene_info <- my_friend$extractSampleAndGeneInfo(expr)$rds_gene_info
  genes <- c("TP53")
  res <- my_friend$getGeneIdcs(genes, gene_info)
  expect_true(length(res) > 0)

  # Test functionality
  data <- matrix(0, nrow=5, ncol=dim(gene_info)[2])
  data[,12] <- c("A", "B", "C", "B", "D")
  colnames(data) <- colnames(gene_info)
  expect_equal(my_friend$getGeneIdcs(c("A", "B"), data.frame(data)), c(1, 2))
})
