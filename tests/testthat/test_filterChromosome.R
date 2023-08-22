context("[NetSciDataCompanion] Testing logNormalization functions ... ")

test_that("Testing logTPMNormalization format",{
  
  my_friend = NetSciDataCompanion::CreateNetSciDataCompanionObject()
  human_projects <- available_projects()
  proj_info <- subset(
    human_projects,
    project == "SRP009615" & project_type == "data_sources"
  )
  expr <- create_rse(proj_info)
  
  # Test integration with rds
  gene_info <- my_friend$extractSampleAndGeneInfo(expr)$rds_gene_info
  chroms <- c("chr6", "chr7")
  res <- my_friend$filterChromosome(gene_info, chroms)
  expect_true(length(res) > 0)
  
  # Test functionality
  data <- matrix(0, nrow=5, ncol=dim(gene_info)[2])
  data[,1] <- c("chr1", "chr6", "chr5", "chr6", "chr7")
  colnames(data) <- colnames(gene_info)
  expect_equal(my_friend$filterChromosome(data.frame(data), chroms), c(2, 4, 5))
})