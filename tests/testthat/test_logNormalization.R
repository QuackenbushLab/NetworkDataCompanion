context("[NetSciDataCompanion] Testing logNormalization functions ... ")

test_that("Testing logCPMNormalization format",{
  
  my_friend = NetSciDataCompanion::CreateNetSciDataCompanionObject()
  nrow <- 300
  ncol <- 20
  expr <- matrix(rpois(36,5),nrow=nrow, ncol = ncol)
  normalized <- my_friend$logCPMNormalization(expr)
  
  # Test dimensionality
  expect_equal(dim(normalized$counts), dim(expr))
  expect_equal(dim(normalized$CPM), dim(expr))
  expect_equal(dim(normalized$logCPM), dim(expr))
  
  # Test counts
  expect_equal(normalized$counts, expr)
  
  # Test logs
  expect_equal(normalized$logCPM, log(normalized$CPM + 1))
})

test_that("Testing logTPMNormalization format",{
  
  my_friend = NetSciDataCompanion::CreateNetSciDataCompanionObject()
  human_projects <- available_projects()
  proj_info <- subset(
    human_projects,
    project == "SRP009615" & project_type == "data_sources"
  )
  expr <- create_rse(proj_info)
  
  normalized <- my_friend$logTPMNormalization(expr)
  
  # Test dimensionality
  expect_equal(dim(normalized$counts), dim(expr))
  expect_equal(dim(normalized$TPM), dim(expr))
  expect_equal(dim(normalized$logTPM), dim(expr))
  
  # Test logs
  expect_equal(normalized$logTPM, log(normalized$TPM + 1))
})