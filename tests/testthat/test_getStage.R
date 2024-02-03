context("[NetSciDataCompanion] Testing getStage function ... ")

test_that("Testing getStage",{
  
  my_friend = NetSciDataCompanion::CreateNetSciDataCompanionObject("./tests/testthat/testdata/clinical_patient.csv")
  barcode1 = "TCGA-A1-1234"
  barcode2 = "TCGA-A1-1236"
  barcodes = c(barcode1, barcode2)
  expect_equal(my_friend$getStage(barcodes), c(3, 1))
  
  # Add non-existing sample to the query
  barcodes <- c(barcodes, "TCGA-A1-0000")
  expect_equal(my_friend$getStage(barcodes), c(3, 1, NA))
})