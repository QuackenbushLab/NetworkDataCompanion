context("[NetSciDataCompanion] Testing getSex function ... ")

test_that("Testing getSex",{
  
  my_friend = NetSciDataCompanion::CreateNetSciDataCompanionObject("./tests/testthat/testdata/clinical_patient.csv")
  barcode1 = "TCGA-A1-1234"
  barcode2 = "TCGA-A1-1236"
  barcodes = c(barcode1, barcode2)
  expect_equal(my_friend$getSex(barcodes), c("M", "F"))
  
  # Add non-existing sample to the query
  barcodes <- c(barcodes, "TCGA-A1-0000")
  expect_equal(my_friend$getSex(barcodes), c("M", "F", NA))
})