context("[NetworkDataCompanion] Testing findDuplicates function ... ")

test_that("Testing findDuplicates",{

  my_friend = NetworkDataCompanion::CreateNetworkDataCompanionObject()
  barcode1 = "TCGA-A1-1234-01A"
  barcode2 = "TCGA-A1-1234-11A"
  barcodes = c(barcode1, barcode2)
  expect_equal(length(my_friend$findDuplicates(barcodes)), 0)
  barcode3 = "TCGA-A1-1234-01A"
  barcode4 = "TCGA-A1-1234-21A"
  barcode5 = "TCGA-A1-1234-01A"
  barcode6 = "TCGA-A1-1234-21A"
  barcodes = c(barcode1, barcode2, barcode3, barcode4, barcode5, barcode6)
  print(my_friend$findDuplicates(barcodes))
  expect_equal(my_friend$findDuplicates(barcodes), c(rep("TCGA-A1-1234-01A", 2), "TCGA-A1-1234-21A"))
})
